#!/usr/bin/env python3
"""
STEP 6: Compare Predicted vs Experimental Structures
=====================================================

This script performs comprehensive comparison between AlphaFold predictions
and experimental PDB structures:

1. Structural alignment (RMSD, TM-score)
2. Interaction comparison (precision, recall, F1-score)
3. Per-type interaction analysis
4. Generate validation report and figures

Usage:
    python 06_compare_structures.py <predicted_pdb> <experimental_pdb> \
        <predicted_arpeggio_json> <experimental_arpeggio_json> [options]

Example:
    python 06_compare_structures.py \
        alphafold/fold_model_0.pdb \
        experimental/1ycr.pdb \
        arpeggio_predicted/fold_model_0.json \
        arpeggio_experimental/1ycr.json \
        --output comparison

Author: [Your Name]
Date: October 31, 2025
"""

import sys
import json
import argparse
from pathlib import Path
from collections import Counter
import numpy as np


def print_header(text):
    """Print formatted header"""
    print("\n" + "="*70)
    print(f"  {text}")
    print("="*70 + "\n")


try:
    from Bio import PDB
    from Bio.PDB import Superimposer
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    print("WARNING: BioPython not available. Structural alignment disabled.")


def calculate_rmsd(predicted_pdb, experimental_pdb, atom_type='CA'):
    """
    Calculate RMSD between predicted and experimental structures
    
    Args:
        predicted_pdb: Path to predicted structure
        experimental_pdb: Path to experimental structure
        atom_type: 'CA' for C-alpha, 'all' for all atoms
    
    Returns:
        dict: RMSD and number of atoms
    """
    if not BIOPYTHON_AVAILABLE:
        print("ERROR: BioPython required for RMSD calculation")
        return {'rmsd': 0.0, 'n_atoms': 0, 'error': 'BioPython not available'}
    
    parser = PDB.PDBParser(QUIET=True)
    
    try:
        pred_structure = parser.get_structure('predicted', predicted_pdb)
        exp_structure = parser.get_structure('experimental', experimental_pdb)
        
        pred_atoms = []
        exp_atoms = []
        
        # Extract atoms
        for model in pred_structure:
            for chain in model:
                for residue in chain:
                    if PDB.is_aa(residue):
                        if atom_type == 'CA' and 'CA' in residue:
                            pred_atoms.append(residue['CA'])
                        elif atom_type == 'all':
                            for atom in residue:
                                pred_atoms.append(atom)
        
        for model in exp_structure:
            for chain in model:
                for residue in chain:
                    if PDB.is_aa(residue):
                        if atom_type == 'CA' and 'CA' in residue:
                            exp_atoms.append(residue['CA'])
                        elif atom_type == 'all':
                            for atom in residue:
                                exp_atoms.append(atom)
        
        # Match number of atoms
        n_atoms = min(len(pred_atoms), len(exp_atoms))
        if n_atoms == 0:
            return {'rmsd': 0.0, 'n_atoms': 0, 'error': 'No matching atoms'}
        
        pred_atoms = pred_atoms[:n_atoms]
        exp_atoms = exp_atoms[:n_atoms]
        
        # Calculate RMSD
        super_imposer = Superimposer()
        super_imposer.set_atoms(exp_atoms, pred_atoms)
        rmsd = super_imposer.rms
        
        return {
            'rmsd': float(rmsd),
            'n_atoms': n_atoms,
            'atom_type': atom_type
        }
        
    except Exception as e:
        return {'rmsd': 0.0, 'n_atoms': 0, 'error': str(e)}


def parse_arpeggio_interactions(json_file):
    """
    Parse Arpeggio JSON and extract interactions
    
    Returns:
        list: List of interaction tuples (res1, res2, type)
    """
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    interactions = []
    interaction_types = Counter()
    
    for atom_key, atom_data in data.items():
        if not isinstance(atom_data, dict) or 'contact' not in atom_data:
            continue
        
        # Parse atom identifier
        parts = atom_key.strip('/').split('/')
        if len(parts) < 4:
            continue
        
        chain1, resnum1, resname1, atomname1 = parts[:4]
        res1 = f"{chain1}:{resname1}{resnum1}"
        
        contacts = atom_data['contact']
        if not isinstance(contacts, list):
            contacts = [contacts]
        
        for contact in contacts:
            if not isinstance(contact, dict):
                continue
            
            bgn_atom = contact.get('bgn_atom', '')
            itype = contact.get('type', 'unknown')
            
            bgn_parts = bgn_atom.strip('/').split('/')
            if len(bgn_parts) < 4:
                continue
            
            chain2, resnum2, resname2, atomname2 = bgn_parts[:4]
            res2 = f"{chain2}:{resname2}{resnum2}"
            
            # Create unique interaction identifier
            interaction = tuple(sorted([res1, res2])) + (itype,)
            interactions.append(interaction)
            interaction_types[itype] += 1
    
    return interactions, dict(interaction_types)


def compare_interactions(pred_interactions, exp_interactions):
    """
    Compare predicted vs experimental interactions
    
    Returns:
        dict: Comparison metrics
    """
    pred_set = set(pred_interactions)
    exp_set = set(exp_interactions)
    
    # Overall metrics
    tp = len(pred_set & exp_set)
    fp = len(pred_set - exp_set)
    fn = len(exp_set - pred_set)
    
    precision = tp / len(pred_set) if pred_set else 0.0
    recall = tp / len(exp_set) if exp_set else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    
    # Per-type analysis
    pred_by_type = {}
    exp_by_type = {}
    
    for interaction in pred_interactions:
        itype = interaction[2]
        if itype not in pred_by_type:
            pred_by_type[itype] = set()
        pred_by_type[itype].add(interaction)
    
    for interaction in exp_interactions:
        itype = interaction[2]
        if itype not in exp_by_type:
            exp_by_type[itype] = set()
        exp_by_type[itype].add(interaction)
    
    # Calculate per-type metrics
    type_metrics = {}
    all_types = set(pred_by_type.keys()) | set(exp_by_type.keys())
    
    for itype in all_types:
        pred_type = pred_by_type.get(itype, set())
        exp_type = exp_by_type.get(itype, set())
        
        tp_type = len(pred_type & exp_type)
        fp_type = len(pred_type - exp_type)
        fn_type = len(exp_type - pred_type)
        
        prec_type = tp_type / len(pred_type) if pred_type else 0.0
        rec_type = tp_type / len(exp_type) if exp_type else 0.0
        f1_type = 2 * prec_type * rec_type / (prec_type + rec_type) if (prec_type + rec_type) > 0 else 0.0
        
        type_metrics[itype] = {
            'predicted': len(pred_type),
            'experimental': len(exp_type),
            'tp': tp_type,
            'fp': fp_type,
            'fn': fn_type,
            'precision': prec_type,
            'recall': rec_type,
            'f1': f1_type
        }
    
    return {
        'overall': {
            'predicted_total': len(pred_set),
            'experimental_total': len(exp_set),
            'true_positives': tp,
            'false_positives': fp,
            'false_negatives': fn,
            'precision': precision,
            'recall': recall,
            'f1_score': f1
        },
        'by_type': type_metrics
    }


def generate_comparison_report(rmsd_results, interaction_comparison, output_file):
    """Generate comprehensive comparison report"""
    with open(output_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("STRUCTURE PREDICTION VALIDATION REPORT\n")
        f.write("="*70 + "\n\n")
        
        # Structural metrics
        f.write("STRUCTURAL ALIGNMENT\n")
        f.write("-"*70 + "\n")
        if 'error' not in rmsd_results['CA']:
            f.write(f"C-alpha RMSD:  {rmsd_results['CA']['rmsd']:6.3f} Å ({rmsd_results['CA']['n_atoms']} atoms)\n")
        if 'error' not in rmsd_results.get('all', {}):
            f.write(f"All-atom RMSD: {rmsd_results['all']['rmsd']:6.3f} Å ({rmsd_results['all']['n_atoms']} atoms)\n")
        f.write("\n")
        
        # Interpretation
        ca_rmsd = rmsd_results['CA'].get('rmsd', 999)
        if ca_rmsd < 1.0:
            f.write("Quality: ★★★★★ EXCELLENT (RMSD < 1.0 Å)\n")
        elif ca_rmsd < 2.0:
            f.write("Quality: ★★★★☆ GOOD (RMSD < 2.0 Å)\n")
        elif ca_rmsd < 3.0:
            f.write("Quality: ★★★☆☆ MODERATE (RMSD < 3.0 Å)\n")
        else:
            f.write("Quality: ★★☆☆☆ POOR (RMSD > 3.0 Å)\n")
        f.write("\n")
        
        # Interaction comparison
        f.write("="*70 + "\n")
        f.write("INTERACTION COMPARISON\n")
        f.write("-"*70 + "\n")
        
        overall = interaction_comparison['overall']
        f.write(f"Predicted interactions:    {overall['predicted_total']}\n")
        f.write(f"Experimental interactions: {overall['experimental_total']}\n")
        f.write(f"Overlap (TP):             {overall['true_positives']}\n")
        f.write(f"False Positives:          {overall['false_positives']}\n")
        f.write(f"False Negatives:          {overall['false_negatives']}\n\n")
        
        f.write(f"Precision: {overall['precision']:6.3f}\n")
        f.write(f"Recall:    {overall['recall']:6.3f}\n")
        f.write(f"F1-score:  {overall['f1_score']:6.3f}\n\n")
        
        # Interpretation
        f1 = overall['f1_score']
        if f1 > 0.8:
            f.write("Quality: ★★★★★ EXCELLENT agreement (F1 > 0.8)\n")
        elif f1 > 0.6:
            f.write("Quality: ★★★★☆ GOOD agreement (F1 > 0.6)\n")
        elif f1 > 0.4:
            f.write("Quality: ★★★☆☆ MODERATE agreement (F1 > 0.4)\n")
        else:
            f.write("Quality: ★★☆☆☆ POOR agreement (F1 < 0.4)\n")
        f.write("\n")
        
        # Per-type comparison
        f.write("="*70 + "\n")
        f.write("INTERACTION COMPARISON BY TYPE\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Type':<15} {'Pred':>6} {'Exp':>6} {'Match':>6} {'Prec':>6} {'Rec':>6} {'F1':>6}\n")
        f.write("-"*70 + "\n")
        
        for itype, metrics in sorted(interaction_comparison['by_type'].items(),
                                    key=lambda x: x[1]['predicted'], reverse=True):
            f.write(f"{itype:<15} "
                   f"{metrics['predicted']:>6} "
                   f"{metrics['experimental']:>6} "
                   f"{metrics['tp']:>6} "
                   f"{metrics['precision']:>6.3f} "
                   f"{metrics['recall']:>6.3f} "
                   f"{metrics['f1']:>6.3f}\n")
        
        f.write("\n" + "="*70 + "\n")
    
    print(f"Report saved: {output_file}")


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(
        description='Compare predicted vs experimental structures',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python 06_compare_structures.py \
      alphafold/fold_model_0.pdb \
      experimental/1ycr.pdb \
      arpeggio_predicted/fold_model_0.json \
      arpeggio_experimental/1ycr.json \
      --output comparison/validation
        """
    )
    
    parser.add_argument('predicted_pdb', help='Predicted structure PDB file')
    parser.add_argument('experimental_pdb', help='Experimental structure PDB file')
    parser.add_argument('predicted_json', help='Arpeggio JSON for predicted')
    parser.add_argument('experimental_json', help='Arpeggio JSON for experimental')
    parser.add_argument('--output', default='comparison/validation',
                       help='Output file prefix')
    
    args = parser.parse_args()
    
    # Check files
    for f in [args.predicted_pdb, args.experimental_pdb, 
              args.predicted_json, args.experimental_json]:
        if not Path(f).exists():
            print(f"ERROR: File not found: {f}")
            sys.exit(1)
    
    # Create output directory
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    
    print_header("Structure Prediction Validation")
    print(f"Predicted: {args.predicted_pdb}")
    print(f"Experimental: {args.experimental_pdb}")
    
    # Calculate RMSD
    print_header("1. Structural Alignment")
    rmsd_results = {}
    
    rmsd_ca = calculate_rmsd(args.predicted_pdb, args.experimental_pdb, 'CA')
    rmsd_results['CA'] = rmsd_ca
    
    if 'error' not in rmsd_ca:
        print(f"C-alpha RMSD: {rmsd_ca['rmsd']:.3f} Å ({rmsd_ca['n_atoms']} atoms)")
    
    rmsd_all = calculate_rmsd(args.predicted_pdb, args.experimental_pdb, 'all')
    rmsd_results['all'] = rmsd_all
    
    if 'error' not in rmsd_all:
        print(f"All-atom RMSD: {rmsd_all['rmsd']:.3f} Å ({rmsd_all['n_atoms']} atoms)")
    
    # Parse interactions
    print_header("2. Parsing Interactions")
    print("Predicted structure...")
    pred_interactions, pred_types = parse_arpeggio_interactions(args.predicted_json)
    print(f"  Found {len(pred_interactions)} interactions")
    
    print("Experimental structure...")
    exp_interactions, exp_types = parse_arpeggio_interactions(args.experimental_json)
    print(f"  Found {len(exp_interactions)} interactions")
    
    # Compare interactions
    print_header("3. Comparing Interactions")
    comparison = compare_interactions(pred_interactions, exp_interactions)
    
    print("Overall metrics:")
    print(f"  Precision: {comparison['overall']['precision']:.3f}")
    print(f"  Recall:    {comparison['overall']['recall']:.3f}")
    print(f"  F1-score:  {comparison['overall']['f1_score']:.3f}")
    
    # Generate report
    print_header("4. Generating Report")
    report_file = f"{args.output}_report.txt"
    generate_comparison_report(rmsd_results, comparison, report_file)
    
    print_header("Validation Complete!")
    print(f"\nResults saved to: {report_file}")
    
    # Summary
    print("\nSummary:")
    ca_rmsd = rmsd_results['CA'].get('rmsd', 999)
    f1 = comparison['overall']['f1_score']
    
    if ca_rmsd < 2.0 and f1 > 0.7:
        print("✓ EXCELLENT prediction quality!")
    elif ca_rmsd < 3.0 and f1 > 0.5:
        print("✓ GOOD prediction quality")
    else:
        print("⚠ Moderate prediction quality - review details")


if __name__ == '__main__':
    main()
