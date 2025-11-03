#!/usr/bin/env python3
"""
STEP 4: Analyze Arpeggio Interaction Results
=============================================

This script parses Arpeggio JSON output and provides detailed analysis:
- Interaction counts by type
- Hot spot residue identification
- Interface characteristics
- Per-residue breakdown
- Export to CSV for further analysis

Usage:
    python 04_analyze_interactions.py <arpeggio_json> [options]

Example:
    python 04_analyze_interactions.py arpeggio_results/fold_model_0.json \
        --output analysis --csv --hot-spots

Author: [Your Name]
Date: October 31, 2025
"""

import json
import sys
from pathlib import Path
from collections import defaultdict, Counter
import argparse


def print_header(text):
    """Print formatted header"""
    print("\n" + "="*70)
    print(f"  {text}")
    print("="*70 + "\n")


def parse_arpeggio_json(json_file):
    """
    Parse Arpeggio JSON file
    
    Returns:
        dict: Parsed interaction data
    """
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data


def extract_interactions(arpeggio_data):
    """
    Extract interactions from Arpeggio data
    
    Returns:
        list: List of interaction dictionaries
    """
    interactions = []
    
    for atom_key, atom_data in arpeggio_data.items():
        if not isinstance(atom_data, dict):
            continue
        
        # Parse atom identifier
        # Format: /chain/residue_number/residue_name/atom_name
        parts = atom_key.strip('/').split('/')
        if len(parts) < 4:
            continue
        
        chain1, resnum1, resname1, atomname1 = parts[:4]
        
        # Get contacts
        if 'contact' not in atom_data:
            continue
        
        contacts = atom_data['contact']
        if not isinstance(contacts, list):
            contacts = [contacts]
        
        for contact in contacts:
            if not isinstance(contact, dict):
                continue
            
            # Extract contact information
            bgn_atom = contact.get('bgn_atom', '')
            interaction_type = contact.get('type', 'unknown')
            distance = contact.get('distance', 0.0)
            
            # Parse contact atom identifier
            bgn_parts = bgn_atom.strip('/').split('/')
            if len(bgn_parts) < 4:
                continue
            
            chain2, resnum2, resname2, atomname2 = bgn_parts[:4]
            
            # Store interaction
            interactions.append({
                'chain1': chain1,
                'resnum1': resnum1,
                'resname1': resname1,
                'atom1': atomname1,
                'chain2': chain2,
                'resnum2': resnum2,
                'resname2': resname2,
                'atom2': atomname2,
                'type': interaction_type,
                'distance': float(distance) if distance else 0.0
            })
    
    return interactions


def analyze_interactions(interactions):
    """
    Analyze interaction patterns
    
    Returns:
        dict: Analysis results
    """
    # Count by type
    type_counts = Counter(inter['type'] for inter in interactions)
    
    # Count by residue pair
    residue_pair_counts = defaultdict(list)
    for inter in interactions:
        res1 = f"{inter['chain1']}:{inter['resname1']}{inter['resnum1']}"
        res2 = f"{inter['chain2']}:{inter['resname2']}{inter['resnum2']}"
        key = tuple(sorted([res1, res2]))
        residue_pair_counts[key].append(inter['type'])
    
    # Count interactions per residue
    residue_counts = defaultdict(int)
    for inter in interactions:
        res1 = f"{inter['chain1']}:{inter['resname1']}{inter['resnum1']}"
        res2 = f"{inter['chain2']}:{inter['resname2']}{inter['resnum2']}"
        residue_counts[res1] += 1
        residue_counts[res2] += 1
    
    # Distance statistics
    distances = [inter['distance'] for inter in interactions if inter['distance'] > 0]
    
    return {
        'total_interactions': len(interactions),
        'type_counts': dict(type_counts),
        'residue_pair_counts': dict(residue_pair_counts),
        'residue_counts': dict(residue_counts),
        'unique_residue_pairs': len(residue_pair_counts),
        'distance_stats': {
            'mean': sum(distances) / len(distances) if distances else 0,
            'min': min(distances) if distances else 0,
            'max': max(distances) if distances else 0
        }
    }


def identify_hot_spots(residue_counts, min_interactions=3):
    """
    Identify hot spot residues (high interaction count)
    
    Args:
        residue_counts: Dictionary of residue interaction counts
        min_interactions: Minimum interactions to be considered hot spot
    
    Returns:
        list: Sorted list of (residue, count) tuples
    """
    hot_spots = [
        (res, count) for res, count in residue_counts.items()
        if count >= min_interactions
    ]
    return sorted(hot_spots, key=lambda x: x[1], reverse=True)


def export_to_csv(interactions, output_file):
    """Export interactions to CSV file"""
    import csv
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Chain1', 'Residue1', 'Atom1',
            'Chain2', 'Residue2', 'Atom2',
            'Type', 'Distance'
        ])
        
        for inter in interactions:
            writer.writerow([
                inter['chain1'],
                f"{inter['resname1']}{inter['resnum1']}",
                inter['atom1'],
                inter['chain2'],
                f"{inter['resname2']}{inter['resnum2']}",
                inter['atom2'],
                inter['type'],
                f"{inter['distance']:.2f}"
            ])
    
    print(f"Exported {len(interactions)} interactions to {output_file}")


def generate_report(analysis, output_file):
    """Generate text report"""
    with open(output_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("ARPEGGIO INTERACTION ANALYSIS REPORT\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Total interactions: {analysis['total_interactions']}\n")
        f.write(f"Unique residue pairs: {analysis['unique_residue_pairs']}\n\n")
        
        f.write("INTERACTIONS BY TYPE\n")
        f.write("-"*70 + "\n")
        for itype, count in sorted(analysis['type_counts'].items(), 
                                   key=lambda x: x[1], reverse=True):
            f.write(f"  {itype:20s} : {count:5d}\n")
        
        f.write("\nDISTANCE STATISTICS\n")
        f.write("-"*70 + "\n")
        f.write(f"  Mean:    {analysis['distance_stats']['mean']:.2f} Å\n")
        f.write(f"  Min:     {analysis['distance_stats']['min']:.2f} Å\n")
        f.write(f"  Max:     {analysis['distance_stats']['max']:.2f} Å\n")
        
        f.write("\n" + "="*70 + "\n")
    
    print(f"Report saved to {output_file}")


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(
        description='Analyze Arpeggio interaction results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis
  python 04_analyze_interactions.py arpeggio_results/structure.json
  
  # Full analysis with exports
  python 04_analyze_interactions.py arpeggio_results/structure.json \
      --output analysis --csv --report --hot-spots
        """
    )
    
    parser.add_argument('json_file', help='Arpeggio JSON output file')
    parser.add_argument('--output', default='analysis',
                       help='Output file prefix (default: analysis)')
    parser.add_argument('--csv', action='store_true',
                       help='Export to CSV')
    parser.add_argument('--report', action='store_true',
                       help='Generate text report')
    parser.add_argument('--hot-spots', action='store_true',
                       help='Identify hot spot residues')
    parser.add_argument('--min-interactions', type=int, default=3,
                       help='Minimum interactions for hot spot (default: 3)')
    
    args = parser.parse_args()
    
    # Check if file exists
    if not Path(args.json_file).exists():
        print(f"ERROR: File not found: {args.json_file}")
        sys.exit(1)
    
    print_header("Arpeggio Interaction Analysis")
    print(f"Input file: {args.json_file}")
    print(f"Output prefix: {args.output}")
    
    # Parse JSON
    print("\nParsing Arpeggio output...")
    arpeggio_data = parse_arpeggio_json(args.json_file)
    
    # Extract interactions
    print("Extracting interactions...")
    interactions = extract_interactions(arpeggio_data)
    print(f"Found {len(interactions)} interactions")
    
    # Analyze
    print("\nAnalyzing interaction patterns...")
    analysis = analyze_interactions(interactions)
    
    # Display summary
    print_header("Summary Statistics")
    print(f"Total interactions: {analysis['total_interactions']}")
    print(f"Unique residue pairs: {analysis['unique_residue_pairs']}")
    print(f"\nInteractions by type:")
    for itype, count in sorted(analysis['type_counts'].items(), 
                               key=lambda x: x[1], reverse=True):
        print(f"  {itype:20s} : {count:5d}")
    
    print(f"\nDistance statistics:")
    print(f"  Mean: {analysis['distance_stats']['mean']:.2f} Å")
    print(f"  Range: {analysis['distance_stats']['min']:.2f} - {analysis['distance_stats']['max']:.2f} Å")
    
    # Hot spots
    if args.hot_spots:
        print_header("Hot Spot Residues")
        hot_spots = identify_hot_spots(
            analysis['residue_counts'],
            args.min_interactions
        )
        
        if hot_spots:
            print(f"Residues with ≥{args.min_interactions} interactions:")
            for residue, count in hot_spots:
                print(f"  {residue:20s} : {count:5d} interactions")
        else:
            print(f"No residues with ≥{args.min_interactions} interactions found")
    
    # Export CSV
    if args.csv:
        csv_file = f"{args.output}_interactions.csv"
        export_to_csv(interactions, csv_file)
    
    # Generate report
    if args.report:
        report_file = f"{args.output}_report.txt"
        generate_report(analysis, report_file)
    
    print_header("Analysis Complete!")
    
    # Next steps
    print("\nNext steps:")
    print("1. Download experimental PDB structure:")
    print("   python scripts/05_download_pdb.py 1YCR")
    print("\n2. Run Arpeggio on experimental structure:")
    print("   bash scripts/03_run_arpeggio.sh 1ycr.pdb A B arpeggio_experimental")
    print("\n3. Compare results:")
    print("   python scripts/06_compare_structures.py")


if __name__ == '__main__':
    main()
