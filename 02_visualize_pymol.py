#!/usr/bin/env python3
"""
STEP 2: PyMOL Visualization of Protein-Protein Interface
=========================================================

This script creates beautiful visualizations of protein-protein interactions
using PyMOL. It generates:
- Interface visualization with chains colored differently
- Surface representation showing binding interface
- Detailed view of key interactions
- High-quality images for publications
- PyMOL session file for interactive exploration

Usage:
    python 02_visualize_pymol.py <pdb_file> [options]

Example:
    python 02_visualize_pymol.py ../examples/mdm2_p53/results/alphafold/fold_model_0.pdb \
        --chain1 A --chain2 B --output ../examples/mdm2_p53/results/figures

Author: AyehBlk
Date: October 31, 2025
"""

import sys
import os
from pathlib import Path
import argparse

try:
    import pymol
    from pymol import cmd, stored
    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False
    print("WARNING: PyMOL not available. Will generate PyMOL script instead.")


def print_header(text):
    """Print formatted header"""
    print("\n" + "="*70)
    print(f"  {text}")
    print("="*70 + "\n")


def visualize_interface_pymol(pdb_file, chain1='A', chain2='B', output_dir='figures', 
                              interface_distance=5.0):
    """
    Create interface visualization using PyMOL
    
    Args:
        pdb_file: Path to PDB file
        chain1: First chain ID
        chain2: Second chain ID
        output_dir: Output directory for images
        interface_distance: Distance to define interface (Angstroms)
    """
    if not PYMOL_AVAILABLE:
        print("ERROR: PyMOL is not installed.")
        print("Install with: conda install -c conda-forge pymol-open-source")
        return
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    pdb_name = Path(pdb_file).stem
    
    print_header(f"Visualizing: {pdb_name}")
    
    # Initialize PyMOL
    pymol.finish_launching(['pymol', '-c'])  # -c for command-line mode
    
    # Load structure
    print(f"Loading structure: {pdb_file}")
    cmd.load(pdb_file, 'protein')
    
    # =================================================================
    # View 1: Overview with both chains
    # =================================================================
    print("\n1. Creating overview image...")
    
    cmd.hide('everything', 'protein')
    cmd.show('cartoon', 'protein')
    cmd.color('cyan', f'chain {chain1}')
    cmd.color('magenta', f'chain {chain2}')
    cmd.bg_color('white')
    cmd.set('ray_shadows', 0)
    cmd.set('antialias', 2)
    
    cmd.orient('protein')
    cmd.zoom('protein', 5)
    
    output_file = output_dir / f'{pdb_name}_overview.png'
    cmd.ray(2400, 2400)
    cmd.png(str(output_file), dpi=300)
    print(f"   Saved: {output_file}")
    
    # =================================================================
    # View 2: Interface close-up
    # =================================================================
    print("\n2. Creating interface close-up...")
    
    # Define interface residues
    cmd.select('interface1', f'chain {chain1} within {interface_distance} of chain {chain2}')
    cmd.select('interface2', f'chain {chain2} within {interface_distance} of chain {chain1}')
    cmd.select('interface', 'interface1 or interface2')
    
    # Show interface in sticks
    cmd.hide('everything', 'protein')
    cmd.show('cartoon', 'protein')
    cmd.show('sticks', 'interface')
    cmd.color('cyan', f'chain {chain1}')
    cmd.color('magenta', f'chain {chain2}')
    cmd.color('yellow', 'interface and elem C')
    
    cmd.orient('interface')
    cmd.zoom('interface', 8)
    
    output_file = output_dir / f'{pdb_name}_interface.png'
    cmd.ray(2400, 2400)
    cmd.png(str(output_file), dpi=300)
    print(f"   Saved: {output_file}")
    
    # =================================================================
    # View 3: Surface representation
    # =================================================================
    print("\n3. Creating surface view...")
    
    cmd.hide('everything', 'protein')
    cmd.show('surface', f'chain {chain1}')
    cmd.show('cartoon', f'chain {chain2}')
    cmd.color('cyan', f'chain {chain1}')
    cmd.color('magenta', f'chain {chain2}')
    cmd.set('transparency', 0.3, f'chain {chain1}')
    
    cmd.orient('protein')
    cmd.zoom('protein', 5)
    
    output_file = output_dir / f'{pdb_name}_surface.png'
    cmd.ray(2400, 2400)
    cmd.png(str(output_file), dpi=300)
    print(f"   Saved: {output_file}")
    
    # =================================================================
    # View 4: Interface with labeled key residues (optional)
    # =================================================================
    print("\n4. Creating labeled interface...")
    
    cmd.hide('everything', 'protein')
    cmd.show('cartoon', 'protein')
    cmd.show('sticks', 'interface')
    cmd.color('cyan', f'chain {chain1}')
    cmd.color('magenta', f'chain {chain2}')
    cmd.color('yellow', 'interface and elem C')
    
    # Label interface residues
    cmd.label('interface and name CA', '"%s%s" % (resn, resi)')
    cmd.set('label_size', 20)
    cmd.set('label_color', 'black')
    
    cmd.orient('interface')
    cmd.zoom('interface', 8)
    
    output_file = output_dir / f'{pdb_name}_interface_labeled.png'
    cmd.ray(2400, 2400)
    cmd.png(str(output_file), dpi=300)
    print(f"   Saved: {output_file}")
    
    # =================================================================
    # Save PyMOL session
    # =================================================================
    print("\n5. Saving PyMOL session...")
    
    cmd.hide('labels')  # Hide labels for session
    session_file = output_dir / f'{pdb_name}_session.pse'
    cmd.save(str(session_file))
    print(f"   Saved: {session_file}")
    
    # Clean up
    cmd.delete('all')
    
    print_header("Visualization Complete!")
    print(f"All files saved to: {output_dir}/")
    print(f"\nTo explore interactively:")
    print(f"  pymol {session_file}")


def generate_pymol_script(pdb_file, chain1='A', chain2='B', output_file='visualize.pml',
                         interface_distance=5.0):
    """
    Generate PyMOL script file for manual visualization
    
    Args:
        pdb_file: Path to PDB file
        chain1: First chain ID
        chain2: Second chain ID
        output_file: Output PyMOL script file
        interface_distance: Distance to define interface
    """
    print_header("Generating PyMOL Script")
    
    script = f"""# PyMOL Visualization Script
# Generated for: {pdb_file}
# Date: October 31, 2025

# Load structure
load {pdb_file}, protein

# Basic setup
bg_color white
set ray_shadows, 0
set antialias, 2

# ============================================================================
# View 1: Overview
# ============================================================================
hide everything, protein
show cartoon, protein
color cyan, chain {chain1}
color magenta, chain {chain2}
orient protein
zoom protein, 5

# Save image
ray 2400, 2400
png overview.png, dpi=300

# ============================================================================
# Export and View 2: Interface
# ============================================================================
select interface1, chain {chain1} within {interface_distance} of chain {chain2}
select interface2, chain {chain2} within {interface_distance} of chain {chain1}
select interface, interface1 or interface2

##Save 'interface' object to a Text File

# Create a list and save to file
stored.residues = []
iterate interface and name CA, stored.residues.append(f"{chain}\t{resn}\t{resi}")

# Remove duplicates and sort
stored.residues = sorted(set(stored.residues))

# Save to file on your PC
with open('interface_residues.txt', 'w') as f:
    f.write("Chain\tResidue\tNumber\n")
    for res in stored.residues:
        f.write(res + "\n")

print(f"Saved {len(stored.residues)} interface residues")

##view and export image

hide everything, protein
show cartoon, protein
show sticks, interface
color cyan, chain {chain1}
color magenta, chain {chain2}
color yellow, interface and elem C

orient interface
zoom interface, 8

# Save image
ray 2400, 2400
png interface.png, dpi=300

# ============================================================================
# View 3: Surface
# ============================================================================
hide everything, protein
show surface, chain {chain1}
show cartoon, chain {chain2}
color cyan, chain {chain1}
color magenta, chain {chain2}
set transparency, 0.3, chain {chain1}

orient protein
zoom protein, 5

# Save image
ray 2400, 2400
png surface.png, dpi=300

# ============================================================================
# Measurements (optional)
# ============================================================================
# Measure distance between key residues (example)
# distance dist1, chain {chain1} and resi 19 and name CA, chain {chain2} and resi 67 and name CA

# ============================================================================
# Color by B-factor (AlphaFold confidence)
# ============================================================================
# For AlphaFold predictions, B-factor = pLDDT
hide everything
show cartoon, protein
spectrum b, blue_white_red, protein, minimum=50, maximum=100
orient protein

# Save confidence view
ray 2400, 2400
png confidence.png, dpi=300

# ============================================================================
# Interactive commands
# ============================================================================
# Uncomment and modify as needed:
# label interface and name CA, "%s%s" % (resn, resi)
# show sticks, resi 19+23+26 and chain {chain2}  # p53 hot spots
# show sticks, resi 54+57+67 and chain {chain1}  # MDM2 binding pocket

# Save session
save protein_session.pse

print "Visualization complete!"
print "Images saved: overview.png, interface.png, surface.png, confidence.png"
print "Session saved: protein_session.pse"
"""
    
    with open(output_file, 'w') as f:
        f.write(script)
    
    print(f"PyMOL script saved: {output_file}")
    print(f"\nTo use this script:")
    print(f"  1. Open PyMOL")
    print(f"  2. File -> Run Script -> {output_file}")
    print(f"  Or from command line:")
    print(f"     pymol {output_file}")


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(
        description='Visualize protein-protein interface with PyMOL',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python 02_visualize_pymol.py structure.pdb
  
  # Specify chains and output
  python 02_visualize_pymol.py structure.pdb --chain1 A --chain2 B --output figures/
  
  # Generate script only (if PyMOL not installed)
  python 02_visualize_pymol.py structure.pdb --script-only
        """
    )
    
    parser.add_argument('pdb_file', help='Input PDB file')
    parser.add_argument('--chain1', default='A', help='First chain ID (default: A)')
    parser.add_argument('--chain2', default='B', help='Second chain ID (default: B)')
    parser.add_argument('--output', default='figures', help='Output directory (default: figures)')
    parser.add_argument('--interface-distance', type=float, default=5.0,
                       help='Interface distance cutoff in Angstroms (default: 5.0)')
    parser.add_argument('--script-only', action='store_true',
                       help='Generate PyMOL script only (don\'t run PyMOL)')
    
    args = parser.parse_args()
    
    # Check if PDB file exists
    if not os.path.exists(args.pdb_file):
        print(f"ERROR: PDB file not found: {args.pdb_file}")
        sys.exit(1)
    
    print_header("PyMOL Protein Interface Visualization")
    print(f"Structure: {args.pdb_file}")
    print(f"Chain 1: {args.chain1}")
    print(f"Chain 2: {args.chain2}")
    print(f"Interface distance: {args.interface_distance} Å")
    print(f"Output directory: {args.output}")
    
    if args.script_only or not PYMOL_AVAILABLE:
        # Generate script file
        script_file = Path(args.output) / 'visualize.pml'
        Path(args.output).mkdir(parents=True, exist_ok=True)
        generate_pymol_script(
            args.pdb_file,
            args.chain1,
            args.chain2,
            str(script_file),
            args.interface_distance
        )
    else:
        # Run PyMOL visualization
        visualize_interface_pymol(
            args.pdb_file,
            args.chain1,
            args.chain2,
            args.output,
            args.interface_distance
        )
    
    print_header("Done!")


if __name__ == '__main__':
    main()
