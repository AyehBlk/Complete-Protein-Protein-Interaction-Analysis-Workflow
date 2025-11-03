#!/usr/bin/env python3
"""
STEP 5: Download Experimental PDB Structure
============================================

This script downloads experimental structures from the PDB database
for comparison with AlphaFold predictions.

Usage:
    python 05_download_pdb.py <pdb_id> [output_dir]

Example:
    python 05_download_pdb.py 1YCR experimental_structures/

Author: [Your Name]
Date: October 31, 2025
"""

import sys
import argparse
from pathlib import Path
from urllib.request import urlretrieve
import gzip
import shutil


def print_header(text):
    """Print formatted header"""
    print("\n" + "="*70)
    print(f"  {text}")
    print("="*70 + "\n")


def download_pdb(pdb_id, output_dir='.', format='pdb'):
    """
    Download PDB structure
    
    Args:
        pdb_id: 4-letter PDB code
        output_dir: Output directory
        format: 'pdb' or 'cif'
    
    Returns:
        str: Path to downloaded file
    """
    pdb_id = pdb_id.upper()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Downloading PDB: {pdb_id}")
    print(f"Output directory: {output_dir}")
    
    # Determine URL and filename
    if format == 'pdb':
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb.gz"
        output_file_gz = output_dir / f"{pdb_id}.pdb.gz"
        output_file = output_dir / f"{pdb_id.lower()}.pdb"
    elif format == 'cif':
        url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
        output_file_gz = output_dir / f"{pdb_id}.cif.gz"
        output_file = output_dir / f"{pdb_id.lower()}.cif"
    else:
        raise ValueError(f"Unknown format: {format}")
    
    try:
        # Download compressed file
        print(f"Fetching from: {url}")
        urlretrieve(url, output_file_gz)
        print(f"✓ Downloaded: {output_file_gz}")
        
        # Decompress
        print(f"Decompressing...")
        with gzip.open(output_file_gz, 'rb') as f_in:
            with open(output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"✓ Decompressed: {output_file}")
        
        # Remove compressed file
        output_file_gz.unlink()
        
        return str(output_file)
        
    except Exception as e:
        print(f"ERROR: Failed to download {pdb_id}: {e}")
        return None


def get_pdb_info(pdb_id):
    """
    Get information about PDB entry
    
    Args:
        pdb_id: 4-letter PDB code
    """
    import json
    from urllib.request import urlopen
    
    pdb_id = pdb_id.upper()
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    try:
        with urlopen(url) as response:
            data = json.load(response)
        
        print(f"\nPDB {pdb_id} Information:")
        print("-" * 70)
        
        # Title
        if 'struct' in data and 'title' in data['struct']:
            print(f"Title: {data['struct']['title']}")
        
        # Release date
        if 'rcsb_accession_info' in data:
            if 'deposit_date' in data['rcsb_accession_info']:
                print(f"Deposit date: {data['rcsb_accession_info']['deposit_date']}")
        
        # Resolution
        if 'rcsb_entry_info' in data:
            if 'resolution_combined' in data['rcsb_entry_info']:
                res = data['rcsb_entry_info']['resolution_combined']
                print(f"Resolution: {res[0]:.2f} Å")
        
        # Experimental method
        if 'exptl' in data and len(data['exptl']) > 0:
            method = data['exptl'][0].get('method', 'Unknown')
            print(f"Method: {method}")
        
        print()
        
    except Exception as e:
        print(f"Could not fetch info for {pdb_id}: {e}")


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(
        description='Download experimental PDB structure',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download MDM2/p53 structure (1YCR)
  python 05_download_pdb.py 1YCR
  
  # Download to specific directory
  python 05_download_pdb.py 1YCR experimental_structures/
  
  # Download CIF format
  python 05_download_pdb.py 1YCR --format cif
        """
    )
    
    parser.add_argument('pdb_id', help='4-letter PDB code')
    parser.add_argument('output_dir', nargs='?', default='.',
                       help='Output directory (default: current directory)')
    parser.add_argument('--format', choices=['pdb', 'cif'], default='pdb',
                       help='File format (default: pdb)')
    parser.add_argument('--info', action='store_true',
                       help='Show PDB entry information')
    
    args = parser.parse_args()
    
    print_header("PDB Structure Download")
    
    # Show info if requested
    if args.info:
        get_pdb_info(args.pdb_id)
    
    # Download structure
    output_file = download_pdb(args.pdb_id, args.output_dir, args.format)
    
    if output_file:
        print_header("Download Complete!")
        print(f"Structure saved to: {output_file}")
        print(f"\nNext steps:")
        print(f"1. Visualize structure:")
        print(f"   pymol {output_file}")
        print(f"\n2. Run Arpeggio analysis:")
        print(f"   bash scripts/03_run_arpeggio.sh {output_file} A B arpeggio_experimental")
        print(f"\n3. Compare with prediction:")
        print(f"   python scripts/06_compare_structures.py")
    else:
        print("\nDownload failed!")
        sys.exit(1)


if __name__ == '__main__':
    main()
