#!/bin/bash
#
# STEP 7 (BONUS): Protein-Ligand Interaction Analysis with Arpeggio
# ==================================================================
#
# This script demonstrates how to use Arpeggio to analyze interactions
# between a protein and a bound ligand (small molecule).
#
# Example: Analyze MDM2-Nutlin interaction (cancer drug)
#
# Usage:
#   bash 07_protein_ligand_example.sh <pdb_file> [ligand_name] [output_dir]
#
# Example:
#   bash 07_protein_ligand_example.sh 4HG7.pdb NUX arpeggio_ligand
#
# Note: PDB 4HG7 is MDM2 bound to Nutlin-3a (anticancer drug)
#

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
NC='\033[0m'

echo -e "${MAGENTA}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${MAGENTA}║  STEP 7 (BONUS): Protein-Ligand Interaction Analysis         ║${NC}"
echo -e "${MAGENTA}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Parse arguments
PDB_FILE=${1}
LIGAND_NAME=${2:-""}  # If not specified, will analyze all ligands
OUTPUT_DIR=${3:-"arpeggio_ligand"}

if [ -z "$PDB_FILE" ]; then
    echo -e "${RED}ERROR: PDB file required${NC}"
    echo "Usage: bash 07_protein_ligand_example.sh <pdb_file> [ligand_name] [output_dir]"
    echo ""
    echo "Example (MDM2-Nutlin):"
    echo "  # First download PDB structure with ligand"
    echo "  python scripts/05_download_pdb.py 4HG7 ."
    echo "  # Then analyze protein-ligand interactions"
    echo "  bash scripts/07_protein_ligand_example.sh 4hg7.pdb NUX"
    exit 1
fi

if [ ! -f "$PDB_FILE" ]; then
    echo -e "${RED}ERROR: PDB file not found: $PDB_FILE${NC}"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo -e "${GREEN}Input PDB:${NC} $PDB_FILE"
if [ -n "$LIGAND_NAME" ]; then
    echo -e "${GREEN}Ligand:${NC} $LIGAND_NAME"
else
    echo -e "${GREEN}Ligand:${NC} All ligands in structure"
fi
echo -e "${GREEN}Output directory:${NC} $OUTPUT_DIR"
echo ""

# ============================================================================
# Detect ligands in PDB file
# ============================================================================
echo -e "${BLUE}Detecting ligands in PDB file...${NC}"

if command -v python3 &> /dev/null; then
    LIGANDS=$(python3 << EOF
import sys
ligands = set()
with open('$PDB_FILE', 'r') as f:
    for line in f:
        if line.startswith('HETATM'):
            resname = line[17:20].strip()
            # Exclude common water and ions
            if resname not in ['HOH', 'WAT', 'NA', 'CL', 'MG', 'CA', 'K', 'SO4', 'PO4']:
                ligands.add(resname)

if ligands:
    print(' '.join(sorted(ligands)))
else:
    print('NONE')
EOF
    )
    
    if [ "$LIGANDS" != "NONE" ]; then
        echo -e "${GREEN}Found ligands:${NC} $LIGANDS"
        echo ""
        
        # If ligand not specified, use first one
        if [ -z "$LIGAND_NAME" ]; then
            LIGAND_NAME=$(echo $LIGANDS | awk '{print $1}')
            echo -e "${YELLOW}Using ligand: $LIGAND_NAME${NC}"
            echo ""
        fi
    else
        echo -e "${RED}No ligands detected in PDB file${NC}"
        echo "This PDB may not contain any ligands (small molecules)"
        exit 1
    fi
else
    echo -e "${YELLOW}Python not available - skipping ligand detection${NC}"
    if [ -z "$LIGAND_NAME" ]; then
        echo -e "${RED}ERROR: Ligand name required${NC}"
        exit 1
    fi
fi

# ============================================================================
# Run Arpeggio for protein-ligand interactions
# ============================================================================
echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Running Arpeggio for Protein-Ligand Interactions            ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Check Arpeggio
if ! command -v arpeggio &> /dev/null; then
    echo -e "${RED}ERROR: Arpeggio not found${NC}"
    echo "Install: git clone https://github.com/PDBeurope/arpeggio.git && cd arpeggio && python setup.py install"
    exit 1
fi

# Copy PDB to output directory
PDB_BASENAME=$(basename "$PDB_FILE")
cp "$PDB_FILE" "$OUTPUT_DIR/"

# Build selection for ligand
# Format: Select ligand and all nearby protein residues
# Arpeggio will find interactions between them
SELECTION="/${LIGAND_NAME}/"

echo -e "${YELLOW}Selection: $SELECTION${NC}"
echo -e "${YELLOW}This will analyze ALL interactions involving ligand $LIGAND_NAME${NC}"
echo ""

cd "$OUTPUT_DIR"

echo -e "${GREEN}Executing Arpeggio...${NC}"
arpeggio "$PDB_BASENAME" -s "$SELECTION" 2>&1 | tee arpeggio_ligand.log

cd - > /dev/null

# ============================================================================
# Analyze results
# ============================================================================
JSON_FILE="$OUTPUT_DIR/${PDB_BASENAME%.pdb}.json"

if [ ! -f "$JSON_FILE" ]; then
    echo -e "${RED}✗ Analysis may have failed - JSON file not found${NC}"
    exit 1
fi

echo ""
echo -e "${GREEN}✓ Analysis complete!${NC}"
echo ""

# ============================================================================
# Parse and display protein-ligand interactions
# ============================================================================
echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Protein-Ligand Interaction Summary                           ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""

if command -v python3 &> /dev/null; then
    python3 << EOF
import json
from collections import Counter, defaultdict

# Parse JSON
with open('$JSON_FILE', 'r') as f:
    data = json.load(f)

# Extract protein-ligand interactions
ligand_interactions = []
interaction_types = Counter()
residue_interactions = defaultdict(list)

for atom_key, atom_data in data.items():
    if not isinstance(atom_data, dict):
        continue
    
    parts = atom_key.strip('/').split('/')
    if len(parts) < 4:
        continue
    
    chain, resnum, resname, atomname = parts[:4]
    
    # Check if this is ligand or protein
    is_ligand = (resname == '$LIGAND_NAME')
    
    if 'contact' not in atom_data:
        continue
    
    contacts = atom_data['contact']
    if not isinstance(contacts, list):
        contacts = [contacts]
    
    for contact in contacts:
        if not isinstance(contact, dict):
            continue
        
        partner_atom = contact.get('bgn_atom', '')
        itype = contact.get('type', 'unknown')
        distance = contact.get('distance', 0.0)
        
        partner_parts = partner_atom.strip('/').split('/')
        if len(partner_parts) < 4:
            continue
        
        p_chain, p_resnum, p_resname, p_atomname = partner_parts[:4]
        is_partner_ligand = (p_resname == '$LIGAND_NAME')
        
        # We want protein-ligand interactions
        if is_ligand and not is_partner_ligand:
            # Ligand atom contacting protein
            residue = f"{p_chain}:{p_resname}{p_resnum}"
            ligand_atom = f"$LIGAND_NAME-{atomname}"
            protein_atom = f"{residue}-{p_atomname}"
            
            ligand_interactions.append({
                'protein_residue': residue,
                'protein_atom': p_atomname,
                'ligand_atom': atomname,
                'type': itype,
                'distance': float(distance) if distance else 0
            })
            
            interaction_types[itype] += 1
            residue_interactions[residue].append(itype)
        
        elif not is_ligand and is_partner_ligand:
            # Protein atom contacting ligand
            residue = f"{chain}:{resname}{resnum}"
            ligand_atom = f"$LIGAND_NAME-{p_atomname}"
            protein_atom = f"{residue}-{atomname}"
            
            ligand_interactions.append({
                'protein_residue': residue,
                'protein_atom': atomname,
                'ligand_atom': p_atomname,
                'type': itype,
                'distance': float(distance) if distance else 0
            })
            
            interaction_types[itype] += 1
            residue_interactions[residue].append(itype)

# Display results
print(f"Total protein-ligand interactions: {len(ligand_interactions)}")
print()

print("Interaction types:")
print("-" * 50)
for itype, count in sorted(interaction_types.items(), key=lambda x: x[1], reverse=True):
    print(f"  {itype:20s} : {count:5d}")
print()

print("Interacting protein residues:")
print("-" * 50)
for residue, types in sorted(residue_interactions.items(), key=lambda x: len(x[1]), reverse=True):
    print(f"  {residue:20s} : {len(types):2d} interactions ({', '.join(set(types))})")
print()

# Export CSV
import csv
csv_file = '$OUTPUT_DIR/protein_ligand_interactions.csv'
with open(csv_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['protein_residue', 'protein_atom', 'ligand_atom', 'type', 'distance'])
    writer.writeheader()
    writer.writerows(ligand_interactions)

print(f"✓ Interactions exported to: {csv_file}")
print()
EOF
else
    echo -e "${YELLOW}Python not available - raw JSON only${NC}"
fi

# ============================================================================
# Output files
# ============================================================================
echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Output Files                                                  ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${GREEN}Generated files:${NC}"
[ -f "$JSON_FILE" ] && echo -e "  ✓ $JSON_FILE"
[ -f "$OUTPUT_DIR/protein_ligand_interactions.csv" ] && echo -e "  ✓ $OUTPUT_DIR/protein_ligand_interactions.csv"
[ -f "$OUTPUT_DIR/arpeggio_ligand.log" ] && echo -e "  ✓ $OUTPUT_DIR/arpeggio_ligand.log"
echo ""

# ============================================================================
# Visualization suggestions
# ============================================================================
echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Visualize Interactions (PyMOL)                               ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""

cat > "$OUTPUT_DIR/visualize_ligand.pml" << 'PYMOL_SCRIPT'
# PyMOL script to visualize protein-ligand interactions
# Load structure
load {PDB_FILE}

# Clean up
remove solvent
remove resn HOH+WAT

# Color protein
color gray80, polymer
show cartoon, polymer

# Show and color ligand
show sticks, resn {LIGAND_NAME}
color yellow, resn {LIGAND_NAME} and elem C
color red, resn {LIGAND_NAME} and elem O
color blue, resn {LIGAND_NAME} and elem N
util.cnc resn {LIGAND_NAME}

# Show binding site residues
select binding_site, polymer within 5 of resn {LIGAND_NAME}
show sticks, binding_site
color cyan, binding_site and elem C

# Label binding site
label binding_site and name CA, "%s%s" % (resn, resi)

# Center view
orient resn {LIGAND_NAME}
zoom resn {LIGAND_NAME}, 8

# Nice rendering
bg_color white
set ray_shadows, 0
set antialias, 2

# Find hydrogen bonds
distance hbonds, polymer, resn {LIGAND_NAME}, 3.5, mode=2
color yellow, hbonds

# Save image
ray 2400, 2400
png protein_ligand_complex.png, dpi=300

# Save session
save protein_ligand_session.pse

print "Visualization complete!"
PYMOL_SCRIPT

# Substitute variables
sed -i "s|{PDB_FILE}|$PDB_FILE|g" "$OUTPUT_DIR/visualize_ligand.pml"
sed -i "s|{LIGAND_NAME}|$LIGAND_NAME|g" "$OUTPUT_DIR/visualize_ligand.pml"

echo -e "PyMOL script generated: ${GREEN}$OUTPUT_DIR/visualize_ligand.pml${NC}"
echo ""
echo -e "Run visualization:"
echo -e "  ${YELLOW}pymol $OUTPUT_DIR/visualize_ligand.pml${NC}"
echo ""

# ============================================================================
# Example: MDM2-Nutlin
# ============================================================================
if [[ "$PDB_FILE" == *"4hg7"* ]] || [[ "$PDB_FILE" == *"4HG7"* ]]; then
    echo -e "${MAGENTA}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${MAGENTA}║  MDM2-Nutlin-3a Complex (Anticancer Drug)                    ║${NC}"
    echo -e "${MAGENTA}╚════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    echo -e "${GREEN}About this complex:${NC}"
    echo -e "  • Nutlin-3a = MDM2 inhibitor (cancer drug)"
    echo -e "  • Mimics p53 residues Phe19, Trp23, Leu26"
    echo -e "  • Blocks MDM2-p53 interaction → reactivates p53"
    echo -e "  • Clinical analogs: Idasanutlin (Phase III trials)"
    echo ""
    echo -e "${GREEN}Expected key interactions:${NC}"
    echo -e "  • Hydrophobic contacts with Leu54, Met62, Ile61, Val93"
    echo -e "  • Multiple Van der Waals contacts (20+)"
    echo -e "  • Possible hydrogen bond with Leu54 backbone"
    echo ""
fi

echo -e "${GREEN}✓ Protein-ligand analysis complete!${NC}"
echo ""
