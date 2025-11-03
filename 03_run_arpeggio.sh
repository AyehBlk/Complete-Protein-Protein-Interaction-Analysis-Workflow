#!/bin/bash
#
# STEP 3: Arpeggio Analysis for Protein-Protein Interactions
# ===========================================================
#
# This script runs Arpeggio to detect all interactions between two protein chains.
# Arpeggio detects:
#   - Hydrogen bonds
#   - Van der Waals contacts
#   - Aromatic interactions
#   - Ionic interactions
#   - Hydrophobic contacts
#   - Carbonyl interactions
#   - And more...
#
# Usage:
#   bash 03_run_arpeggio.sh <pdb_file> [chain1] [chain2] [output_dir]
#
# Example:
#   bash 03_run_arpeggio.sh ../examples/mdm2_p53/results/alphafold/fold_model_0.pdb A B \
#       ../examples/mdm2_p53/results/arpeggio_predicted
#

set -e  # Exit on error

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  STEP 3: Arpeggio Protein-Protein Interaction Analysis        ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Parse arguments
PDB_FILE=${1}
CHAIN1=${2:-"A"}
CHAIN2=${3:-"B"}
OUTPUT_DIR=${4:-"arpeggio_results"}

if [ -z "$PDB_FILE" ]; then
    echo -e "${RED}ERROR: PDB file required${NC}"
    echo "Usage: bash 03_run_arpeggio.sh <pdb_file> [chain1] [chain2] [output_dir]"
    exit 1
fi

if [ ! -f "$PDB_FILE" ]; then
    echo -e "${RED}ERROR: PDB file not found: $PDB_FILE${NC}"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo -e "${GREEN}Input PDB:${NC} $PDB_FILE"
echo -e "${GREEN}Chains:${NC} $CHAIN1 and $CHAIN2"
echo -e "${GREEN}Output directory:${NC} $OUTPUT_DIR"
echo ""

# Check if Arpeggio is installed
if ! command -v arpeggio &> /dev/null; then
    echo -e "${RED}ERROR: Arpeggio not found in PATH${NC}"
    echo ""
    echo "Please install Arpeggio:"
    echo "  git clone https://github.com/PDBeurope/arpeggio.git"
    echo "  cd arpeggio"
    echo "  python setup.py install"
    echo ""
    exit 1
fi

echo -e "${BLUE}Checking Arpeggio installation...${NC}"
ARPEGGIO_VERSION=$(arpeggio --version 2>&1 || echo "unknown")
echo -e "  Version: $ARPEGGIO_VERSION"
echo ""

# ============================================================================
# Run Arpeggio
# ============================================================================
echo -e "${BLUE}Running Arpeggio analysis...${NC}"
echo ""

# Build selection string for inter-chain interactions
# Format: /A/,/B/ to analyze interactions between chains A and B
SELECTION="/${CHAIN1}/,/${CHAIN2}/"

echo -e "${YELLOW}Selection: $SELECTION${NC}"
echo -e "${YELLOW}This will analyze interactions BETWEEN chains $CHAIN1 and $CHAIN2${NC}"
echo ""

# Copy PDB to output directory
PDB_BASENAME=$(basename "$PDB_FILE")
cp "$PDB_FILE" "$OUTPUT_DIR/"

# Run Arpeggio
# -s: selection (which chains/residues to analyze)
# Output files will be created in current directory, so cd to output dir
cd "$OUTPUT_DIR"

echo -e "${GREEN}Executing Arpeggio...${NC}"
arpeggio "$PDB_BASENAME" -s "$SELECTION" 2>&1 | tee arpeggio.log

cd - > /dev/null  # Return to original directory

# Check if output files were created
JSON_FILE="$OUTPUT_DIR/${PDB_BASENAME%.pdb}.json"
CONTACTS_FILE="$OUTPUT_DIR/${PDB_BASENAME%.pdb}.contacts"

if [ -f "$JSON_FILE" ]; then
    echo -e "${GREEN}✓ Analysis complete!${NC}"
else
    echo -e "${RED}✗ Analysis may have failed - JSON file not found${NC}"
    echo -e "Check log file: $OUTPUT_DIR/arpeggio.log"
    exit 1
fi

# ============================================================================
# Display results summary
# ============================================================================
echo ""
echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Output Files                                                  ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""

ls -lh "$OUTPUT_DIR" | grep -E "\.(json|contacts|rings|amide)" || true

echo ""
echo -e "${GREEN}Generated files:${NC}"
[ -f "$JSON_FILE" ] && echo -e "  ✓ ${GREEN}$JSON_FILE${NC} - Complete interaction data (JSON format)"
[ -f "$CONTACTS_FILE" ] && echo -e "  ✓ ${GREEN}$CONTACTS_FILE${NC} - Contact list (text format)"
[ -f "$OUTPUT_DIR/${PDB_BASENAME%.pdb}.rings" ] && echo -e "  ✓ Aromatic ring data"
[ -f "$OUTPUT_DIR/${PDB_BASENAME%.pdb}.amide" ] && echo -e "  ✓ Amide plane data"
[ -f "$OUTPUT_DIR/arpeggio.log" ] && echo -e "  ✓ Analysis log"
echo ""

# ============================================================================
# Quick summary using Python (if available)
# ============================================================================
if command -v python3 &> /dev/null && [ -f "$JSON_FILE" ]; then
    echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  Quick Summary                                                 ║${NC}"
    echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    
    python3 << EOF
import json
import sys
from collections import Counter

try:
    with open('$JSON_FILE', 'r') as f:
        data = json.load(f)
    
    # Count interactions by type
    interaction_counts = Counter()
    total_interactions = 0
    
    for atom_key, atom_data in data.items():
        if isinstance(atom_data, dict) and 'contact' in atom_data:
            contacts = atom_data['contact']
            for contact in contacts:
                if isinstance(contact, dict) and 'type' in contact:
                    interaction_type = contact['type']
                    interaction_counts[interaction_type] += 1
                    total_interactions += 1
    
    print(f"Total interactions found: {total_interactions}")
    print()
    print("Breakdown by type:")
    print("-" * 50)
    for itype, count in sorted(interaction_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"  {itype:20s} : {count:5d}")
    print()
    
except Exception as e:
    print(f"Could not parse JSON: {e}")
    sys.exit(0)
EOF
else
    echo -e "${YELLOW}Install Python to see interaction summary${NC}"
fi

# ============================================================================
# Next steps
# ============================================================================
echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Next Steps                                                    ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "1. ${GREEN}Analyze interactions in detail:${NC}"
echo -e "   python scripts/04_analyze_interactions.py $JSON_FILE"
echo ""
echo -e "2. ${GREEN}Download PDB structure for comparison:${NC}"
echo -e "   python scripts/05_download_pdb.py 1YCR"
echo ""
echo -e "3. ${GREEN}Run Arpeggio on experimental structure:${NC}"
echo -e "   bash scripts/03_run_arpeggio.sh 1ycr.pdb $CHAIN1 $CHAIN2 arpeggio_experimental"
echo ""
echo -e "4. ${GREEN}Compare predicted vs experimental:${NC}"
echo -e "   python scripts/06_compare_structures.py"
echo ""

# ============================================================================
# Interpretation guide
# ============================================================================
echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Understanding Arpeggio Output                                 ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${GREEN}Interaction Types:${NC}"
echo -e "  • ${YELLOW}clash${NC}          - Steric clash (bad)"
echo -e "  • ${YELLOW}covalent${NC}       - Covalent bond"
echo -e "  • ${YELLOW}vdw_clash${NC}      - Van der Waals clash"
echo -e "  • ${YELLOW}vdw${NC}            - Van der Waals contact (good)"
echo -e "  • ${YELLOW}proximal${NC}       - Close proximity"
echo -e "  • ${YELLOW}hbond${NC}          - Hydrogen bond (important!)"
echo -e "  • ${YELLOW}weak_hbond${NC}     - Weak hydrogen bond"
echo -e "  • ${YELLOW}xbond${NC}          - Halogen bond"
echo -e "  • ${YELLOW}ionic${NC}          - Ionic interaction (salt bridge)"
echo -e "  • ${YELLOW}metal_complex${NC}  - Metal coordination"
echo -e "  • ${YELLOW}aromatic${NC}       - Aromatic interaction (π-π)"
echo -e "  • ${YELLOW}hydrophobic${NC}    - Hydrophobic contact"
echo -e "  • ${YELLOW}carbonyl${NC}       - Carbonyl interaction"
echo ""
echo -e "${GREEN}Good interfaces typically have:${NC}"
echo -e "  • Multiple hydrogen bonds (5-10+)"
echo -e "  • Many Van der Waals contacts (20-50+)"
echo -e "  • Some aromatic interactions"
echo -e "  • Hydrophobic core"
echo -e "  • Few or no clashes"
echo ""

echo -e "${GREEN}✓ Arpeggio analysis complete!${NC}"
echo ""
