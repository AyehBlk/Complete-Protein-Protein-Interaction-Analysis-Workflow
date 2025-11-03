#!/bin/bash
#
# STEP 1: AlphaFold3 Heterodimer Prediction
# ==========================================
#
# This script prepares input for AlphaFold3 and provides commands
# to run prediction for a heterodimer (two protein chains).
#
# Example: MDM2/p53 complex
#
# Usage:
#   bash 01_run_alphafold.sh <seq1_fasta> <seq2_fasta> <output_dir>
#
# Example:
#   bash 01_run_alphafold.sh ../examples/mdm2_p53/sequences/mdm2.fasta \
#                            ../examples/mdm2_p53/sequences/p53.fasta \
#                            ../examples/mdm2_p53/results/alphafold
#

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  STEP 1: AlphaFold3 Heterodimer Prediction                    ║${NC}"
echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Parse arguments
SEQ1_FASTA=${1:-"../examples/mdm2_p53/sequences/mdm2.fasta"}
SEQ2_FASTA=${2:-"../examples/mdm2_p53/sequences/p53.fasta"}
OUTPUT_DIR=${3:-"../examples/mdm2_p53/results/alphafold"}

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo -e "${GREEN}Input sequences:${NC}"
echo -e "  Chain A: $SEQ1_FASTA"
echo -e "  Chain B: $SEQ2_FASTA"
echo -e "${GREEN}Output directory:${NC} $OUTPUT_DIR"
echo ""

# Check if files exist
if [ ! -f "$SEQ1_FASTA" ]; then
    echo -e "${RED}ERROR: Sequence file not found: $SEQ1_FASTA${NC}"
    exit 1
fi

if [ ! -f "$SEQ2_FASTA" ]; then
    echo -e "${RED}ERROR: Sequence file not found: $SEQ2_FASTA${NC}"
    exit 1
fi

# Read sequences
echo -e "${BLUE}Reading sequences...${NC}"
SEQ1=$(grep -v "^>" "$SEQ1_FASTA" | tr -d '\n')
SEQ2=$(grep -v "^>" "$SEQ2_FASTA" | tr -d '\n')

echo -e "  Chain A length: ${#SEQ1} residues"
echo -e "  Chain B length: ${#SEQ2} residues"
echo ""

# Create combined FASTA for AlphaFold3
COMBINED_FASTA="$OUTPUT_DIR/heterodimer_input.fasta"
echo -e "${BLUE}Creating combined input file: $COMBINED_FASTA${NC}"

cat > "$COMBINED_FASTA" << EOF
>ChainA
$SEQ1
>ChainB
$SEQ2
EOF

echo -e "${GREEN}✓ Combined FASTA created${NC}"
echo ""

# Create JSON input for AlphaFold3 (if using local installation)
JSON_INPUT="$OUTPUT_DIR/alphafold_input.json"
echo -e "${BLUE}Creating JSON input: $JSON_INPUT${NC}"

cat > "$JSON_INPUT" << EOF
{
  "name": "heterodimer_prediction",
  "sequences": [
    {
      "proteinChain": {
        "sequence": "$SEQ1",
        "count": 1
      }
    },
    {
      "proteinChain": {
        "sequence": "$SEQ2",
        "count": 1
      }
    }
  ],
  "modelSeeds": [1, 2, 3, 4, 5]
}
EOF

echo -e "${GREEN}✓ JSON input created${NC}"
echo ""

# ============================================================================
# Option 1: AlphaFold Server (Recommended for most users)
# ============================================================================
echo -e "${YELLOW}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${YELLOW}║  OPTION 1: AlphaFold Server (Recommended)                     ║${NC}"
echo -e "${YELLOW}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "1. Go to: ${BLUE}https://alphafoldserver.com/${NC}"
echo -e "2. Sign in (free for academic use)"
echo -e "3. Click 'New Prediction'"
echo -e "4. Enter:"
echo -e "   - Name: ${GREEN}MDM2_p53_heterodimer${NC}"
echo -e "   - Paste sequences from: ${GREEN}$COMBINED_FASTA${NC}"
echo -e "   - Or enter manually (use $COMBINED_FASTA content)"
echo -e "5. Click 'Submit'"
echo -e "6. Wait for prediction (~minutes to hours)"
echo -e "7. Download results:"
echo -e "   - ${GREEN}fold_model_0.pdb${NC} (best model)"
echo -e "   - ${GREEN}confidences.json${NC} (pLDDT, PAE)"
echo -e "8. Save to: ${GREEN}$OUTPUT_DIR/${NC}"
echo ""

# ============================================================================
# Option 2: Local AlphaFold3 (Requires installation)
# ============================================================================
echo -e "${YELLOW}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${YELLOW}║  OPTION 2: Local AlphaFold3 Installation                      ║${NC}"
echo -e "${YELLOW}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${RED}Note: Requires AlphaFold3 installation and GPU${NC}"
echo ""
echo -e "If you have AlphaFold3 installed locally, run:"
echo ""
echo -e "${GREEN}# Using Docker (recommended)${NC}"
echo -e "docker run --gpus all \\"
echo -e "  -v \$(pwd)/$OUTPUT_DIR:/data \\"
echo -e "  alphafold3 \\"
echo -e "  --json_path=/data/alphafold_input.json \\"
echo -e "  --output_dir=/data \\"
echo -e "  --model_dir=/models"
echo ""
echo -e "${GREEN}# Or using Python directly${NC}"
echo -e "python run_alphafold.py \\"
echo -e "  --json_path=$JSON_INPUT \\"
echo -e "  --output_dir=$OUTPUT_DIR \\"
echo -e "  --model_dir=/path/to/models"
echo ""

# ============================================================================
# Option 3: ColabFold (AlphaFold2 - still useful)
# ============================================================================
echo -e "${YELLOW}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${YELLOW}║  OPTION 3: ColabFold (AlphaFold2 - Alternative)               ║${NC}"
echo -e "${YELLOW}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "For quick testing (uses AlphaFold2):"
echo -e "1. Go to: ${BLUE}https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb${NC}"
echo -e "2. Upload: $COMBINED_FASTA"
echo -e "3. Run notebook"
echo -e "4. Download results"
echo ""

# ============================================================================
# What to expect
# ============================================================================
echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Expected Output Files                                        ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "After prediction completes, you should have:"
echo -e "  ${GREEN}fold_model_0.pdb${NC}        - Predicted structure (best model)"
echo -e "  ${GREEN}fold_model_1-4.pdb${NC}      - Alternative models"
echo -e "  ${GREEN}confidences.json${NC}        - Quality metrics"
echo -e "  ${GREEN}ranking_scores.json${NC}     - Model ranking"
echo -e ""
echo -e "Quality indicators:"
echo -e "  ${GREEN}pLDDT > 90${NC}              - Very high confidence"
echo -e "  ${GREEN}pLDDT 70-90${NC}             - Good confidence"
echo -e "  ${GREEN}PAE < 5 Å${NC}               - Excellent inter-chain prediction"
echo -e "  ${GREEN}PAE 5-10 Å${NC}              - Good inter-chain prediction"
echo ""

# ============================================================================
# Next steps
# ============================================================================
echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Next Steps                                                    ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "Once you have the predicted structure:"
echo -e ""
echo -e "1. ${GREEN}Visualize with PyMOL:${NC}"
echo -e "   python scripts/02_visualize_pymol.py $OUTPUT_DIR/fold_model_0.pdb"
echo -e ""
echo -e "2. ${GREEN}Analyze interactions with Arpeggio:${NC}"
echo -e "   bash scripts/03_run_arpeggio.sh $OUTPUT_DIR/fold_model_0.pdb"
echo -e ""
echo -e "3. ${GREEN}Run complete workflow:${NC}"
echo -e "   bash scripts/run_complete_workflow.sh"
echo ""

# ============================================================================
# Example for MDM2/p53
# ============================================================================
echo -e "${YELLOW}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${YELLOW}║  Example: MDM2/p53 Complex                                    ║${NC}"
echo -e "${YELLOW}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "Expected results for MDM2/p53:"
echo -e "  • Interface pLDDT: ${GREEN}~87${NC} (high confidence)"
echo -e "  • Interface PAE: ${GREEN}~5.2 Å${NC} (excellent)"
echo -e "  • RMSD vs PDB 1YCR: ${GREEN}~0.85 Å${NC} (excellent)"
echo -e "  • Key residues: ${GREEN}Phe19, Trp23, Leu26${NC} (p53)"
echo ""
echo -e "${GREEN}✓ Preparation complete!${NC}"
echo -e "${YELLOW}Now run AlphaFold3 prediction using one of the options above.${NC}"
echo ""
