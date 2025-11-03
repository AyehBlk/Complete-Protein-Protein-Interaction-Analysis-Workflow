#!/bin/bash
#
# Complete Workflow: Protein-Protein Interaction Analysis
# =========================================================
#
# This master script runs the complete workflow from AlphaFold prediction
# to validation against experimental structure.
#
# Usage:
#   bash run_complete_workflow.sh
#
# Or with custom files:
#   bash run_complete_workflow.sh \
#       --seq1 protein1.fasta \
#       --seq2 protein2.fasta \
#       --pdb-id 1ABC \
#       --output results

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m'

# Default parameters (MDM2/p53 example)
SEQ1="../examples/mdm2_p53/sequences/mdm2.fasta"
SEQ2="../examples/mdm2_p53/sequences/p53.fasta"
PDB_ID="1YCR"
CHAIN1="A"
CHAIN2="B"
OUTPUT_BASE="../examples/mdm2_p53/results"
PREDICTED_PDB=""  # User needs to provide this from AlphaFold

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --seq1)
            SEQ1="$2"
            shift 2
            ;;
        --seq2)
            SEQ2="$2"
            shift 2
            ;;
        --pdb-id)
            PDB_ID="$2"
            shift 2
            ;;
        --predicted)
            PREDICTED_PDB="$2"
            shift 2
            ;;
        --output)
            OUTPUT_BASE="$2"
            shift 2
            ;;
        --chain1)
            CHAIN1="$2"
            shift 2
            ;;
        --chain2)
            CHAIN2="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

mkdir -p "$OUTPUT_BASE"

echo -e "${CYAN}╔══════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${CYAN}║                                                                  ║${NC}"
echo -e "${CYAN}║     COMPLETE PROTEIN-PROTEIN INTERACTION WORKFLOW                ║${NC}"
echo -e "${CYAN}║     From Prediction to Validation                                ║${NC}"
echo -e "${CYAN}║                                                                  ║${NC}"
echo -e "${CYAN}╚══════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# ============================================================================
# STEP 1: AlphaFold3 Prediction Setup
# ============================================================================
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BLUE}  STEP 1: AlphaFold3 Prediction Setup${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

bash 01_run_alphafold.sh "$SEQ1" "$SEQ2" "$OUTPUT_BASE/alphafold"

echo ""
echo -e "${YELLOW}⚠ MANUAL STEP REQUIRED:${NC}"
echo -e "  1. Run AlphaFold3 prediction using the instructions above"
echo -e "  2. Download the predicted structure (fold_model_0.pdb)"
echo -e "  3. Place it in: ${GREEN}$OUTPUT_BASE/alphafold/fold_model_0.pdb${NC}"
echo ""
echo -e "${YELLOW}Press Enter when you have the predicted structure ready...${NC}"

# If predicted PDB not provided, wait for user
if [ -z "$PREDICTED_PDB" ]; then
    read -p ""
    PREDICTED_PDB="$OUTPUT_BASE/alphafold/fold_model_0.pdb"
    
    if [ ! -f "$PREDICTED_PDB" ]; then
        echo -e "${RED}ERROR: Predicted structure not found: $PREDICTED_PDB${NC}"
        echo -e "${YELLOW}Continuing with remaining steps that don't require prediction...${NC}"
        SKIP_PREDICTION=true
    fi
else
    if [ ! -f "$PREDICTED_PDB" ]; then
        echo -e "${RED}ERROR: Predicted structure not found: $PREDICTED_PDB${NC}"
        SKIP_PREDICTION=true
    fi
fi

# ============================================================================
# STEP 2: Visualize with PyMOL
# ============================================================================
if [ "$SKIP_PREDICTION" != "true" ]; then
    echo ""
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}  STEP 2: PyMOL Visualization${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    
    python3 02_visualize_pymol.py "$PREDICTED_PDB" \
        --chain1 "$CHAIN1" \
        --chain2 "$CHAIN2" \
        --output "$OUTPUT_BASE/figures_predicted" \
        --script-only
    
    echo -e "${GREEN}✓ Visualization script created${NC}"
fi

# ============================================================================
# STEP 3: Arpeggio Analysis (Predicted)
# ============================================================================
if [ "$SKIP_PREDICTION" != "true" ]; then
    echo ""
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}  STEP 3: Arpeggio Analysis (Predicted Structure)${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    
    bash 03_run_arpeggio.sh "$PREDICTED_PDB" "$CHAIN1" "$CHAIN2" \
        "$OUTPUT_BASE/arpeggio_predicted"
    
    echo -e "${GREEN}✓ Arpeggio analysis complete (predicted)${NC}"
fi

# ============================================================================
# STEP 4: Analyze Interactions (Predicted)
# ============================================================================
if [ "$SKIP_PREDICTION" != "true" ]; then
    echo ""
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}  STEP 4: Analyze Interactions (Predicted)${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    
    JSON_PRED="$OUTPUT_BASE/arpeggio_predicted/$(basename ${PREDICTED_PDB%.pdb}).json"
    
    python3 04_analyze_interactions.py "$JSON_PRED" \
        --output "$OUTPUT_BASE/analysis_predicted" \
        --csv --report --hot-spots
    
    echo -e "${GREEN}✓ Interaction analysis complete (predicted)${NC}"
fi

# ============================================================================
# STEP 5: Download Experimental Structure
# ============================================================================
echo ""
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BLUE}  STEP 5: Download Experimental Structure${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

python3 05_download_pdb.py "$PDB_ID" "$OUTPUT_BASE/experimental" --info

EXPERIMENTAL_PDB="$OUTPUT_BASE/experimental/${PDB_ID,,}.pdb"
echo -e "${GREEN}✓ Experimental structure downloaded${NC}"

# ============================================================================
# STEP 6: Visualize Experimental
# ============================================================================
echo ""
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BLUE}  STEP 6: Visualize Experimental Structure${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

python3 02_visualize_pymol.py "$EXPERIMENTAL_PDB" \
    --chain1 "$CHAIN1" \
    --chain2 "$CHAIN2" \
    --output "$OUTPUT_BASE/figures_experimental" \
    --script-only

echo -e "${GREEN}✓ Visualization script created${NC}"

# ============================================================================
# STEP 7: Arpeggio Analysis (Experimental)
# ============================================================================
echo ""
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BLUE}  STEP 7: Arpeggio Analysis (Experimental Structure)${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

bash 03_run_arpeggio.sh "$EXPERIMENTAL_PDB" "$CHAIN1" "$CHAIN2" \
    "$OUTPUT_BASE/arpeggio_experimental"

echo -e "${GREEN}✓ Arpeggio analysis complete (experimental)${NC}"

# ============================================================================
# STEP 8: Analyze Interactions (Experimental)
# ============================================================================
echo ""
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BLUE}  STEP 8: Analyze Interactions (Experimental)${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

JSON_EXP="$OUTPUT_BASE/arpeggio_experimental/$(basename ${EXPERIMENTAL_PDB%.pdb}).json"

python3 04_analyze_interactions.py "$JSON_EXP" \
    --output "$OUTPUT_BASE/analysis_experimental" \
    --csv --report --hot-spots

echo -e "${GREEN}✓ Interaction analysis complete (experimental)${NC}"

# ============================================================================
# STEP 9: Compare Predicted vs Experimental
# ============================================================================
if [ "$SKIP_PREDICTION" != "true" ]; then
    echo ""
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}  STEP 9: Compare Predicted vs Experimental${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    
    python3 06_compare_structures.py \
        "$PREDICTED_PDB" \
        "$EXPERIMENTAL_PDB" \
        "$JSON_PRED" \
        "$JSON_EXP" \
        --output "$OUTPUT_BASE/comparison/validation"
    
    echo -e "${GREEN}✓ Comparison complete${NC}"
fi

# ============================================================================
# Summary
# ============================================================================
echo ""
echo -e "${CYAN}╔══════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${CYAN}║                                                                  ║${NC}"
echo -e "${CYAN}║                   WORKFLOW COMPLETE! ✓                           ║${NC}"
echo -e "${CYAN}║                                                                  ║${NC}"
echo -e "${CYAN}╚══════════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${GREEN}All results saved to:${NC} $OUTPUT_BASE/"
echo ""
echo -e "${YELLOW}Generated files:${NC}"
echo -e "  📂 alphafold/              - AlphaFold prediction files"
echo -e "  📂 experimental/           - PDB structure"
echo -e "  📂 figures_predicted/      - PyMOL visualizations (predicted)"
echo -e "  📂 figures_experimental/   - PyMOL visualizations (experimental)"
echo -e "  📂 arpeggio_predicted/     - Arpeggio analysis (predicted)"
echo -e "  📂 arpeggio_experimental/  - Arpeggio analysis (experimental)"
echo -e "  📂 analysis_predicted/     - Interaction analysis (predicted)"
echo -e "  📂 analysis_experimental/  - Interaction analysis (experimental)"
if [ "$SKIP_PREDICTION" != "true" ]; then
    echo -e "  📂 comparison/             - Validation report ⭐"
fi
echo ""

if [ "$SKIP_PREDICTION" != "true" ]; then
    echo -e "${BLUE}Key results to check:${NC}"
    echo -e "  1. Validation report: ${GREEN}$OUTPUT_BASE/comparison/validation_report.txt${NC}"
    echo -e "  2. Interaction CSV: ${GREEN}$OUTPUT_BASE/analysis_predicted_interactions.csv${NC}"
    echo -e "  3. PyMOL sessions: ${GREEN}$OUTPUT_BASE/figures_*/*.pse${NC}"
fi

echo ""
echo -e "${YELLOW}Next steps:${NC}"
echo -e "  • Review validation report for quality metrics"
echo -e "  • Visualize structures in PyMOL"
echo -e "  • Analyze hot spot residues"
echo -e "  • Interpret biological significance"
echo ""

echo -e "${MAGENTA}🎉 Analysis complete! Happy researching! 🧬${NC}"
echo ""
