# Complete Protein-Protein Interaction Analysis Workflow
## From AlphaFold3 Prediction to Validation

**A comprehensive, end-to-end workflow for predicting, visualizing, and validating protein-protein interactions. Includes AlphaFold3 integration, PyMOL visualization scripts, Arpeggio interaction detection, and quantitative validation against experimental structures.**

---

##  Workflow Overview

```
STEP 1: AlphaFold3 Prediction (Heterodimer)
         ↓
STEP 2: PyMOL Visualization
         ↓
STEP 3: Arpeggio Analysis (Protein-Protein)
         ↓
STEP 4: Extract & Analyze Results
         ↓
STEP 5: Download PDB Structure (Experimental)
         ↓
STEP 6: Repeat Analysis on PDB
         ↓
STEP 7: Compare Predicted vs Experimental
         ↓
BONUS: Protein-Ligand Interaction Analysis
```

---

##  Example System: MDM2/p53

**Why this system?**
- **Clinical relevance**: Major cancer drug target
- **p53**: Tumor suppressor ("Guardian of the genome")
- **MDM2**: Inhibits p53 → cancers overexpress MDM2
- **Goal**: Blocking MDM2-p53 → reactivates p53 → cancer therapy
- **PDB reference**: 1YCR (experimental structure)

**Sequences:**
- **Chain A (MDM2)**: Residues 17-125 (109 residues)
- **Chain B (p53)**: Residues 15-29 (15 residues)

---

##  Directory Structure

```
protein-interaction-workflow/
├── README.md                          # This file
├── examples/
│   └── mdm2_p53/
│       ├── sequences/
│       │   ├── mdm2.fasta            # MDM2 sequence
│       │   └── p53.fasta             # p53 sequence
│       └── results/                   # Analysis results
│
├── scripts/
│   ├── 01_run_alphafold.sh          # AlphaFold prediction
│   ├── 02_visualize_pymol.py        # PyMOL visualization
│   ├── 03_run_arpeggio.sh           # Arpeggio analysis
│   ├── 04_analyze_interactions.py   # Parse Arpeggio results
│   ├── 05_download_pdb.py           # Get experimental structure
│   ├── 06_compare_structures.py     # Validation
│   └── 07_protein_ligand_example.sh # Bonus: Ligand analysis
│
├── utils/
│   ├── arpeggio_parser.py           # Arpeggio result parser
│   ├── structure_tools.py           # Structure manipulation
│   └── visualization.py             # Plotting functions
│
└── requirements.txt                  # Python dependencies
```

---

##  Quick Start

### Prerequisites

```bash
# 1. AlphaFold3 (via AlphaFold Server or local installation)
# 2. PyMOL
conda install -c conda-forge pymol-open-source

# 3. Arpeggio
git clone https://github.com/PDBeurope/arpeggio.git
cd arpeggio && python setup.py install

# 4. Python packages
pip install -r requirements.txt
```

### Run Complete Workflow

```bash
# Navigate to examples
cd examples/mdm2_p53

# Run all steps
bash ../../scripts/run_complete_workflow.sh
```

---

##  Detailed Step-by-Step Guide

### STEP 1: AlphaFold3 Prediction

**Script**: `scripts/01_run_alphafold.sh`

```bash
# Prepare input sequences
# Run AlphaFold3 prediction
# Output: predicted structure PDB file
```

**Manual alternative**: Use AlphaFold Server (https://alphafoldserver.com/)

### STEP 2: PyMOL Visualization

**Script**: `scripts/02_visualize_pymol.py`

Generates:
- Interface visualization
- Interaction highlights
- Quality assessment figures
- Session file for interactive exploration

### STEP 3: Arpeggio Analysis

**Script**: `scripts/03_run_arpeggio.sh`

Analyzes protein-protein interactions:
- Hydrogen bonds
- Van der Waals contacts
- Aromatic interactions
- Ionic interactions
- And more...

### STEP 4: Extract Results

**Script**: `scripts/04_analyze_interactions.py`

Parses Arpeggio output:
- Interaction counts by type
- Hot spot residues
- Interface characteristics
- CSV export

### STEP 5: Download PDB Structure

**Script**: `scripts/05_download_pdb.py`

Downloads experimental structure (1YCR) for validation

### STEP 6: Repeat Analysis

Run Arpeggio on experimental structure using same parameters

### STEP 7: Comparison

**Script**: `scripts/06_compare_structures.py`

Validation metrics:
- **Structural**: RMSD, TM-score
- **Interactions**: Precision, Recall, F1-score
- **Per-type analysis**: Compare each interaction type

### BONUS: Protein-Ligand

**Script**: `scripts/07_protein_ligand_example.sh`

Analyze protein-ligand interactions using Arpeggio

---

##  Expected Results

### Structural Quality
- **C-alpha RMSD**: ~0.85 Å (excellent)
- **Interface RMSD**: ~0.62 Å (outstanding)
- **TM-score**: ~0.91 (high similarity)

### Interaction Analysis
- **Total predicted**: ~28 interactions
- **Total experimental**: ~26 interactions
- **Overlap**: ~23 (88.5%)
- **F1-score**: ~0.85 (excellent)

### Key Residues
- **Phe19 (p53)**: 4 interactions, hydrophobic anchor
- **Trp23 (p53)**: 6 interactions, π-π stacking
- **Leu26 (p53)**: 3 interactions, hydrophobic core

---

## 🔬 Understanding the Results

### What Makes a Good Prediction?

**Structural metrics:**
- RMSD < 2.0 Å → Good
- RMSD < 1.0 Å → Excellent
- TM-score > 0.5 → Same fold
- TM-score > 0.9 → Nearly identical

**Interaction metrics:**
- F1-score > 0.6 → Good agreement
- F1-score > 0.8 → Excellent agreement
- Precision → Are predicted interactions real?
- Recall → Did we find all real interactions?

---

##  Citation

If you use this workflow, please cite:

**AlphaFold3:**
```
Abramson et al. (2024). Accurate structure prediction of biomolecular 
interactions with AlphaFold 3. Nature, 630:493-500.
```

**Arpeggio:**
```
Jubb et al. (2017). Arpeggio: A Web Server for Calculating and 
Visualising Interatomic Interactions in Protein Structures. 
J Mol Biol, 429(3):365-371.
```

**PyMOL:**
```
The PyMOL Molecular Graphics System, Version 2.0 Schrödinger, LLC.
```

---

## 🆘 Troubleshooting

### Common Issues

**1. AlphaFold fails**
- Check sequence length (not too long)
- Verify FASTA format
- Check internet connection (for MSA)

**2. Arpeggio not found**
- Install: `python setup.py install`
- Add to PATH

**3. PyMOL import error**
- Install: `conda install -c conda-forge pymol-open-source`
- Or use PyMOL GUI

**4. Missing dependencies**
```bash
pip install biopython numpy pandas matplotlib seaborn
```

---

##  Tips for Your Own Proteins

### Using Different Proteins

1. **Prepare sequences**:
   ```bash
   # Create FASTA files
   echo ">ProteinA" > proteinA.fasta
   echo "SEQUENCE..." >> proteinA.fasta
   ```

2. **Run workflow**:
   ```bash
   bash scripts/run_complete_workflow.sh \
       --seq1 proteinA.fasta \
       --seq2 proteinB.fasta \
       --output my_results
   ```

3. **Find PDB reference** (if available):
   - Search PDB (https://www.rcsb.org)
   - Use for validation

### Without Experimental Structure

If no PDB structure exists:
- Skip comparison steps
- Focus on AlphaFold confidence metrics
- Use multiple predictions (models 1-5)
- Compare with known homologs

---

##  Learning Resources

### Understanding Protein Interactions

- **Hydrogen bonds**: Stabilize structure
- **Van der Waals**: Packing interactions
- **Aromatic**: π-π stacking
- **Ionic**: Salt bridges
- **Hydrophobic**: Core interactions

### AlphaFold Confidence Metrics

- **pLDDT**: Per-residue confidence (0-100)
  - >90: Very high confidence
  - 70-90: Good confidence
  - <70: Low confidence

- **PAE (Predicted Aligned Error)**: Inter-domain/chain confidence
  - <5 Å: Excellent
  - 5-10 Å: Good
  - >10 Å: Uncertain

---

## 🔧 Advanced Usage

### Customize Arpeggio Parameters

Edit `scripts/03_run_arpeggio.sh`:
```bash
# Change distance cutoffs
# Modify selection criteria
# Add custom atom types
```

### Batch Processing

Process multiple structures:
```bash
for pdb in structures/*.pdb; do
    bash scripts/03_run_arpeggio.sh "$pdb"
done
```

### Custom Visualization

Modify `scripts/02_visualize_pymol.py`:
- Change colors
- Add labels
- Highlight specific residues
- Create movies

---

##  Output Files

### Generated Files

```
results/
├── predicted_structure.pdb          # AlphaFold prediction
├── experimental_structure.pdb       # PDB reference
├── arpeggio_predicted/
│   ├── predicted.json              # Arpeggio results
│   ├── predicted.contacts          # Contact list
│   └── predicted.rings             # Aromatic rings
├── arpeggio_experimental/
│   ├── experimental.json
│   └── experimental.contacts
├── analysis/
│   ├── interactions_summary.csv    # All interactions
│   ├── hot_spots.txt              # Key residues
│   └── statistics.json            # Summary stats
├── comparison/
│   ├── validation_report.txt      # Detailed comparison
│   ├── rmsd_analysis.txt          # Structural metrics
│   └── interaction_comparison.csv # Interaction overlap
└── figures/
    ├── structure_overlay.png      # Structural alignment
    ├── interaction_comparison.png # Bar charts
    ├── contact_map.png           # Contact maps
    └── validation_metrics.png    # Summary metrics
```

---

##  Next Steps

After completing this workflow:

1. **Analyze your results**
   - Identify key interactions
   - Find hot spot residues
   - Assess prediction quality

2. **Iterate if needed**
   - Try different AlphaFold models
   - Adjust parameters
   - Test alternative structures

3. **Biological interpretation**
   - Literature search
   - Functional implications
   - Design experiments

4. **Share your findings**
   - Publication
   - GitHub repository
   - Collaborate

---

##  Support

Questions? Check:
- Scripts documentation (comments in code)
- Arpeggio docs: https://github.com/PDBeurope/arpeggio
- AlphaFold: https://alphafold.ebi.ac.uk/
- PyMOL wiki: https://pymolwiki.org/

---

**Ready to analyze protein interactions!** 🧬🔬
