# Quick Start Guide

Get started with the complete protein-protein interaction workflow in 10 minutes!

---

##  Installation (5 minutes)

### 1. Prerequisites

```bash
# Python 3.8+
python --version

# Create conda environment (recommended)
conda create -n protein-workflow python=3.9
conda activate protein-workflow
```

### 2. Install Arpeggio

```bash
git clone https://github.com/PDBeurope/arpeggio.git
cd arpeggio
python setup.py install
cd ..
```

### 3. Install PyMOL (optional, for visualization)

```bash
conda install -c conda-forge pymol-open-source
```

### 4. Install Python packages

```bash
cd protein-interaction-workflow
pip install -r requirements.txt
```

---

##  Run Example: MDM2/p53 (5 minutes)

### Option A: Step-by-Step (Learn the workflow)

```bash
cd scripts

# STEP 1: Prepare AlphaFold input
bash 01_run_alphafold.sh ../examples/mdm2_p53/sequences/mdm2.fasta \
                         ../examples/mdm2_p53/sequences/p53.fasta \
                         ../examples/mdm2_p53/results/alphafold

# → Follow instructions to run AlphaFold3
# → Download prediction to ../examples/mdm2_p53/results/alphafold/fold_model_0.pdb

# STEP 2: Visualize (generates PyMOL script)
python 02_visualize_pymol.py ../examples/mdm2_p53/results/alphafold/fold_model_0.pdb \
    --chain1 A --chain2 B \
    --output ../examples/mdm2_p53/results/figures \
    --script-only

# STEP 3: Run Arpeggio
bash 03_run_arpeggio.sh ../examples/mdm2_p53/results/alphafold/fold_model_0.pdb A B \
    ../examples/mdm2_p53/results/arpeggio_predicted

# STEP 4: Analyze interactions
python 04_analyze_interactions.py \
    ../examples/mdm2_p53/results/arpeggio_predicted/fold_model_0.json \
    --output ../examples/mdm2_p53/results/analysis \
    --csv --report --hot-spots

# STEP 5: Download experimental structure
python 05_download_pdb.py 1YCR ../examples/mdm2_p53/results/experimental

# STEP 6: Repeat Arpeggio on experimental
bash 03_run_arpeggio.sh ../examples/mdm2_p53/results/experimental/1ycr.pdb A B \
    ../examples/mdm2_p53/results/arpeggio_experimental

# STEP 7: Compare predicted vs experimental
python 06_compare_structures.py \
    ../examples/mdm2_p53/results/alphafold/fold_model_0.pdb \
    ../examples/mdm2_p53/results/experimental/1ycr.pdb \
    ../examples/mdm2_p53/results/arpeggio_predicted/fold_model_0.json \
    ../examples/mdm2_p53/results/arpeggio_experimental/1ycr.json \
    --output ../examples/mdm2_p53/results/comparison/validation
```

### Option B: Automated (Run entire workflow)

```bash
cd scripts
bash run_complete_workflow.sh
```

**Note:** You'll need to manually run AlphaFold and provide the predicted structure.

---

##  BONUS: Protein-Ligand Example

Analyze MDM2-Nutlin (anticancer drug) interactions:

```bash
cd scripts

# Download structure with ligand
python 05_download_pdb.py 4HG7 ligand_example/

# Analyze protein-ligand interactions
bash 07_protein_ligand_example.sh ligand_example/4hg7.pdb NUX \
    ligand_example/arpeggio

# Visualize in PyMOL
pymol ligand_example/arpeggio/visualize_ligand.pml
```

---

##  Expected Results (MDM2/p53)

### Structural Quality
- **C-α RMSD**: ~0.85 Å (excellent)
- **Interface RMSD**: ~0.62 Å (outstanding)

### Interaction Analysis
- **Total interactions**: ~28 predicted, ~26 experimental
- **F1-score**: ~0.85 (excellent agreement)

### Key Residues
- **Phe19 (p53)**: Hydrophobic anchor
- **Trp23 (p53)**: π-π stacking with Tyr67
- **Leu26 (p53)**: Hydrophobic core

---

##  Output Structure

```
results/
├── alphafold/
│   └── fold_model_0.pdb           # Predicted structure
├── experimental/
│   └── 1ycr.pdb                   # Experimental structure
├── figures_predicted/
│   ├── overview.png               # Structure overview
│   ├── interface.png              # Interface close-up
│   └── visualize.pml              # PyMOL script
├── arpeggio_predicted/
│   ├── fold_model_0.json          # Interaction data
│   └── fold_model_0.contacts      # Contact list
├── arpeggio_experimental/
│   ├── 1ycr.json
│   └── 1ycr.contacts
├── analysis/
│   ├── analysis_interactions.csv  # All interactions
│   └── analysis_report.txt        # Summary statistics
└── comparison/
    └── validation_report.txt      # ⭐ Main validation report
```

---

## 🔧 Troubleshooting

### "Arpeggio not found"
```bash
# Check installation
which arpeggio

# If not found, reinstall
cd arpeggio
python setup.py install --user
```

### "BioPython ImportError"
```bash
pip install biopython
```

### "PyMOL not available"
```bash
# Option 1: Conda
conda install -c conda-forge pymol-open-source

# Option 2: Use --script-only flag
python 02_visualize_pymol.py structure.pdb --script-only
# Then open script in PyMOL GUI
```

### "AlphaFold prediction failed"
- Check sequence length (not too long)
- Verify FASTA format
- Try AlphaFold Server: https://alphafoldserver.com/
- Alternative: Use ColabFold

---

##  Using Your Own Proteins

### 1. Prepare Sequences

```bash
# Create FASTA files
echo ">ProteinA" > proteinA.fasta
echo "MKTAYIAKQRQISFVKSHFSRQLE..." >> proteinA.fasta

echo ">ProteinB" > proteinB.fasta
echo "PQITLWQRPLVTIKIGGQLKE..." >> proteinB.fasta
```

### 2. Run Workflow

```bash
cd scripts
bash run_complete_workflow.sh \
    --seq1 ../proteinA.fasta \
    --seq2 ../proteinB.fasta \
    --output ../my_results
```

### 3. If No Experimental Structure Exists

Skip comparison steps:
- Run steps 1-4 only
- Focus on AlphaFold confidence metrics (pLDDT, PAE)
- Compare multiple AlphaFold models (models 0-4)

---

##  Understanding the Results

### Structural Metrics

**RMSD (Root Mean Square Deviation)**
- < 1.0 Å: Excellent
- 1-2 Å: Good
- 2-3 Å: Moderate
- > 3 Å: Poor

**TM-score**
- > 0.9: Nearly identical
- 0.5-0.9: Same fold
- < 0.5: Different fold

### Interaction Metrics

**Precision**: Are predicted interactions real?
**Recall**: Did we find all real interactions?
**F1-score**: Overall agreement
- > 0.8: Excellent
- 0.6-0.8: Good
- 0.4-0.6: Moderate
- < 0.4: Poor

### AlphaFold Confidence

**pLDDT (per-residue confidence)**
- > 90: Very high confidence (blue)
- 70-90: Good confidence (cyan)
- 50-70: Low confidence (yellow)
- < 50: Very low confidence (orange/red)

**PAE (Predicted Aligned Error)**
- < 5 Å: Excellent inter-chain prediction
- 5-10 Å: Good
- > 10 Å: Uncertain

---

##  Next Steps

### For Analysis
1. Review validation report
2. Identify hot spot residues
3. Compare with literature
4. Design experiments

### For Visualization
```bash
# Open PyMOL session
pymol results/figures_predicted/*_session.pse

# Create custom views
# Edit visualize.pml script
```

### For Publication
- Use high-resolution images (2400x2400 DPI)
- Export interaction tables to Excel
- Include validation metrics
- Cite Arpeggio and AlphaFold3

---

##  Citations

**Arpeggio:**
```
Jubb et al. (2017). Arpeggio: A Web Server for Calculating and 
Visualising Interatomic Interactions in Protein Structures. 
J Mol Biol, 429(3):365-371.
```

**AlphaFold3:**
```
Abramson et al. (2024). Accurate structure prediction of biomolecular 
interactions with AlphaFold 3. Nature, 630:493-500.
```

**PyMOL:**
```
The PyMOL Molecular Graphics System, Version 2.0 Schrödinger, LLC.
```

---

## 🆘 Getting Help

- **Issues**: Check troubleshooting section above
- **Arpeggio docs**: https://github.com/PDBeurope/arpeggio
- **AlphaFold**: https://alphafold.ebi.ac.uk/
- **PyMOL wiki**: https://pymolwiki.org/

---

**Ready to analyze protein interactions!** 🧬🔬

For detailed documentation, see: [README.md](README.md)
