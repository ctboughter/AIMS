# AIMS Immune - Complete API Documentation

## Overview

This directory contains comprehensive API documentation for the **AIMS Immune** package - a software suite for the analysis of immune repertoires. AIMS includes advanced tools for sequence encoding, biophysical analysis, clustering, and statistical characterization of immunological sequences.

---

## Quick Navigation

### Core Modules

1. **[aims_loader.md](aims_loader.md)** - Sequence Loading & I/O
   - Load FASTA and CSV sequence files
   - Auto-format detection
   - Duplicate removal
   - Sequence subsetting

2. **[aims_analysis.md](aims_analysis.md)** - Core Analysis Engine
   - Sequence encoding matrices
   - Information theory calculations (Shannon entropy, MI)
   - Biophysical property analysis
   - Statistical testing
   - MHC and TCR-specific functions

3. **[aims_classification.md](aims_classification.md)** - Machine Learning & Classification
   - Feature matrix generation
   - Dimensionality reduction
   - Classification algorithms (MDA, SVM, Random Forest)
   - Cross-validation strategies
   - Feature selection methods

4. **[aims_utils.md](aims_utils.md)** - Utility Functions for Notebooks
   - Data preparation pipelines
   - Visualization functions
   - Clustering and metadata management
   - Property analysis and statistics
   - Sequence output utilities

5. **[aims_cli.md](aims_cli.md)** - Command-Line Interface
   - Complete analysis pipeline arguments
   - Usage examples
   - Output file descriptions
   - Workflow documentation

---

## Architecture Overview

### Data Flow Pipeline

```
Input Files (FASTA/CSV)
    ↓
[aims_loader] → Load & preprocess sequences
    ↓
[aims_analysis] → Encode to matrices
    ↓
[aims_classification] → Generate features
    ↓
[aims_utils/aims_cli] → Reduce dimensions & cluster
    ↓
Visualization & Statistics
    ↓
Output (figures, data files)
```

### Module Dependencies

```
aims_loader          [No internal dependencies]
    ↓
aims_analysis        [Requires: aims_loader]
    ↓
aims_classification  [Requires: aims_analysis]
    ↓
aims_utils           [Requires: all above]
    ↓
aims_cli             [Requires: all above]
```

---

## Key Concepts

### Sequence Encoding

Sequences are converted to numeric matrices where:
- **0** = Gap/padding
- **1-20** = Amino acids (alphabetical order: A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V)

### Biophysical Properties

62 biophysical properties per amino acid:
- **0-1**: Custom encodings
- **2-15**: Traditional properties (Charge, Hydropathy, Bulk, Flexibility, etc.)
- **16-61**: Hot-spot prediction properties

### Matrix Dimensions

For typical immunoglobulin sequences:
- **6 CDR loops** (CDR1L, CDR2L, CDR3L, CDR1H, CDR2H, CDR3H)
- Variable lengths per loop (center-aligned)
- Total matrix: sequences × positions × properties

### Alignment Methods

- **center**: Center-align shorter sequences
- **left**: Left-align with right padding
- **right**: Right-align with left padding
- **bulge**: Germline-aligned with bulge centered (CDR3)

---

## Function Categories by Purpose

### Loading & Data Preparation
- `aims_loader.seq_loader()` - Load sequences
- `aims_loader.get_msa_sub()` - Extract subsequences
- `aims_utils.loadDat()` - Multi-file loading
- `aims_utils.get_keys()` - Define amino acid ordering

### Sequence Encoding
- `aims_analysis.gen_tcr_matrix()` - TCR/Ig encoding
- `aims_analysis.gen_1Chain_matrix()` - Single-chain encoding
- `aims_analysis.gen_peptide_matrix()` - Peptide encoding
- `aims_analysis.gen_MSA_matrix()` - MSA encoding with gaps

### Feature Generation
- `aims_classification.getBig()` - Biophysical properties
- `aims_classification.get_bigass_matrix()` - Full feature matrix
- `aims_utils.grab_big()` - Integrated feature pipeline

### Analysis
- `aims_analysis.calculate_shannon()` - Shannon entropy
- `aims_analysis.calculate_MI()` - Mutual information
- `aims_analysis.gen_dset_props()` - Property averages
- `aims_analysis.prop_patterning()` - Discriminative features

### Classification
- `aims_classification.do_classy_mda()` - Cross-validated classification
- `aims_classification.do_linear_split()` - Direct comparison
- `aims_classification.classy_apply()` - Apply trained model

### Statistics
- `aims_analysis.do_statistics()` - Permutation testing
- `aims_utils.netStats()` - Property significance
- `aims_utils.posAvg_stats()` - Position-wise testing

### Visualization
- `aims_utils.plot_clusters()` - Projection plots
- `aims_utils.plotSubsets()` - Sequence heatmaps
- `aims_utils.plot_quantClust()` - Composition bar plots
- `aims_utils.do_netAvg()` - Property bar plots
- `aims_utils.do_posAvg()` - Position-sensitive plots
- `aims_utils.pos_Shannon()` - Entropy plots
- `aims_utils.pos_MI()` - MI heatmaps

---

## Common Workflows

### Workflow 1: Simple Clustering Analysis

```python
from aims_immune import aims_loader, aims_analysis, aims_utils

# 1. Load sequences
seq = aims_loader.seq_loader('data.fasta', 'dataset')

# 2. Get matrix size
mat_size = aims_analysis.get_sequence_dimension(seq)

# 3. Encode sequences
encoded = aims_analysis.gen_tcr_matrix(seq.values, giveSize=mat_size)

# 4. Generate features
features = utils.grab_big(seq, seq.values, encoded, 
                         ['dataset'], normalize='msuv', nCores=-1)

# 5. Cluster
clust_input, clusters = utils.do_cluster(features[1], encoded, 'pca', 'kmean', NClust=5)

# 6. Visualize
metadata = utils.get_metadata(token_df, ['dataset'], clusters)
fig, ax = utils.plot_clusters(clust_input, metadata, 'kmean')
```

### Workflow 2: Binary Classification with Statistics

```python
# 1. Load two datasets
seq1 = aims_loader.seq_loader('healthy.csv', 'Healthy')
seq2 = aims_loader.seq_loader('disease.csv', 'Disease')

# 2. Encode and generate features (for both)
# ... [encoding steps]

# 3. Perform classification
accuracies = classy.do_classy_mda(seq1_matrix, seq2_matrix, 
                                   feat_sel='max_diff', classif='forest')

# 4. Get discriminative properties
full_big, weights, acc, _, _, top_names = classy.do_linear_split(
    seq1_matrix, seq2_matrix, ridCorr=True)

# 5. Test significance
p_vals = aims.do_statistics(data1, data2, test='average', num_reps=10000)
```

### Workflow 3: Command-Line Pipeline

```bash
aims-cli -m Ig \
  -f dataset1.csv dataset2.csv \
  -n "Healthy" "Disease" \
  -a bulge -nl 6 \
  -c kmean -cs 3 \
  -ds True -db True \
  -pp True \
  -od ./results
```

---

## Parameter Guide

### Common Parameter Meanings

| Parameter | Type | Common Values | Purpose |
|-----------|------|---------------|---------|
| `AA_key` | list | 20-char string | Amino acid ordering |
| `mat_size` | int/list | 20-100 | Feature matrix dimensions |
| `alignment` | str | center/left/right/bulge | Sequence alignment |
| `normalize` | str | msuv/zscore/0to1 | Property normalization |
| `feat_sel` | str | PCA/kbest/max_diff/none | Feature selection method |
| `classif` | str | mda/svm/logReg/forest | Classification algorithm |
| `xVal` | str | kfold/loo/strat_kfold | Cross-validation strategy |
| `drop_dups` | bool | True/False | Remove duplicate sequences |
| `ridCorr` | bool | True/False | Remove correlated features |
| `nCores` | int | -1/1/2/4... | Parallel threads (-1=all) |

---

## Data Format Specifications

### FASTA Format
```
>sequence_id description
MKVLWAALLVTFLAGCAKAKAVS...
>sequence_id2 description
MKVLWAALLVTFLAGCAKAKAVS...
```

### CSV Format (Ig/Peptide)
```
CDR1L,CDR2L,CDR3L,CDR1H,CDR2H,CDR3H
AARY,ATYYDY,CARXXXXX,ARSY,GSYY,CARQ
...
```

### CSV Format (MSA)
```
Sequence
MKVLWAALLVTFLAGCAKAKAVS...
MKVLWAALLVTFLAGCAKAKAVS...
...
```

---

## Installation & Setup

### Requirements
```
Python ≥ 3.7
numpy ≥ 1.18.1
pandas ≥ 1.0.3
scipy ≥ 1.4.1
scikit-learn ≥ 0.22.1
matplotlib ≥ 3.1.3
biopython ≥ 1.76
umap-learn ≥ 0.5.3
numba (for parallel processing)
```

### Installation
```bash
pip install aims_immune
# OR
pip install -e . # (from source directory)
```

---

## Performance Notes

### Computational Complexity

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| Sequence loading | O(n) | n = sequences |
| Matrix encoding | O(n×p) | p = positions |
| Feature generation | O(n×p×64) | 64 properties |
| Distance calculation | O(n²) | Scales quadratically |
| Clustering | O(n log n) to O(n²) | Algorithm-dependent |
| MI calculation | O(n×p²) | Most expensive |

### Optimization Tips

- Use `nCores=-1` for parallel processing
- Use `cd='parse'` (correlation-filtered) instead of `cd='full'`
- For large datasets: use UMAP + DBSCAN instead of PCA + KMeans
- Skip MI calculation for >10,000 sequences

---

## Citations & References

When publishing research using AIMS, please cite:

- **Core AIMS papers**: [Check GitHub repository for citations]
- **Biophysical properties**: Liu et al. "Hot spot prediction in protein-protein interactions by an ensemble system", BMC Systems Biology
- **Visualization methods**: Inspired by information-theoretic approaches in immunology

---

## Troubleshooting

### Common Issues

**"ImportError: No module named aims_immune"**
- Solution: Install package with `pip install -e .`

**"ValueError: could not convert string to float"**
- Solution: Check data format; ensure numeric columns for classifiers

**"MemoryError during feature generation"**
- Solution: Reduce dataset size or use `cd='parse'` to filter features

**"Mismatched dimensions in matrix operations"**
- Solution: Verify all sequences have same loop structure

---

## API Consistency

### Naming Conventions
- `mono_` prefix: monoreactive sequences
- `poly_` prefix: polyreactive sequences
- `pre_` prefix: preliminary/intermediate data
- `fin_` prefix: final output

### Return Value Patterns
- Single item: direct type
- Multiple items: tuple of types
- Optional items: included if parameter=True

### Error Handling
- Invalid parameters: raise descriptive exceptions
- Data mismatches: print warnings, attempt recovery
- Resource limits: raise MemoryError or similar

---

## Contact & Support

- **Documentation**: See individual module files
- **Issues/Bugs**: Check GitHub issues
- **Contributing**: See CONTRIBUTING.md in repository

---

## Version History

- **0.9.4**: Current stable release
- **0.9.x**: Ongoing development
- **0.8.x**: Previous stable versions (archived)

---

## License

[See LICENSE file in repository]

---

## Quick Reference Card

### Essential Imports
```python
from aims_immune import aims_loader as aimsLoad
from aims_immune import aims_analysis as aims
from aims_immune import aims_classification as classy
from aims_immune import aims_utils as utils
```

### Essential Functions
```python
# Load
seq = aimsLoad.seq_loader('file.fasta', 'label')

# Encode
mat = aims.gen_tcr_matrix(seq.values, giveSize=mat_size)

# Features
features = classy.get_bigass_matrix(seq.values)

# Cluster
clust_input, clusters = utils.do_cluster(features, 'pca', 'kmean')

# Plot
fig, ax = utils.plot_clusters(clust_input, metadata, 'kmean')
```

### Essential Parameters
```
alignment='center'     # center, left, right, bulge
feat_sel='max_diff'   # PCA, kPCA, kbest, max_diff, none
classif='mda'         # mda, svm, logReg, forest
normalize='msuv'      # msuv, zscore, 0to1
```

---

Last Updated: 2026-09-09
Documentation Version: 1.0
