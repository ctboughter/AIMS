# AIMS Immune - aims_cli Module API

## Module Overview
The `aims_cli` module provides a unified command-line interface for the complete AIMS analysis pipeline, allowing users to perform comprehensive immune repertoire analysis without writing Python code.

---

## Main Functions

### `main()`
Parses command-line arguments and returns argument namespace.

**Returns:**
- `args` (argparse.Namespace): Parsed command-line arguments

---

### `run()`
Executes the complete AIMS analysis pipeline with user-provided parameters.

**Workflow:**
1. Parse arguments
2. Load sequence data
3. Generate AIMS encoding matrix
4. Calculate biophysical property matrices
5. Perform dimensionality reduction
6. Cluster sequences
7. Incorporate metadata
8. Visualize results
9. Perform subset analysis
10. Calculate statistics (optional)

---

## Command-Line Arguments

### Data Input Arguments

```bash
-dd, --datDir DIRECTORY
    Data directory containing input files (default: './')

-od, --outputDir DIRECTORY
    Output directory for results (default: 'AIMS_out')

-m, --molecule {Ig,MSA,Peptide}
    Molecule type for analysis (required)

-f, --fileNames FILE1 FILE2 ...
    Input filenames (required; space-separated list)

-n, --datNames NAME1 NAME2 ...
    Dataset labels (default: auto-generated as dat0, dat1, etc.)
```

### Sequence Processing Arguments

```bash
-a, --align {center,left,right,bulge}
    Sequence alignment method for Ig (default: 'center')

-nl, --numLoop INT
    Number of loops for Ig analysis (default: 1)

-dp, --dropDup {True,False}
    Drop duplicate sequences (default: False)

-s, --subset {True,False}
    Enable subset extraction (default: False)

-ss, --subStart INT INT ...
    Start positions for subset extraction

-se, --subEnd INT INT ...
    End positions for subset extraction

-bp, --bulgePad {2,4,6,8}
    Bulge region padding for bulge alignment (default: 8)
```

### Feature Engineering Arguments

```bash
-np, --normProp {msuv,zscore,0to1}
    Property normalization method (default: 'msuv')

-rn, --REnorm {True,False}
    Renormalize by Shannon entropy (default: True)

-cd, --clustData {full,parse,avg}
    Data format for clustering (default: 'parse')
    - 'full': All features
    - 'parse': Correlation-filtered
    - 'avg': Per-sequence property averages
```

### Dimensionality Reduction & Clustering

```bash
-pa, --projAlg {pca,umap}
    Projection algorithm (default: 'pca')

-us, --umapSeed INT
    Random seed for UMAP (default: unset = random)

-c, --clustAlg {kmean,optics,dbscan}
    Clustering algorithm (default: 'optics')

-cs, --clustSize INT
    Cluster size parameter:
    - KMeans: number of clusters (default: 10)
    - OPTICS/DBSCAN: minimum samples (default: 10)
```

### Visualization Arguments

```bash
-sp, --showProj {2d,3d,both}
    Projection dimensions to display (default: 'both')

-sc, --showClust {cluster,metadata,both}
    Visualization coloring scheme (default: 'both')

-nb, --normBar {True,False}
    Normalize bar plots to fractions (default: True)

-sf, --saveFmt {png,pdf}
    Figure save format (default: 'pdf')
```

### Analysis Selection Arguments

```bash
-as, --analysisSel {cluster,metadata}
    Subset type for detailed analysis (default: 'cluster')

-sd, --selDat INT INT ...
    Specific subsets to analyze (default: [0,1])

-lo, --seqlogo {True,False}
    Generate seqlogo plots (default: False)

-ln, --logoNum INT
    Seqlogo sequence length (default: 14)

-sv, --saveSeqs {True,False}
    Save clustered sequences separately (default: False)

-sac, --saveAllClusts {True,False}
    Save sequences for all clusters (default: False)
```

### Biophysical Property Analysis

```bash
-p1, --prop1 INT
    First property to analyze (1=Charge, 2=Hydrophobicity) (default: 1)

-p2, --prop2 INT
    Second property to analyze (default: 2)

-pp, --Plotprops {True,False}
    Plot biophysical properties for clusters (default: False)
```

### Statistical Analysis Arguments

```bash
-ds, --DOstats {True,False}
    Calculate statistical significance (default: False)

-db, --DOboot {True,False}
    Use bootstrap resampling (default: False)

-bt, --boots INT
    Bootstrap replicates (default: 1000)

-mi, --MIboots INT
    Mutual information bootstrap replicates (default: 10)
```

### Distance-Based Analysis

```bash
-gd, --GETdist {True,False}
    Calculate sequence distances (AIMS distance) (default: False)

-pd, --PARdist {True,False}
    Parallel distance calculation (default: False)
```

### Metadata & Miscellaneous

```bash
-mf, --metaForm {category,quant}
    Metadata format (default: 'category')

-mn, --metaName STRING
    Metadata column name (default: 'meta')

-aa, --AAorder STRING
    Custom amino acid order (20-letter string; default: standard)

-cc, --colors COLOR1 COLOR2 ...
    Colors for plots (default: ['purple','orange'])

-p, --parallel {True,False}
    Use parallel processing (default: False)

-ms, --matSize INT
    Feature matrix size for LDA (default: 10)

-sl, --showLabel {True,False}
    Annotate points with labels (default: False)

-co, --clustOnly {True,False}
    Stop after clustering (default: False)
```

---

## Output Files

### Main Figures
- `AIMS_mat.png/pdf` - Encoded sequence matrix visualization
- `AIMS_projections.png` - 2D/3D dimensionality reduction plots
- `AIMS_clusterQuant.png/pdf` - Stacked bar chart of cluster composition
- `AIMS_clusterPurity.png/pdf` - Cluster purity and significance heatmap
- `AIMS_[subset_sel]_subViz.png/pdf` - Selected subset visualization
- `AIMS_posSensAvg.png/pdf` - Position-sensitive property averages
- `AIMS_netAvgProp.png/pdf` - Network-averaged properties
- `AIMS_entropy.png/pdf` - Shannon entropy by position
- `AIMS_MI.png/pdf` - Mutual information matrices
- `AIMS_MIdiff.png/pdf` - Differential MI between subsets
- `AIMS_freq_sides.png/pdf` - Amino acid frequency heatmaps

### Cluster Characterization
- `AIMS_*_charge.png/pdf` - Charge per sequence/position
- `AIMS_*_hydropathy.png/pdf` - Hydropathy per sequence/position
- `cluster[#]_all.txt` - Sequences in cluster (if --saveAllClusts)

### Statistical Results
- `bar_sig.csv` - P-values for property bar plots
- `AIMS_PosSens_pval.png/pdf` - Position-sensitive p-values
- `aims_entropy_pval.png/pdf` - Entropy p-values
- `dist_clust.csv` - Distance-based cluster assignments

### Sequence Files
- `[label]_all.txt` - All sequences in subset
- `[label]_logo.txt` - Sequences of specified length

---

## Common Usage Examples

### Immunoglobulin (Ig) Analysis
```bash
aims-cli -m Ig \
  -f human_antibody.csv mouse_antibody.csv \
  -n Human Mouse \
  -nl 6 \
  -a bulge \
  -dd ./data \
  -od ./output \
  -c optics \
  -cs 10
```

### Multiple Sequence Alignment (MSA)
```bash
aims-cli -m MSA \
  -f hla_alleles.fasta \
  -n "HLA-A" \
  -dd ./data \
  -od ./msa_results \
  -c kmean \
  -cs 5
```

### Peptide Analysis
```bash
aims-cli -m Peptide \
  -f peptides.csv \
  -n "MHC_bound" \
  -cd avg \
  -c hdbscan \
  -ds True \
  -db True \
  -bt 10000
```

### Statistical Comparison
```bash
aims-cli -m Ig \
  -f dataset1.csv dataset2.csv \
  -n "Healthy" "Disease" \
  -c kmean -cs 3 \
  -sd 0 1 \
  -ds True \
  -db True \
  -bt 1000 \
  -pp True
```

### Distance-Based Analysis
```bash
aims-cli -m Ig \
  -f sequences.fasta \
  -n "TCR_repertoire" \
  -c kmean -cs 5 \
  -gd True \
  -pd True \
  -cd full
```

---

## Workflow Steps

### 1. Data Loading
- Reads sequences from specified files
- Auto-detects FASTA vs CSV format
- Removes duplicates if requested
- Extracts subsequences (optional)

### 2. Encoding
- Converts amino acid sequences to numeric encoding
- Applies selected alignment method
- Generates position-specific matrices

### 3. Feature Generation
- Calculates biophysical properties for all positions
- Removes highly correlated features (>0.75)
- Removes all-zero features
- Normalizes properties

### 4. Dimensionality Reduction
- Applies PCA or UMAP
- Generates 2D and/or 3D projections

### 5. Clustering
- Applies selected algorithm
- Returns cluster assignments

### 6. Visualization
- Creates projection plots
- Generates cluster composition charts
- Produces per-subset heatmaps

### 7. Characterization
- Calculates Shannon entropy
- Computes mutual information
- Analyzes biophysical properties

### 8. Statistics (Optional)
- Permutation tests on properties
- P-value calculations
- Bootstrap confidence intervals

---

## Important Notes

### Data Format Requirements
- **FASTA**: Standard format with >header lines
- **CSV**: Header row required; one sequence per row
- **Columns**: For Ig/Peptide: one per loop; for MSA: single sequence per row

### Molecule Types
- **Ig**: Immunoglobulin (antibodies) or TCR; multiple CDR loops
- **MSA**: Multiple sequence alignment; allows gaps
- **Peptide**: Single sequences; 8-20 amino acids

### Alignment Methods
- **center**: Equal padding on both sides
- **left**: Pad on right only
- **right**: Pad on left only
- **bulge**: Germline regions aligned, CDR3 bulge centered

### Performance Considerations
- Parallel processing (-p True) recommended for large datasets
- UMAP slower but often better for immune data than PCA
- HDBSCAN requires scikit-learn ≥1.3
- Distance calculations scale as O(n²)

### Output Organization
- All figures saved in --outputDir
- Named with prefixes indicating analysis type
- PNG format: high resolution (600 dpi)
- PDF format: vector graphics, smaller file size

---

## Error Handling

Common error messages and solutions:

```
"ERROR: Custom Key Wrong Length!"
- Solution: --AAorder must be exactly 20 characters

"Wrong number of loops selected"
- Solution: Verify --numLoop matches sequence structure

"Wrong # loops, go back and redefine"
- Solution: CSV missing required loop columns; check format
```

---

## Environment Setup

Required environment variables (optional):
```bash
# Limit OpenMP threads (useful for multi-agent runs)
export OMP_NUM_THREADS=2

# Set matplotlib backend if running headless
export MPLBACKEND=Agg
```

---

## Reproducing Analysis

To ensure reproducible results:
```bash
# Use fixed random seeds
-us 42  # UMAP seed
# For other algorithms, set in configuration

# Note: Some randomness may remain in some algorithms
# Distance-based clustering shows good reproducibility
```

---

## Exit Codes

- `0`: Success
- `1`: Argument parsing error
- `Other`: Runtime error (check console output)
