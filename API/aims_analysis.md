# AIMS Immune - aims_analysis Module API

## Module Overview
The `aims_analysis` module is the core analysis engine for immune repertoire characterization, providing sequence encoding, matrix generation, entropy calculations, and statistical analysis functions.

---

## Sequence Dimension & Matrix Functions

### `get_sequence_dimension(re_poly)`
Calculates maximum sequence length per loop/region.

**Parameters:**
- `re_poly` (pandas.DataFrame): Sequences (rows=loops, cols=sequences)

**Returns:**
- `seqlenF` (list): Maximum sequence length per loop (+3 spacing buffer)

---

### `gen_tcr_matrix(pre_poly, AA_key=AA_key, key=AA_num_key_new, binary=False, pre_mono=[], giveSize='', return_Size=False, manuscript_arrange=False, alignment='center', bulge_pad=8)`
Generates encoded sequence matrix for TCR/immunoglobulin sequences.

**Parameters:**
- `pre_poly` (numpy.ndarray/list): Polyreactive sequences
- `AA_key` (list): Amino acid alphabet
- `key` (numpy.ndarray): Numerical encoding for amino acids
- `binary` (bool): If True, also process pre_mono for comparison
- `pre_mono` (numpy.ndarray): Monoreactive sequences (if binary=True)
- `giveSize` (list/int): Explicit matrix dimensions; calculated if empty
- `return_Size` (bool): Also return calculated dimensions
- `manuscript_arrange` (bool): Rearrange loops [0,1,2,5,4,3]
- `alignment` (str): 'center'/'left'/'right'/'bulge'
- `bulge_pad` (int): Bulge region size (2, 4, 6, or 8)

**Returns:**
- `final_poly` (numpy.ndarray): Encoded polyreactive matrix
- `final_mono` (numpy.ndarray): Encoded monoreactive (if binary=True)
- `max_lenp` (list): Dimensions (if return_Size=True)

**Alignment Methods:**
- `'center'`: Center-align shorter sequences
- `'left'`: Left-align sequences
- `'right'`: Right-align sequences
- `'bulge'`: Center-align CDR3 bulge, cap ends

---

### `gen_1Chain_matrix(pre_poly, AA_key=AA_key, key=AA_num_key_new, binary=False, pre_mono=[], giveSize='', return_Size=False)`
Generates matrix for single-chain molecules (HLA, single domain antibodies).

**Parameters:**
- Same as gen_tcr_matrix (3 loops instead of 6)

**Returns:**
- Single-chain encoded matrices

---

### `gen_peptide_matrix(pre_pep1, AA_key=AA_key, key=AA_num_key_new, binary=False, pre_pep2=[])`
Generates matrix for peptide sequences.

**Parameters:**
- `pre_pep1` (numpy.ndarray): Peptide set 1
- `AA_key` (list): Amino acid alphabet
- `key` (numpy.ndarray): Numerical encoding
- `binary` (bool): Compare two peptide sets
- `pre_pep2` (numpy.ndarray): Peptide set 2 (if binary=True)

**Returns:**
- Encoded peptide matrices (14 positions fixed)

**Features:**
- Positions 2 and last treated as anchors
- Central region ("bulge") center-aligned
- Bulge size = total_length - 8

---

### `gen_MSA_matrix(pre_poly, AA_key_dash=AA_key_dash, key=AA_num_key_new, binary=False, pre_mono=[], giveSize='', return_Size=False)`
Generates matrix for multiple sequence alignments (allows gaps).

**Parameters:**
- `pre_poly` (numpy.ndarray): MSA sequences
- `AA_key_dash` (list): Amino acids + gap character
- `key` (numpy.ndarray): Numerical encoding (21 values)
- `binary` (bool): Compare two MSAs
- `pre_mono` (numpy.ndarray): Second MSA
- `giveSize` (list/int): Matrix dimensions
- `return_Size` (bool): Return dimensions

**Returns:**
- MSA-encoded matrices with gap handling

---

## Information Theory Functions

### `calculate_shannon(poly_PCA, num_threads=-1)`
Numba-accelerated Shannon entropy calculation.

**Parameters:**
- `poly_PCA` (numpy.ndarray): Encoded sequences (clones × positions)
- `num_threads` (int): Parallel threads (-1 = all)

**Returns:**
- `shannon_poly` (numpy.ndarray): Shannon entropy per position
- `poly_count` (numpy.ndarray): AA frequency per position
- `coverage` (numpy.ndarray): Coverage (non-zero) per position

**Entropy Calculation:**
- H = -Σ(p_i × log₂(p_i)) for each amino acid
- Typically ranges 0-4.32 bits (20 amino acids)

---

### `calculate_MI(poly_PCA)`
Numba-accelerated mutual information calculation.

**Parameters:**
- `poly_PCA` (numpy.ndarray): Encoded sequences

**Returns:**
- `MI_final_poly` (numpy.ndarray): Mutual information matrix (positions × positions)
- `poly_count_cond` (numpy.ndarray): Conditional probabilities
- `poly_count` (numpy.ndarray): Marginal probabilities

**Formula:**
- MI(X,Y) = H(X) - Σ P(X|Y)H(X|Y)
- Quantifies covariation between sequence positions

---

### `joint_prob(poly_PCA)`
Calculates joint probability distributions between positions.

**Parameters:**
- `poly_PCA` (numpy.ndarray): Encoded sequences

**Returns:**
- `joint` (numpy.ndarray): Joint probabilities (positions × AAs × positions × AAs)

---

## Biophysical Property Analysis

### `gen_dset_props(poly_PCA, props=properties, stdev=False)`
Averages biophysical properties across sequences.

**Parameters:**
- `poly_PCA` (numpy.ndarray): Encoded sequences
- `props` (numpy.ndarray): Property matrix
- `stdev` (bool): Also calculate standard deviations

**Returns:**
- `poly_prop` (numpy.ndarray): Properties × positions
- `poly_prop_stdev` (numpy.ndarray): Standard deviations (if stdev=True)

---

### `gen_clone_props(poly_PCA)`
Calculates biophysical properties averaged per sequence.

**Parameters:**
- `poly_PCA` (numpy.ndarray): Encoded sequences (AAs × positions)

**Returns:**
- `poly_prop_pca` (numpy.ndarray): Properties × sequences

---

### `get_props()`
Retrieves amino acid keys and property matrices.

**Returns:**
- `AA_key` (list): 20 standard amino acids
- `AA_num_key` (numpy.ndarray): Numeric encoding [1-20]
- `AA_num_key_new` (numpy.ndarray): Normalized encoding
- `props` (numpy.ndarray): Full property matrix (62 × 20)

---

### `prop_patterning(mono_PCA, poly_PCA, mat_size=100, props=properties[1:], ridZero=False, win_size=3, returnBig=False)`
Identifies most discriminative biophysical properties between two datasets.

**Parameters:**
- `mono_PCA` (numpy.ndarray): First encoded dataset
- `poly_PCA` (numpy.ndarray): Second encoded dataset
- `mat_size` (int): Top features to retain
- `props` (numpy.ndarray): Property matrix
- `ridZero` (bool): Remove all-zero columns
- `win_size` (int): Averaging window size
- `returnBig` (bool): Return full property masks

**Returns:**
- `new_mat_mono`, `new_mat_poly` (numpy.ndarray): Discriminative features
- `max_diffs` (numpy.ndarray): Feature locations and magnitudes
- `poly_prop_masks`, `mono_prop_masks` (optional): Full property arrays

---

### `prop_pairing(ALL_mono, ALL_poly, mat_size=100, props=properties[1:], win_size=3)`
Analyzes loop-to-loop interactions and property pairings.

**Parameters:**
- `ALL_mono` (numpy.ndarray): Monoreactive sequences
- `ALL_poly` (numpy.ndarray): Polyreactive sequences
- `mat_size` (int): Features to analyze
- `props` (numpy.ndarray): Property subset
- `win_size` (int): Window for averaging

**Returns:**
- `poly_pair_score`, `mono_pair_score`: Loop interaction scores
- `poly_prop_masks`, `mono_prop_masks`: Property distributions

---

## Statistical Functions

### `parse_props(X_train, y_train, mat_size=100)`
Extracts most discriminative features via difference analysis.

**Parameters:**
- `X_train` (numpy.ndarray): Features (sequences × properties)
- `y_train` (numpy.ndarray): Labels (1=mono, 2=poly)
- `mat_size` (int): Number of top features

**Returns:**
- `max_diffs` (numpy.ndarray): Top features (value, index, location)

---

### `do_statistics(data1, data2, num_reps=1000, test='average', test_func=None, func_val=None)`
Performs permutation-based statistical testing.

**Parameters:**
- `data1`, `data2` (numpy.ndarray): Datasets to compare
- `num_reps` (int): Permutation replicates
- `test` (str): 'average'/'function'
- `test_func` (callable): Custom function for 'function' test
- `func_val` (int): Index to extract from function output

**Returns:**
- `p_values` (numpy.ndarray): P-values per position

---

## Special Analyses (MHC/TCR)

### `encode_meta(metadat)`
Encodes categorical metadata to numeric IDs.

**Parameters:**
- `metadat` (pandas.DataFrame): Categorical metadata

**Returns:**
- `meta_encoded` (pandas.DataFrame): Numeric encoding

---

### `randomize_tcr_mhc_pair(trv_mat_df, mhc_mat_df, multiOrg=False, mhc_orgs=['Human'])`
Randomizes TCR-MHC pairings while preserving properties.

**Parameters:**
- `trv_mat_df`, `mhc_mat_df` (pandas.DataFrame): Encoded sequences
- `multiOrg` (bool): Handle multiple organisms
- `mhc_orgs` (list): Organism filters

**Returns:**
- Randomized paired matrices

---

### `get_mhcSub(mhc, seq_choice, multiOrg=False, input_loc=[])`
Extracts peptide-binding regions from MHC structures.

**Parameters:**
- `mhc` (str): 'classI', 'classIIa', or 'classIIb'
- `seq_choice` (list): MHC sequences
- `multiOrg` (bool): Track organism origin
- `input_loc` (list): Organism labels

**Returns:**
- Subset sequences for specified MHC class
- Organism mapping (if multiOrg=True)

---

### `get_distClusts(dists, metadat, max_d=5)`
Clusters sequences by computed distance (AIMS distance clustering).

**Parameters:**
- `dists` (numpy.ndarray): Distance matrix
- `metadat` (pandas.DataFrame): Metadata
- `max_d` (float): Distance threshold for clustering

**Returns:**
- `distance_clusters` (pandas.DataFrame): Cluster assignments

---

## Usage Examples

```python
from aims_immune import aims_analysis as aims
import numpy as np

# Encode TCR sequences with center alignment
tcr_matrix = aims.gen_tcr_matrix(
    sequences,
    alignment='center',
    bulge_pad=8
)

# Calculate entropy
shannon, frequencies, coverage = aims.calculate_shannon(tcr_matrix)

# Find discriminative properties
mono_props, poly_props, max_diffs = aims.prop_patterning(
    monoreactive_encoded,
    polyreactive_encoded,
    mat_size=50
)

# Statistical test
p_vals = aims.do_statistics(
    dataset1,
    dataset2,
    num_reps=10000,
    test='average'
)
```

---

## Important Notes

- Sequence encoding: 0=gap/padding, 1-20=amino acids (alphabetical)
- Shannon entropy typically 0-4.32 bits for 20 amino acids
- MI units are in bits; values typically 0-2 for sequence positions
- Properties are normalized (mean-subtracted, unit-length)
- All calculations support parallel processing
- Gap character only used in MSA analysis
