# AIMS Immune - aims_classification Module API

## Module Overview
The `aims_classification` module provides machine learning classification functionality for immune repertoire analysis, including feature selection, dimensionality reduction, and classifier training using biophysical properties of protein sequences.

---

## Core Functions

### `apply_matrix(mono_PCA, max_diffs, mat_size=100, props=properties_global[1:], ridZero=False, win_size=3)`
Applies biophysical property matrices to sequences by selecting discriminative features.

**Parameters:**
- `mono_PCA` (numpy.ndarray): Encoded sequence matrix
- `max_diffs` (numpy.ndarray): Feature locations and properties (from parse_props)
- `mat_size` (int): Number of top discriminative features to retain
- `props` (numpy.ndarray): Biophysical property matrix
- `ridZero` (bool): Remove all-zero columns before processing
- `win_size` (int): Sliding window size for averaging

**Returns:**
- `new_mat_mono` (numpy.ndarray): Property-encoded matrix

---

### `getBig(mono_PCA, properties, norm='msuv', num_threads=-1)`
Numba-accelerated function for computing biophysical property masks across all sequences.

**Parameters:**
- `mono_PCA` (numpy.ndarray): Encoded sequence matrix
- `properties` (numpy.ndarray): Biophysical property matrix (shape: props × AAs)
- `norm` (str): Normalization method ('msuv'/'zscore'/'0to1')
- `num_threads` (int): Number of threads (-1 = all available)

**Returns:**
- `mono_prop_masks` (numpy.ndarray): Property masks (shape: props × sequences × positions)

**Normalization Options:**
- `'msuv'`: Mean-subtracted unit vector (default)
- `'zscore'`: Z-score normalization
- `'0to1'`: Min-max scaling to [0, 1]

---

### `get_bigass_matrix(ALL_mono, AA_key=AA_key, AA_key_dash=AA_key_dash, OneChain=False, giveSize=[], onlyCen=False, bulge_pad=8, prop_parse=False, manuscript_arrange=False, special='', alignment='center', norm='msuv', nCores=-1)`
Creates comprehensive feature matrix combining encoded sequences with biophysical properties.

**Parameters:**
- `ALL_mono` (numpy.ndarray): Input sequences
- `AA_key` (list): Standard amino acid key [20 AAs]
- `AA_key_dash` (list): AA key with gap character [21 symbols]
- `OneChain` (bool): Single chain analysis (vs TCR dual chains)
- `giveSize` (list): Matrix dimensions; auto-calculated if empty
- `onlyCen` (bool): Keep only central regions (remove 4 positions each side)
- `bulge_pad` (int): Padding for bulge alignment (2, 4, 6, or 8)
- `prop_parse` (bool): Use only older properties (16) vs all (62)
- `manuscript_arrange` (bool): Rearrange CDR loops for publication figures
- `special` (str): Special encoding ('peptide', 'MSA', or '')
- `alignment` (str): CDR alignment ('center', 'left', 'right', 'bulge')
- `norm` (str): Property normalization method
- `nCores` (int): Number of threads for parallel processing

**Returns:**
- `mono_pca_stack` (numpy.ndarray): Stacked features (sequences × features)

---

### `do_classy_mda(ALL_mono, ALL_poly, matsize=100, OneChain=False, special='', xVal='kfold', ridCorr=False, feat_sel='none', classif='mda')`
Performs multi-dimensional analysis (MDA) classification with cross-validation.

**Parameters:**
- `ALL_mono` (numpy.ndarray): Monoreactive sequences
- `ALL_poly` (numpy.ndarray): Polyreactive sequences
- `matsize` (int): Matrix size for feature selection
- `OneChain` (bool): Single vs dual chain analysis
- `special` (str): Special encoding mode ('peptide', 'MSA', '')
- `xVal` (str): Cross-validation strategy ('kfold'/'loo'/'strat_kfold')
- `ridCorr` (bool): Remove highly correlated features (>0.75)
- `feat_sel` (str): Feature selection ('none'/'PCA'/'kPCA'/'kbest'/'max_diff')
- `classif` (str): Classifier ('mda'/'svm'/'logReg'/'forest')

**Returns:**
- `acc_fin` (numpy.ndarray): Accuracy scores per fold

**Cross-Validation Methods:**
- `'kfold'`: 10-fold cross-validation
- `'loo'`: Leave-one-out
- `'strat_kfold'`: Stratified 10-fold

**Classifiers:**
- `'mda'`: Linear Discriminant Analysis
- `'svm'`: Support Vector Machine
- `'logReg'`: Logistic Regression
- `'forest'`: Random Forest (500 trees)

---

### `do_linear_split(test_mono, test_poly, ridCorr=True, giveSize=[], matSize=75, prop_parse=False, manuscript_arrange=False, pca_split=False, special='', align='center', got_big=False)`
Linear discriminant analysis without cross-validation for comparing two datasets.

**Parameters:**
- `test_mono` (numpy.ndarray or pandas.DataFrame): First dataset
- `test_poly` (numpy.ndarray or pandas.DataFrame): Second dataset
- `ridCorr` (bool): Remove highly correlated features
- `giveSize` (list): Matrix dimensions
- `matSize` (int): Number of top features
- `prop_parse` (bool): Use subset of properties
- `manuscript_arrange` (bool): Arrange for publication
- `pca_split` (bool): Use PCA instead of difference-based feature selection
- `special` (str): Special encoding ('peptide', 'MSA')
- `align` (str): Alignment method
- `got_big` (bool): Input is already bigass matrix

**Returns:**
- Varies by ridCorr and pca_split parameters; includes:
  - `bigF` (pandas.DataFrame): Full feature matrix
  - `weights` (numpy.ndarray): LDA coefficients
  - `acc_all` (float): Classification accuracy
  - `mda_all` (numpy.ndarray): LDA transformed data
  - `final` (pandas.DataFrame): Processed features
  - `top_names` (list): Top discriminative features

---

### `classy_apply(test_mat, y_test, train_mat, y_train, matsize=100, OneChain=False, ridCorr=False, feat_sel='none', classif='mda')`
Apply pre-trained classifier to new test data.

**Parameters:**
- `test_mat` (numpy.ndarray): Test sequences
- `y_test` (numpy.ndarray): Test labels
- `train_mat` (numpy.ndarray): Training sequences
- `y_train` (numpy.ndarray): Training labels
- `matsize` (int): Feature matrix size
- `OneChain` (bool): Single/dual chain
- `ridCorr` (bool): Remove correlated features
- `feat_sel` (str): Feature selection method
- `classif` (str): Classifier type

**Returns:**
- `acc_all` (float): Classification accuracy

---

### `apply_pretrained_LDA(bigass_mono, top_names, weights, prop_parse=False)`
Apply previously trained LDA model to new feature matrices.

**Parameters:**
- `bigass_mono` (numpy.ndarray): Features for new sequences
- `top_names` (list): Feature names used in training
- `weights` (numpy.ndarray): LDA coefficients
- `prop_parse` (bool): Property subset flag

**Returns:**
- `final_apply` (numpy.ndarray): LDA scores for new data

---

## Biophysical Properties

The module includes 62 total biophysical properties:
- **Properties 0-1**: Custom amino acid encodings
- **Properties 2-15**: Traditional properties (Charge, Phobic1/2, Bulk, Flexibility, etc.)
- **Properties 16-61**: Hot-spot prediction properties (46 orthogonal dimensions)

---

## Usage Examples

```python
from aims_immune import aims_classification as classy
import numpy as np

# Create and train classifier
mono_seqs = np.random.randint(1, 21, (100, 50))
poly_seqs = np.random.randint(1, 21, (150, 50))

accuracies = classy.do_classy_mda(
    mono_seqs, 
    poly_seqs, 
    matsize=50,
    feat_sel='max_diff',
    classif='svm'
)

# Perform linear discriminant analysis
full_big, weights, acc, mda, parsed, top_names = classy.do_linear_split(
    mono_seqs, 
    poly_seqs,
    ridCorr=True,
    pca_split=False
)

# Apply pre-trained model
new_scores = classy.apply_pretrained_LDA(
    new_features, 
    top_names, 
    weights
)
```

---

## Important Notes

- Sequence labels (1-20) represent amino acids in alphabetical order
- Label 0 represents gaps or padding
- Feature selection methods require sufficient data (especially 'max_diff')
- Parallel processing uses all available threads by default (num_threads=-1)
- Correlation threshold for feature removal is fixed at 0.75
- MSA analysis requires AA_key_dash with gap character
