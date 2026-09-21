# AIMS Immune - aims_utils Module API

## Module Overview
The `aims_utils` module provides utility functions specifically designed for Jupyter notebook and command-line analysis workflows. It wraps core analysis functions with enhanced usability and visualization capabilities.

---

## Data Loading & Preparation

### `loadDat(fileName, datName, datDir, drop_duplicates, subStart=[], subEnd=[], subset=False)`
Loads and preprocesses multiple sequence files.

**Parameters:**
- `fileName` (list): List of input filenames
- `datName` (list): Dataset labels/names
- `datDir` (str): Directory containing files
- `drop_duplicates` (bool): Remove duplicate sequences
- `subStart` (list): Start positions for subsetting
- `subEnd` (list): End positions for subsetting
- `subset` (bool): Enable subsetting

**Returns:**
- `seqF` (pandas.DataFrame): Concatenated sequences
- `xtick_loc` (list): X-tick positions for visualization
- `AA_num_key` (numpy.ndarray): Amino acid numeric encoding
- `mat_size` (list): Matrix dimensions per loop

---

### `get_keys(custom_key='')`
Defines amino acid ordering for visualization.

**Parameters:**
- `custom_key` (str): Custom 20-letter AA sequence; empty = default

**Returns:**
- `my_AA_key` (list): Ordered amino acids
- `my_AA_key_dash` (list): AA key + gap character

**Default Key Order:**
- A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V

---

## Feature Matrix Generation

### `grab_big(seqF, dsetF, seq_MIf, datName, normalize, renormalize, molecule, my_AA_key, my_AA_key_dash, mat_size, align, pad, nCores)`
Creates comprehensive feature matrix with property filtering and normalization.

**Parameters:**
- `seqF` (pandas.DataFrame): Sequences
- `dsetF` (numpy.ndarray): Dataset array form
- `seq_MIf` (pandas.DataFrame): Encoded sequence matrix
- `datName` (list): Dataset names
- `normalize` (str): Property normalization ('msuv'/'zscore'/'0to1')
- `renormalize` (bool): Renormalize by entropy
- `molecule` (str): Molecule type ('ig'/'msa'/'peptide')
- `my_AA_key` (list): Amino acid alphabet
- `my_AA_key_dash` (list): AA + gap character
- `mat_size` (list): Matrix dimensions
- `align` (str): Alignment method
- `pad` (int): Bulge padding
- `nCores` (int): Number of processing threads

**Returns:**
- `full_big` (pandas.DataFrame): Complete feature matrix
- `parsed_mat` (pandas.DataFrame): Correlation-filtered features
- `NonNorm_big` (pandas.DataFrame): Unrenormalized features
- `IDed_full_big` (pandas.DataFrame): Features with dataset IDs
- `seq_bigReshape` (numpy.ndarray): Reshaped features (seqs × props × positions)
- `token_df` (pandas.DataFrame): Dataset ID mapping

**Feature Generation Process:**
1. Generate biophysical property masks
2. Remove highly correlated features (>0.75 correlation)
3. Remove all-zero feature columns
4. Optional: Renormalize by Shannon entropy
5. Flatten to feature vectors

---

## Dimensionality Reduction & Clustering

### `do_cluster(chosen_dset, seq_MIf, reduce, clust, NClust=100, min_samples=10, eps=0.5, nComp=3, n_neighbors=25, state=617)`
Performs dimensionality reduction and clustering.

**Parameters:**
- `chosen_dset` (pandas.DataFrame or numpy.ndarray): Input features
- `seq_MIf` (pandas.DataFrame): Original sequence matrix
- `reduce` (str): Reduction method ('pca'/'umap')
- `clust` (str): Clustering algorithm ('kmean'/'optics'/'dbscan'/'hdbscan')
- `NClust` (int): Number of clusters (KMeans only)
- `min_samples` (int): Minimum samples per cluster (OPTICS/DBSCAN)
- `eps` (float): Epsilon distance threshold (DBSCAN)
- `nComp` (int): Number of components (UMAP)
- `n_neighbors` (int): Neighborhood size (UMAP, default 25)
- `state` (int): Random seed for reproducibility

**Returns:**
- `clust_input` (numpy.ndarray): Reduced-dimension coordinates
- `cluster_dset` (pandas.DataFrame): Cluster assignments

**Algorithms:**
- `'pca'`: PCA to 3 components
- `'umap'`: UMAP to specified components
- `'kmean'`: K-Means clustering
- `'optics'`: OPTICS density clustering
- `'dbscan'`: DBSCAN density clustering
- `'hdbscan'`: Hierarchical DBSCAN (scikit-learn 1.3+)

---

## Metadata Management

### `get_metadata(token_df, datName, cluster_dset=[], compile_index=[], metaPath='', meta_form='category', gotMeta=False)`
Incorporates metadata for visualization and analysis.

**Parameters:**
- `token_df` (pandas.DataFrame): Dataset ID tokens
- `datName` (list): Dataset names
- `cluster_dset` (pandas.DataFrame): Cluster assignments (optional)
- `compile_index` (list): Sequence indices
- `metaPath` (str): Path to external metadata CSV
- `meta_form` (str): 'category' (default) or 'quant'
- `gotMeta` (bool): External metadata loaded

**Returns:**
- `clust_map` (pandas.DataFrame): Cluster IDs per sequence
- `clust_leg` (list): Unique cluster IDs
- `clust_name` (str): Cluster column name
- `meta_map` (pandas.DataFrame): Metadata per sequence
- `meta_leg` (list): Unique metadata values
- `meta_name` (str): Metadata column name
- `metadat` (pandas.DataFrame): Original metadata

---

## Visualization Functions

### `plot_clusters(clust_input, plot_metas, clust, show_labels=False, proj_show='both', clust_show='both', need_meta=True, labels=[])`
Creates interactive cluster visualizations (2D and/or 3D).

**Parameters:**
- `clust_input` (numpy.ndarray): Reduced-dimension data
- `plot_metas` (tuple): Metadata from get_metadata()
- `clust` (str): Clustering algorithm used
- `show_labels` (bool): Annotate points with labels
- `proj_show` (str): '2d'/'3d'/'both'
- `clust_show` (str): 'cluster'/'metadata'/'both'
- `need_meta` (bool): Include metadata coloring
- `labels` (list): Sequence labels for annotation

**Returns:**
- `fig3d` (matplotlib.figure): Figure object
- `ax` (list): Subplot axes
- `keep_map` (numpy.ndarray): Color mapping (if need_meta=True)

---

### `plot_quantClust(cluster_dset, keep_map, plot_metas, norm=True, cmap3=colormap)`
Creates stacked bar plots of cluster composition.

**Parameters:**
- `cluster_dset` (pandas.DataFrame): Cluster assignments
- `keep_map` (numpy.ndarray): Color mapping
- `plot_metas` (tuple): Metadata from get_metadata()
- `norm` (bool): Normalize to fractions (True) or counts (False)
- `cmap3` (matplotlib.colormap): Color palette

**Returns:**
- `fig` (matplotlib.figure): Figure object
- `final_breakdown` (pandas.DataFrame): Cluster × metadata cross-tabulation
- `meta_legF` (numpy.ndarray): Metadata legend

---

### `plotSubsets(seq_bigReshape, seq_MIf, plot_metas, mat_size, subset_sel='cluster', plot_props=False, showProp=1, show_lines=False)`
Visualizes sequences organized by cluster or metadata.

**Parameters:**
- `seq_bigReshape` (numpy.ndarray): Reshaped feature matrix
- `seq_MIf` (pandas.DataFrame): Encoded sequence matrix
- `plot_metas` (tuple): Metadata from get_metadata()
- `mat_size` (list): Dimensions per region
- `subset_sel` (str): 'cluster' or 'metadata'
- `plot_props` (bool): Color by property (True) or encoding (False)
- `showProp` (int): Property index (if plot_props=True)
- `show_lines` (bool): Draw boundaries between groups

**Returns:**
- `fig` (matplotlib.figure): Figure object
- `chosen_map` (pandas.DataFrame): Subset assignments
- `chosen_name` (str): Subset column name

---

## Statistical Analysis

### `run_AIMSdist(full_big, seq_MIf, seqF, plot_metas, chosen_map, chosen_name, parallel_dist=True, get_distClusts=True, maxD=5, nThreads=-1, reorder=False, noFig=False)`
Calculates sequence distances and identifies distance-based clusters.

**Parameters:**
- `full_big` (pandas.DataFrame): Feature matrix
- `seq_MIf` (pandas.DataFrame): Encoded sequences
- `seqF` (pandas.DataFrame): Original sequences
- `plot_metas` (tuple): Metadata
- `chosen_map` (pandas.DataFrame): Cluster/metadata assignment
- `chosen_name` (str): Cluster column name
- `parallel_dist` (bool): Use parallel processing
- `get_distClusts` (bool): Cluster by distance
- `maxD` (float): Distance threshold
- `nThreads` (int): Parallel threads
- `reorder` (bool): Sort by cluster
- `noFig` (bool): Skip figure generation

**Returns:**
- `dists` (numpy.ndarray): Pairwise distance matrix
- `distance_clusters` (pandas.DataFrame): Distance-based clusters (optional)
- `fig` (matplotlib.figure): Distance heatmap (optional)

---

### `newParallel(full_big, num_threads=-1)`
Numba-accelerated Euclidean distance calculation.

**Parameters:**
- `full_big` (numpy.ndarray): Feature matrix
- `num_threads` (int): Parallel threads

**Returns:**
- `dist_calc` (numpy.ndarray): Pairwise distance matrix

---

## Biophysical Property Analysis

### `netStats(ref_sub, sub_big, posLen, sub_matF, prop_names=['Charge','Hydrophobicity','Bulkiness','Flexibility'], reps=1000)`
Calculates statistical significance for averaged biophysical properties.

**Parameters:**
- `ref_sub` (pandas.DataFrame): Subset reference with group assignments
- `sub_big` (pandas.DataFrame): Subset features
- `posLen` (int): Number of positions
- `sub_matF` (pandas.DataFrame): Subset encoding matrix
- `prop_names` (list): Properties to test
- `reps` (int): Permutation replicates

**Returns:**
- `p_list` (list): P-values per property

---

### `do_netAvg(sub_big, ref_sub, sub_matF, sub_sels, posLen, label, colors=['Crimson','darkorchid'], bootstrap=False, boots=1000, stats=False)`
Plots network-averaged biophysical properties.

**Parameters:**
- `sub_big` (pandas.DataFrame): Subset features
- `ref_sub` (pandas.DataFrame): Subset assignments
- `sub_matF` (pandas.DataFrame): Encoding matrix
- `sub_sels` (list): Subset indices to plot
- `posLen` (int): Number of positions
- `label` (list): Subset labels
- `colors` (list): Color per subset
- `bootstrap` (bool): Use bootstrap resampling
- `boots` (int): Bootstrap replicates
- `stats` (bool): Include statistics

**Returns:**
- `fig` (matplotlib.figure): Bar plot
- `p_list` (list): P-values (if stats=True)

---

### `do_posAvg(sub_sels, ref_sub, sub_matF, sub_big, posLen, mat_size, label, prop_sel=[1,2], colors=['Crimson','darkorchid'], bootstrap=False, boots=1000)`
Plots position-sensitive property averages.

**Parameters:**
- Similar to do_netAvg
- `prop_sel` (list): Property indices [1=Charge, 2=Hydrophobicity]
- Creates separate plots per property

**Returns:**
- `fig` (matplotlib.figure): Line plots with error bands

---

### `pos_Shannon(seq_MIf, sub_sels, sub_matF, ref_sub, mat_size, label, colors=['Crimson','darkorchid'], bootstrap=False, boots=1000)`
Plots Shannon entropy and coverage by subset.

**Parameters:**
- Similar structure to position analysis functions
- Top panel: Coverage (non-zero positions)
- Bottom panel: Entropy values

**Returns:**
- `fig` (matplotlib.figure): Stacked plot
- `entropy` (list): Entropy per subset
- `frequencies` (list): AA frequencies per subset

---

### `pos_MI(seq_MIf, sub_sels, sub_matF, ref_sub, label, mat_size)`
Visualizes mutual information matrices by subset.

**Parameters:**
- `seq_MIf` (pandas.DataFrame): Encoded matrix
- `sub_sels` (list): Subset indices
- `sub_matF` (pandas.DataFrame): Subset matrix
- `ref_sub` (pandas.DataFrame): Assignments
- `label` (list): Subset labels
- `mat_size` (list): Dimensions

**Returns:**
- `fig` (matplotlib.figure): MI heatmaps
- `MI` (list): MI matrices per subset
- `ent_cond` (list): Conditional entropy

---

## Sequence Output

### `back_decode(seq_MIf, outputDir, my_AA_key, AA_num_key, outName='AIMS_encodedSeqs.csv')`
Converts numeric encoding back to amino acid sequences.

**Parameters:**
- `seq_MIf` (pandas.DataFrame): Encoded matrix
- `outputDir` (str): Output directory
- `my_AA_key` (list): AA alphabet
- `AA_num_key` (numpy.ndarray): Numeric encoding
- `outName` (str): Output filename

**Returns:**
- Writes CSV file with decoded sequences

---

### `viz_subs(dset, seq_MIf, meta_legF, sub_sels, chosen_map, chosen_name, seqlogo=False, save_subSeqs=False, saveAll=False)`
Visualizes selected subsets of sequences.

**Parameters:**
- `dset` (pandas.DataFrame): Original sequences
- `seq_MIf` (pandas.DataFrame): Encoded matrix
- `meta_legF` (list): Metadata labels
- `sub_sels` (list): Subset indices to visualize
- `chosen_map` (pandas.DataFrame): Assignments
- `chosen_name` (str): Assignment column
- `seqlogo` (bool): Generate seqlogo files
- `save_subSeqs` (bool): Save subset sequences
- `saveAll` (bool): Visualize all subsets

**Returns:**
- `fig` (matplotlib.figure): Subset visualization
- `sub_matF` (pandas.DataFrame): Subset encoding
- `label` (list): Subset labels
- `sub_seqF` (pandas.DataFrame): Subset sequences (if save_subSeqs=True)

---

## Usage Examples

```python
from aims_immune import aims_utils as utils

# Load data and create features
seqF, xticks, AA_key, mat_size = utils.loadDat(
    ['file1.fasta', 'file2.fasta'],
    ['Dataset1', 'Dataset2'],
    './data',
    drop_duplicates=True
)

# Generate comprehensive feature matrix
full_big, parsed, _, _, reshaped, token_df = utils.grab_big(
    seqF, seqF.values, seq_MIf,
    ['Dataset1', 'Dataset2'],
    normalize='msuv',
    renormalize=True,
    molecule='ig',
    my_AA_key=get_keys()[0],
    my_AA_key_dash=get_keys()[1],
    mat_size=mat_size,
    align='center',
    pad=8,
    nCores=-1
)

# Cluster data
clust_input, cluster_dset = utils.do_cluster(
    parsed, seq_MIf, 'pca', 'kmean',
    NClust=5
)

# Visualize
metadata = utils.get_metadata(token_df, ['Dataset1', 'Dataset2'], cluster_dset)
fig, ax = utils.plot_clusters(clust_input, metadata, 'kmean')
```

---

## Important Notes

- Most functions expect pandas DataFrames for easy column manipulation
- Parallel processing default: -1 (all available threads)
- Color palettes use matplotlib colormaps
- Bootstrap confidence intervals calculated via resampling
- Distance matrix upper triangle used (symmetric)
- Metadata must align with sequence ordering
