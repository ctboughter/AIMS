# AIMS Immune - aims_loader Module API

## Module Overview
The `aims_loader` module provides functionality for loading and processing immunological sequence data in various formats (FASTA, CSV). It handles sequence subsetting and duplicate removal.

---

## Functions

### `get_msa_sub(seqF, loc_start, loc_end)`
Extracts subsequences from multiple sequence alignment data at specified locations.

**Parameters:**
- `seqF` (pandas.DataFrame): Sequence data with sequences as columns
- `loc_start` (list): Starting positions for subsequence extraction
- `loc_end` (list): Ending positions for subsequence extraction

**Returns:**
- `seqNEW` (pandas.DataFrame): DataFrame with extracted subsequences

**Raises:**
- Prints error if loc_start and loc_end have different lengths

---

### `seq_loader(seqPath, label, drop_dups=False, return_index=False, dataLoc=[], datType='', subset=False, subset_starts=[], subset_ends=[])`
Smart sequence loader that auto-detects file format (FASTA or CSV) and loads sequences with optional subsetting.

**Parameters:**
- `seqPath` (str): Path to sequence file
- `label` (str): Label/identifier for sequences
- `drop_dups` (bool): If True, remove duplicate sequences
- `return_index` (bool): If True, return original sequence indices
- `dataLoc` (list): Column indices to use from CSV (empty = all)
- `datType` (str): Format type ('fasta' or 'csv'); auto-detected if empty
- `subset` (bool): If True, extract subsequences
- `subset_starts` (list): Starting positions for subsetting
- `subset_ends` (list): Ending positions for subsetting

**Returns:**
- `fin_out` (pandas.DataFrame): Processed sequences
- `fin_id` (optional): Original sequence indices if return_index=True

**Supported Formats:**
- FASTA files (detected by '>' in first line)
- CSV files with header row
- Automatically adds '_#' suffix to sequence identifiers

---

### `convert_3Let(inp)`
Converts three-letter amino acid codes to single-letter codes.

**Parameters:**
- `inp` (list/array): Three-letter amino acid codes

**Returns:**
- `sin_final` (numpy.ndarray): Single-letter amino acid codes

**Amino Acid Mapping:**
- ALA→A, GLY→G, ARG→R, LYS→K, ASP→D, GLU→E, ASN→N, GLN→Q
- MET→M, CYS→C, PHE→F, TYR→Y, THR→T, TRP→W, PRO→P, SER→S
- LEU→L, VAL→V, HIS→H, ILE→I

---

## Data Format Requirements

### FASTA Format
```
>sequence_id description
MKVLWAALLVTFLAGCAKAKAVS...
>sequence_id2 description
MKVLWAALLVTFLAGCAKAKAVS...
```

### CSV Format
- First row contains headers
- One sequence per row
- Can contain multiple columns (CDR loops, etc.)
- Empty cells are handled gracefully

---

## Usage Examples

```python
from aims_immune import aims_loader as aimsLoad

# Load FASTA file without duplicates
seq = aimsLoad.seq_loader('sequences.fasta', 'my_dataset', drop_dups=True)

# Load CSV with subset extraction (e.g., CDR loops from positions 24-96)
seq = aimsLoad.seq_loader('data.csv', 'CDR_region', 
                          subset=True, 
                          subset_starts=[24, 50], 
                          subset_ends=[50, 96])

# Load with index preservation
seq, seq_index = aimsLoad.seq_loader('file.fasta', 'label', return_index=True)

# Convert amino acids
single_letter = aimsLoad.convert_3Let(['ALA', 'GLY', 'ARG'])
```

---

## Notes
- Automatic format detection eliminates need to manually specify datType
- The module removes entries with 'X' amino acids (non-standard)
- Incomplete entries (missing data) are filtered
- Sequence identifiers are automatically labeled with dataset identifier and index
