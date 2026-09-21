# Fix Recommendations for aims_analysis Tests

## Overview

This document provides specific recommendations for fixing the 17 failing tests in `test_aims_analysis.py`.

## Critical Fixes (High Priority - 8 tests)

### 1. Type Mismatch Issues - DataFrame vs NumPy Array (6 tests)

**Affected Tests:**
- `TestCalcAimsDist` (3 tests): test_calc_aims_dist_basic, test_calc_aims_dist_identical_matrices, test_calc_aims_dist_different_matrices
- `TestGenCloneProps` (1 test): test_gen_clone_props_basic
- `TestGenDsetProps` (1 test): test_gen_dset_props_basic
- `TestDoStatistics` (1 test): test_do_statistics_basic

**Problem:** Functions expect specific data types (DataFrame, ndarray) but tests provide wrong types or use incompatible operations.

**Solutions:**

#### For calc_AIMSdist tests:
```python
# WRONG (current):
mat1 = np.random.rand(100, 50)
mat2 = np.random.rand(100, 50)
result = analysis.calc_AIMSdist(mat1, mat2)

# RIGHT (fixed):
mat1_df = pd.DataFrame(np.random.rand(100, 50))
mat2_df = pd.DataFrame(np.random.rand(100, 50))
result = analysis.calc_AIMSdist(mat1_df, mat2_df)
```

#### For gen_clone_props and gen_dset_props:
```python
# WRONG (current):
seqF = pd.DataFrame({'seq1': ['AAAA'], 'seq2': ['GGGG']})
result = analysis.gen_clone_props(seqF)

# RIGHT (fixed):
seqF = pd.DataFrame({'seq1': ['AAAA'], 'seq2': ['GGGG']})
pre_poly = seqF.values.astype('U10')  # Convert to numpy array
result = analysis.gen_clone_props(pre_poly)
```

#### For do_statistics:
```python
# WRONG (current):
labels1 = ['A'] * 50
labels2 = ['B'] * 50

# RIGHT (fixed):
labels1 = np.arange(50)  # Use numeric labels
labels2 = np.arange(50, 100)
```

---

### 2. Missing Required Parameters (2 tests)

**Affected Tests:**
- `TestFullAAFreq` (2 tests): test_full_aa_freq_basic, test_full_aa_freq_known_sequences

**Problem:** Function signature requires `my_AA_key` parameter which tests don't provide.

**Function Signature:**
```python
def full_AA_freq(seqF, my_AA_key, ...):
```

**Fix:**
```python
# WRONG (current):
result = analysis.full_AA_freq(seqF)

# RIGHT (fixed):
result = analysis.full_AA_freq(seqF, my_AA_key=analysis.AA_key)
```

---

### 3. Incompatible Test Data (2 tests)

**Affected Tests:**
- `TestGenPeptideMatrix` (2 tests): test_gen_peptide_matrix_basic, test_gen_peptide_matrix_shape
- `TestIntegration` (1 test): test_full_analysis_pipeline

**Problem:** Test data (CDR regions) doesn't match peptide matrix requirements.

**Fix - Option A: Use Real Peptide Data**
```python
# Use real test data
def test_gen_peptide_matrix_basic(self, peptide_data_paths):
    pep_data = loader.seq_loader(str(peptide_data_paths['kidney']), label='peptide')
    pre_pep = [pep_data.iloc[i].tolist() for i in range(min(3, len(pep_data)))]
    result = analysis.gen_peptide_matrix(pre_pep)
    assert isinstance(result, np.ndarray)
```

**Fix - Option B: Create Proper Test Data**
```python
# Create properly formatted peptide data
peptide_sequences = [
    ['MKVLWAALLVGTGSASHSR', 'MEQLTLLLPLATSLALTASM', 'QSISSYDASQHRSTWPPNG'],
    # ... (multiple columns of peptide sequences)
]
result = analysis.gen_peptide_matrix(peptide_sequences)
```

---

## Important Fixes (Medium Priority - 5 tests)

### 4. Function API Misunderstandings (4 tests)

**Affected Tests:**
- `TestGenTcrMatrix` (1 test): test_gen_tcr_matrix_with_size_specification
- `TestGenTcrMatrix` (1 test): test_gen_tcr_matrix_with_binary_mode
- `TestGetSequenceDimension` (2 tests): test_get_sequence_dimension_basic, test_get_sequence_dimension_single_column
- `TestCreateMsaPairs` (1 test): test_create_msa_pairs_basic

#### Fix for gen_tcr_matrix_with_size_specification:
```python
# WRONG (current):
giveSize = 10  # Integer
result = analysis.gen_tcr_matrix(pre_poly, giveSize=giveSize)

# RIGHT (fixed):
giveSize = [10, 10, 10]  # List with size for each loop
result = analysis.gen_tcr_matrix(pre_poly, giveSize=giveSize)
```

#### Fix for gen_tcr_matrix_with_binary_mode:
```python
# WRONG (current):
result = analysis.gen_tcr_matrix(pre_poly, pre_mono=pre_mono, binary=True)
assert isinstance(result, np.ndarray)

# RIGHT (fixed):
result = analysis.gen_tcr_matrix(pre_poly, pre_mono=pre_mono, binary=True)
# Result is tuple when binary=True
assert isinstance(result, tuple)
poly_mat, mono_mat = result
```

#### Fix for get_sequence_dimension tests:
```python
# WRONG (current):
result = analysis.get_sequence_dimension(sample_dataframe)
assert len(result) == 6  # Assumes 6 columns

# RIGHT (fixed):
result = analysis.get_sequence_dimension(sample_dataframe)
assert len(result) == sample_dataframe.shape[0]  # Number of rows
```

#### Fix for create_msa_pairs:
```python
# WRONG (current):
result = analysis.create_msa_pairs(seqF)

# RIGHT (fixed):
# Check function signature and provide all required parameters
result = analysis.create_msa_pairs(
    seqF,
    file2=seqF,  # Second file
    pair_list=[0, 1, 2],  # Pairing information
    names=['pair1', 'pair2', 'pair3']  # Names
)
```

---

### 5. Return Type Mismatches (1 test)

**Affected Test:**
- `TestEncodeMetadata` (1 test): test_encode_meta_basic

**Problem:** Function returns pandas Series, test expects ndarray or list.

**Fix:**
```python
# WRONG (current):
result = analysis.encode_meta(metadata)
assert isinstance(result, (np.ndarray, list))

# RIGHT (fixed):
result = analysis.encode_meta(metadata)
assert isinstance(result, (np.ndarray, list, pd.Series))
```

---

## Implementation Priority

### Phase 1 - Critical Fixes (1-2 hours)
1. Fix type mismatches for calc_AIMSdist, gen_clone_props, gen_dset_props
2. Add missing parameters to full_aa_freq tests
3. Fix do_statistics to use numeric data

### Phase 2 - Important Fixes (1-2 hours)
4. Update size parameter handling for gen_tcr_matrix
5. Fix get_sequence_dimension test assumptions
6. Update return type assertions

### Phase 3 - Cleanup (30 minutes - 1 hour)
7. Replace incompatible test data with real data
8. Update create_msa_pairs with proper parameters
9. Review and finalize all assertions

---

## Testing the Fixes

After implementing fixes, verify with:

```bash
# Run specific fixed test
pytest test_aims_analysis.py::TestCalcAimsDist::test_calc_aims_dist_basic -v

# Run all tests in a class
pytest test_aims_analysis.py::TestCalcAimsDist -v

# Run full module with summary
pytest test_aims_analysis.py -v --tb=short
```

---

## Additional Recommendations

### For Test Robustness:
1. **Use real test data:** Replace synthetic test data with actual data from `app_data/test_data/`
2. **Add parametrized tests:** Use pytest.mark.parametrize for multiple input variations
3. **Document expectations:** Add comments explaining expected inputs and outputs

### For Function Robustness:
1. **Type hints:** Add type annotations to functions
2. **Input validation:** Check input types early and provide clear error messages
3. **Flexible input:** Consider supporting multiple input types (DataFrame, ndarray)

### Example Enhanced Test:
```python
def test_calc_aims_dist_with_real_data(self, abs_data_paths):
    """Test distance calculation with real antibody data"""
    # Load real data
    poly_df = loader.seq_loader(str(abs_data_paths['poly']), label='poly')
    mono_df = loader.seq_loader(str(abs_data_paths['mono']), label='mono')
    
    # Generate matrices (test-specific helper)
    poly_matrix = analysis.gen_peptide_matrix(...)
    mono_matrix = analysis.gen_peptide_matrix(...)
    
    # Test calculation
    distance = analysis.calc_AIMSdist(poly_matrix, mono_matrix)
    
    # Verify result
    assert isinstance(distance, (float, np.floating))
    assert distance >= 0
```

---

## Quick Reference - Severity Levels

| Severity | Count | Impact | Time to Fix |
|----------|-------|--------|------------|
| **High** | 8 | Function not testable as-is | 2-3 hours |
| **Medium** | 5 | Test logic incorrect | 1-2 hours |
| **Low** | 4 | Minor assertion issue | 15-30 min |

---

## Next Steps

1. **Review** this document with the team
2. **Prioritize** fixes based on module dependencies
3. **Implement** Phase 1 fixes first
4. **Verify** each fix with local test runs
5. **Document** any function API changes discovered
6. **Update** function docstrings if needed

---

**Last Updated:** 2026-09-08  
**Test Module:** aims_immune/test/test_aims_analysis.py  
**Failing Tests:** 17/31 (55%)
