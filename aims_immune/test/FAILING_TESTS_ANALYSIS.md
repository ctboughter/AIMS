# aims_analysis Test Failures - Detailed Analysis

## Summary

**Total Tests:** 31  
**Passing:** 14 (45%)  
**Failing:** 17 (55%)  
**Execution Time:** 1.76s

---

## Failing Tests Breakdown

### 1. TestGetSequenceDimension (2 failures)

#### test_get_sequence_dimension_basic ❌
**Location:** aims_immune/test/test_aims_analysis.py:16  
**Status:** FAILED

**Error:**
```
assert len(result) == 6  # 6 columns in sample data
assert 3 == 6
where 3 = len([np.int64(16), np.int64(18), np.int64(19)])
```

**Root Cause:** The sample_dataframe has 3 sequence types (l1, l2, l3, h1, h2, h3) but when transposed for this function, it returns 3 results, not 6. The test assumption about input structure is incorrect.

**Fix Needed:** Adjust test to match actual behavior - the function returns dimensions for each row, not each column.

---

#### test_get_sequence_dimension_single_column ❌
**Location:** aims_immune/test/test_aims_analysis.py:30  
**Status:** FAILED

**Error:**
```
assert len(result) == 1
assert 3 == 1
where 3 = len([np.int64(9), np.int64(9), np.int64(9)])
```

**Root Cause:** Same as above - misunderstanding of function behavior. The function returns one dimension per row in the input DataFrame.

**Fix Needed:** Either pass a properly structured DataFrame or adjust assertion.

---

### 2. TestGenTcrMatrix (2 failures)

#### test_gen_tcr_matrix_with_size_specification ❌
**Location:** aims_immune/test/test_aims_analysis.py:76  
**Status:** FAILED

**Error:**
```
IndexError: list index out of range
File "aims_immune/aims_analysis.py", line 134, in gen_tcr_matrix
    LoopMax = max_len[k]
              ^^^^^^^^^^
```

**Root Cause:** The function expects `giveSize` to be a list with multiple elements, but the test passes a single integer. When the code tries to index max_len[k], it fails.

**Fix Needed:** Either pass a list of sizes or ensure the function handles integer input properly.

---

#### test_gen_tcr_matrix_with_binary_mode ❌
**Location:** aims_immune/test/test_aims_analysis.py:93  
**Status:** FAILED

**Error:**
```
AssertionError: assert False
where False = isinstance((array([...]), array([...])), <class 'numpy.ndarray'>)
```

**Root Cause:** When `binary=True`, the function returns a tuple of arrays (poly_PCA, mono_PCA), not a single ndarray. The test expects a single ndarray.

**Fix Needed:** Adjust assertion to check for tuple or unpack the return values.

---

### 3. TestGenPeptideMatrix (2 failures)

#### test_gen_peptide_matrix_basic ❌
**Location:** aims_immune/test/test_aims_analysis.py:130  
**Status:** FAILED

**Error:**
```
IndexError: index 14 is out of bounds for axis 0 with size 14
File "aims_immune/aims_analysis.py", line 836, in gen_peptide_matrix
    pep_PCA[i][count]=key[j]
    ^^^^^^^^^^^^^^^^^
```

**Root Cause:** The test data (l1, l2, l3) doesn't have enough columns for peptide matrix generation. The function expects peptide sequences in standard format.

**Fix Needed:** Use proper peptide test data or adjust test expectations.

---

#### test_gen_peptide_matrix_shape ❌
**Location:** aims_immune/test/test_aims_analysis.py:137  
**Status:** FAILED

**Error:**
```
Same as test_gen_peptide_matrix_basic
IndexError: index 14 is out of bounds for axis 0 with size 14
```

**Root Cause:** Same as above - incompatible test data structure.

**Fix Needed:** Use proper test data for peptide matrices.

---

### 4. TestCalcAimsDist (3 failures)

#### test_calc_aims_dist_basic ❌
**Location:** aims_immune/test/test_aims_analysis.py:168  
**Status:** FAILED

**Error:**
```
AttributeError: 'numpy.ndarray' object has no attribute 'equals'
File "aims_immune/aims_analysis.py", line 1860, in calc_AIMSdist
    if seqSet1.equals(seqSet2):
       ^^^^^^^^^^^^^^
```

**Root Cause:** The function expects pandas DataFrames but the test passes numpy arrays. The `.equals()` method is a DataFrame method.

**Fix Needed:** Pass pandas DataFrames or convert arrays to DataFrames in the test.

---

#### test_calc_aims_dist_identical_matrices ❌
**Location:** aims_immune/test/test_aims_analysis.py:174  
**Status:** FAILED

**Error:**
```
Same as test_calc_aims_dist_basic
AttributeError: 'numpy.ndarray' object has no attribute 'equals'
```

**Root Cause:** Same as above - type mismatch between expected (DataFrame) and provided (ndarray).

**Fix Needed:** Use DataFrames instead of numpy arrays.

---

#### test_calc_aims_dist_different_matrices ❌
**Location:** aims_immune/test/test_aims_analysis.py:182  
**Status:** FAILED

**Error:**
```
Same as test_calc_aims_dist_basic
AttributeError: 'numpy.ndarray' object has no attribute 'equals'
```

**Root Cause:** Same type mismatch issue.

**Fix Needed:** Convert test data to pandas DataFrames.

---

### 5. TestFullAAFreq (2 failures)

#### test_full_aa_freq_basic ❌
**Location:** aims_immune/test/test_aims_analysis.py:192  
**Status:** FAILED

**Error:**
```
TypeError: full_AA_freq() missing 1 required positional argument: 'my_AA_key'
```

**Root Cause:** The function requires `my_AA_key` parameter (amino acid key) but the test doesn't provide it.

**Function Signature:** `full_AA_freq(seqF, my_AA_key, ...)`

**Fix Needed:** Pass the amino acid key parameter:
```python
analysis.full_AA_freq(seqF, analysis.AA_key)
```

---

#### test_full_aa_freq_known_sequences ❌
**Location:** aims_immune/test/test_aims_analysis.py:198  
**Status:** FAILED

**Error:**
```
Same as test_full_aa_freq_basic
TypeError: full_AA_freq() missing 1 required positional argument: 'my_AA_key'
```

**Root Cause:** Same - missing required parameter.

**Fix Needed:** Add `my_AA_key` parameter.

---

### 6. TestGenCloneProps (1 failure)

#### test_gen_clone_props_basic ❌
**Location:** aims_immune/test/test_aims_analysis.py:211  
**Status:** FAILED

**Error:**
```
KeyError: (np.int64(0), np.int64(0))
File "aims_immune/aims_analysis.py", line 365, in gen_clone_props
    if poly_PCA[j,k]==m:
       ^^^^^^^^^^^^^
pandas/_libs/index.pyx:583: in pandas._libs/index.StringObjectEngine._check_type
```

**Root Cause:** The function expects a numpy array but receives a DataFrame. The indexing `poly_PCA[j,k]` fails because DataFrame doesn't support tuple indexing like arrays.

**Fix Needed:** Convert DataFrame to numpy array before passing or adjust test data format.

---

### 7. TestGenDsetProps (1 failure)

#### test_gen_dset_props_basic ❌
**Location:** aims_immune/test/test_aims_analysis.py:225  
**Status:** FAILED

**Error:**
```
KeyError: (np.int64(0), np.int64(0))
File "aims_immune/aims_analysis.py", line 320, in gen_dset_props
    if poly_PCA[j,k]==m:
```

**Root Cause:** Same as TestGenCloneProps - type mismatch between expected numpy array and provided DataFrame.

**Fix Needed:** Convert test data to numpy arrays or adjust function expectations.

---

### 8. TestDoStatistics (1 failure)

#### test_do_statistics_basic ❌
**Location:** aims_immune/test/test_aims_analysis.py:238  
**Status:** FAILED

**Error:**
```
TypeError: the resolved dtypes are not compatible with add.reduce. 
Resolved (dtype('<U1'), dtype('<U1'), dtype('<U2'))
File "numpy/_core/_methods.py", line 132, in _mean
    ret = umr_sum(arr, axis, dtype, out, keepdims, where=where)
```

**Root Cause:** The labels are strings ('A', 'B') which can't be used for numerical median calculation. The function expects numeric data.

**Fix Needed:** Pass numeric labels or appropriate data structure.

---

### 9. TestEncodeMetadata (1 failure)

#### test_encode_meta_basic ❌
**Location:** aims_immune/test/test_aims_analysis.py:250  
**Status:** FAILED

**Error:**
```
AssertionError: assert False
where False = isinstance(  group\n0     A\n1     B\n2     C, 
                          (<class 'numpy.ndarray'>, <class 'list'>))
```

**Root Cause:** The function returns a pandas Series (or DataFrame repr) rather than a numpy array or list.

**Fix Needed:** Adjust assertion to accept Series or check `isinstance(result, (np.ndarray, list, pd.Series))`.

---

### 10. TestCreateMsaPairs (1 failure)

#### test_create_msa_pairs_basic ❌
**Location:** aims_immune/test/test_aims_analysis.py:260  
**Status:** FAILED

**Error:**
```
TypeError: create_msa_pairs() missing 3 required positional arguments: 
'file2', 'pair_list', and 'names'
```

**Root Cause:** The function requires more parameters than the test provides. It needs:
- `seqF` (provided)
- `file2` (missing)
- `pair_list` (missing)
- `names` (missing)

**Function Signature:** `create_msa_pairs(seqF, file2, pair_list, names, ...)`

**Fix Needed:** Either provide all required parameters or check function documentation.

---

### 11. TestIntegration::test_full_analysis_pipeline (1 failure)

#### test_full_analysis_pipeline ❌
**Location:** aims_immune/test/test_aims_analysis.py:294  
**Status:** FAILED

**Error:**
```
IndexError: index 14 is out of bounds for axis 0 with size 14
File "aims_immune/aims_analysis.py", line 836, in gen_peptide_matrix
    pep_PCA[i][count]=key[j]
```

**Root Cause:** Cascading failure from `gen_peptide_matrix` with incompatible data.

**Fix Needed:** Use proper peptide sequence data.

---

## Summary by Root Cause

| Root Cause | Count | Tests |
|-----------|-------|-------|
| **Incorrect test data structure** | 5 | peptide_matrix (2), gen_clone_props, gen_dset_props, full_analysis_pipeline |
| **Missing required parameters** | 3 | full_aa_freq (2), create_msa_pairs |
| **Type mismatch (array vs DataFrame)** | 4 | calc_AIMSdist (3), gen_clone_props* |
| **Wrong return type assumptions** | 1 | gen_tcr_matrix_with_binary_mode |
| **Incorrect input structure** | 2 | get_sequence_dimension (2) |
| **Function API misunderstanding** | 2 | gen_tcr_matrix_with_size_specification, encode_meta |
| **Data type incompatibility** | 1 | do_statistics |

---

## Recommendations for Fixing

### High Priority (Core Logic Issues)

1. **Type Consistency** - Review and standardize whether functions expect arrays or DataFrames
2. **Parameter Documentation** - Clearly document all required parameters in function signatures
3. **Input Validation** - Add type checking to provide clear error messages

### Medium Priority (Test Refinement)

4. **Test Data** - Use real test data from `app_data/test_data/` for all tests
5. **Return Type Checking** - Update assertions to match actual return types
6. **Documentation** - Add docstrings to explain expected input formats

### Low Priority (Enhancement)

7. **Flexibility** - Consider supporting multiple input types (array, DataFrame)
8. **Error Messages** - Improve error messages when input types don't match

---

## Passing Tests (14 ✅)

The following tests pass and validate correct functionality:

✅ TestGetSequenceDimension::test_get_sequence_dimension_returns_lengths  
✅ TestGetSequenceDimension::test_get_sequence_dimension_long_sequences  
✅ TestGetSequenceDimension::test_get_sequence_dimension_variable_lengths  
✅ TestGenTcrMatrix::test_gen_tcr_matrix_basic  
✅ TestGenTcrMatrix::test_gen_tcr_matrix_dimensions  
✅ TestGenTcrMatrix::test_gen_tcr_matrix_alignment_options  
✅ TestGenTcrMatrix::test_gen_tcr_matrix_invalid_bulge_pad  
✅ TestGetProps::test_get_props_returns_dict  
✅ TestGetProps::test_get_props_length  
✅ TestGetProps::test_get_props_first_element  
✅ TestGenMsaMatrix::test_gen_msa_matrix_basic  
✅ TestGenMsaMatrix::test_gen_msa_matrix_with_gaps  
✅ TestIntegration::test_load_and_generate_tcr_matrix  
✅ TestIntegration::test_load_and_get_sequence_dimensions  

---

## Execution Environment

- **Python:** 3.12.13
- **pytest:** 9.1.1
- **Platform:** Darwin (macOS)
- **Conda Environment:** aims_dev

---

## Test Execution Command

```bash
conda activate aims_dev
python -m pytest aims_immune/test/test_aims_analysis.py -v --tb=short
```

---

**Report Generated:** 2026-09-08  
**Test File:** aims_immune/test/test_aims_analysis.py  
**Module Tested:** aims_immune.aims_analysis
