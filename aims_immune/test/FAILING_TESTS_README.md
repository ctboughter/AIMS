# Failing Tests Output - Complete Reference

## Quick Summary

**File Analyzed:** `aims_immune/test/test_aims_analysis.py`  
**Total Tests:** 31  
**Passing:** 14 ✅  
**Failing:** 17 ❌  
**Success Rate:** 45%  

---

## Output Files Created

This analysis includes 3 detailed output files:

1. **FAILING_TESTS_ANALYSIS.md** - Detailed analysis of each failing test
   - Complete error messages
   - Root cause analysis
   - Organized by test class
   
2. **FIX_RECOMMENDATIONS.md** - Actionable fixes for each failure
   - Code examples showing wrong vs. right approaches
   - Implementation priority (Phase 1, 2, 3)
   - Time estimates
   
3. **failing_tests_summary.csv** - Quick reference table
   - Test name and class
   - Error type
   - Severity level
   - Fix required

4. **FAILING_TESTS_README.md** - This file
   - Navigation guide
   - How to use the output files
   - Quick reference tables

---

## How to Use These Files

### For Quick Overview
**Read:** `failing_tests_summary.csv`
- 2-minute read
- Spreadsheet-friendly format
- All 17 failures in one table

### For Detailed Understanding
**Read:** `FAILING_TESTS_ANALYSIS.md`
- 15-minute read
- Full error messages and stack traces
- Root cause analysis for each failure
- Passing tests listed

### For Implementation
**Read:** `FIX_RECOMMENDATIONS.md`
- 20-minute read
- Code examples for each fix
- Step-by-step implementation plan
- Phased approach (3 phases)

---

## Failure Categories Summary

### By Error Type

| Error Type | Count | Tests |
|-----------|-------|-------|
| **IndexError** | 4 | peptide_matrix (2), gen_tcr_size, full_pipeline |
| **TypeError** | 4 | full_aa_freq (2), create_msa_pairs, do_statistics |
| **AttributeError** | 3 | calc_aims_dist (3) |
| **KeyError** | 3 | gen_clone_props, gen_dset_props, encode_meta |
| **AssertionError** | 3 | get_sequence_dimension (2), gen_tcr_binary |

### By Severity

| Severity | Count | Examples |
|----------|-------|----------|
| **HIGH** | 8 | Type mismatches, incompatible data, missing parameters |
| **MEDIUM** | 5 | Test logic errors, API misunderstandings |
| **LOW** | 4 | Return type assertions, minor fixes |

### By Root Cause

| Root Cause | Count | Typical Fix Time |
|-----------|-------|-----------------|
| Type mismatch (array vs DataFrame) | 6 | 15 min per test |
| Missing parameters | 2 | 10 min per test |
| Incompatible test data | 2 | 30-45 min per test |
| Wrong return type assumptions | 2 | 5 min per test |
| API misunderstandings | 3 | 10-20 min per test |
| Incorrect test logic | 2 | 10 min per test |

---

## Failing Tests at a Glance

### TestGetSequenceDimension (2 failures)
❌ test_get_sequence_dimension_basic  
❌ test_get_sequence_dimension_single_column  
**Issue:** Test assumes wrong input structure  
**Fix Time:** 10 min

### TestGenTcrMatrix (2 failures)
❌ test_gen_tcr_matrix_with_size_specification  
❌ test_gen_tcr_matrix_with_binary_mode  
**Issue:** Wrong parameter type, wrong return type handling  
**Fix Time:** 15 min

### TestGenPeptideMatrix (2 failures)
❌ test_gen_peptide_matrix_basic  
❌ test_gen_peptide_matrix_shape  
**Issue:** Test data incompatible with peptide format  
**Fix Time:** 45 min

### TestCalcAimsDist (3 failures)
❌ test_calc_aims_dist_basic  
❌ test_calc_aims_dist_identical_matrices  
❌ test_calc_aims_dist_different_matrices  
**Issue:** Numpy arrays instead of DataFrames  
**Fix Time:** 20 min

### TestFullAAFreq (2 failures)
❌ test_full_aa_freq_basic  
❌ test_full_aa_freq_known_sequences  
**Issue:** Missing required parameter  
**Fix Time:** 10 min

### TestGenCloneProps (1 failure)
❌ test_gen_clone_props_basic  
**Issue:** DataFrame indexing instead of array indexing  
**Fix Time:** 15 min

### TestGenDsetProps (1 failure)
❌ test_gen_dset_props_basic  
**Issue:** DataFrame indexing instead of array indexing  
**Fix Time:** 15 min

### TestDoStatistics (1 failure)
❌ test_do_statistics_basic  
**Issue:** String labels incompatible with numeric operations  
**Fix Time:** 10 min

### TestEncodeMetadata (1 failure)
❌ test_encode_meta_basic  
**Issue:** Returns Series not array/list  
**Fix Time:** 5 min

### TestCreateMsaPairs (1 failure)
❌ test_create_msa_pairs_basic  
**Issue:** Missing required parameters  
**Fix Time:** 20 min

### TestIntegration (1 failure)
❌ test_full_analysis_pipeline  
**Issue:** Cascading from gen_peptide_matrix  
**Fix Time:** Resolved by fixing gen_peptide_matrix

---

## Recommended Reading Order

### For Developers Fixing Tests
1. Start: `failing_tests_summary.csv` (2 min)
2. Deep dive: `FAILING_TESTS_ANALYSIS.md` - Read relevant test section (10 min)
3. Implementation: `FIX_RECOMMENDATIONS.md` - Find your test class (5 min)
4. Code: Review the code examples in FIX_RECOMMENDATIONS.md (10 min)

**Total Time:** 25-30 minutes per test to understand and fix

### For Project Managers
1. `FAILING_TESTS_ANALYSIS.md` - Read "Summary by Root Cause" section (3 min)
2. `FIX_RECOMMENDATIONS.md` - Read "Implementation Priority" section (2 min)
3. `failing_tests_summary.csv` - Review severity column (1 min)

**Total Time:** 5-10 minutes to understand scope

### For Team Leads
1. `FIX_RECOMMENDATIONS.md` - Read all (20 min)
2. `FAILING_TESTS_ANALYSIS.md` - Skim for understanding (10 min)
3. Plan implementation phases (15 min)

**Total Time:** 45 minutes for full understanding

---

## Key Metrics

### Test Status Breakdown
```
Passing Tests:     14 ✅ (45%)
├── Sequence Dimensions:  3 tests
├── TCR Matrix:           4 tests  
├── Properties:           3 tests
├── MSA Matrix:           2 tests
└── Integration:          2 tests

Failing Tests:     17 ❌ (55%)
├── High Priority:   8 tests (2-3 hours)
├── Medium Priority: 5 tests (1-2 hours)
└── Low Priority:    4 tests (15-30 min)

Total Fix Time:    4-5.5 hours for all
```

### By Test Class

| Test Class | Total | Passing | Failing | Pass Rate |
|-----------|-------|---------|---------|-----------|
| TestGetSequenceDimension | 5 | 3 | 2 | 60% |
| TestGenTcrMatrix | 6 | 4 | 2 | 67% |
| TestGetProps | 3 | 3 | 0 | 100% ✅ |
| TestGenPeptideMatrix | 2 | 0 | 2 | 0% |
| TestGenMsaMatrix | 2 | 2 | 0 | 100% ✅ |
| TestCalcAimsDist | 3 | 0 | 3 | 0% |
| TestFullAAFreq | 2 | 0 | 2 | 0% |
| TestGenCloneProps | 1 | 0 | 1 | 0% |
| TestGenDsetProps | 1 | 0 | 1 | 0% |
| TestDoStatistics | 1 | 0 | 1 | 0% |
| TestEncodeMetadata | 1 | 0 | 1 | 0% |
| TestCreateMsaPairs | 1 | 0 | 1 | 0% |
| TestIntegration | 3 | 2 | 1 | 67% |

---

## Quick Commands

### Run Only Failing Tests
```bash
pytest test_aims_analysis.py -k "test_get_sequence_dimension_basic or test_gen_tcr_matrix_with_size"
```

### Run Tests by Class
```bash
pytest test_aims_analysis.py::TestCalcAimsDist -v
```

### Run with Output Capture
```bash
pytest test_aims_analysis.py -v -s --tb=short
```

### Generate Comparison Report
```bash
# Before fixes
pytest test_aims_analysis.py -v > before_fixes.txt

# After fixes
pytest test_aims_analysis.py -v > after_fixes.txt

# Compare
diff before_fixes.txt after_fixes.txt
```

---

## Next Steps

1. **Review** the summary files (this doc + summary.csv)
2. **Understand** root causes (FAILING_TESTS_ANALYSIS.md)
3. **Plan** fixes (FIX_RECOMMENDATIONS.md)
4. **Implement** Phase 1 fixes (highest priority)
5. **Test** and verify
6. **Move** to Phase 2 and 3

---

## File Locations

All output files are in:
```
/Users/cboughter/Desktop/AIMS/aims_immune/test/
```

| File | Size | Purpose |
|------|------|---------|
| FAILING_TESTS_ANALYSIS.md | ~8 KB | Detailed analysis |
| FIX_RECOMMENDATIONS.md | ~6 KB | Implementation guide |
| failing_tests_summary.csv | ~2 KB | Quick reference |
| FAILING_TESTS_README.md | ~7 KB | This file |

---

## Additional Resources

- **Test Results:** Run `pytest test_aims_analysis.py -v --tb=short`
- **Function Docs:** Check docstrings in `aims_analysis.py`
- **Test Data:** Located in `app_data/test_data/`
- **Fixtures:** Defined in `conftest.py`

---

## Questions?

For specific test failures:
1. Find test in `failing_tests_summary.csv`
2. Look up detailed error in `FAILING_TESTS_ANALYSIS.md`
3. Find fix in `FIX_RECOMMENDATIONS.md`
4. Review code example
5. Implement and test

---

**Generated:** 2026-09-08  
**Module:** aims_immune.aims_analysis  
**Test File:** test_aims_analysis.py  
**Environment:** Python 3.12.13, pytest 9.1.1
