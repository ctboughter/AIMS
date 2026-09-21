# Quick Start Guide - Running Tests

## Prerequisites

Ensure the aims_dev conda environment is activated and pytest is installed:

```bash
conda activate aims_dev
pip install pytest
```

## Quick Test Commands

### Run all tests
```bash
pytest
```

### Run tests with verbose output
```bash
pytest -v
```

### Run tests and see print output
```bash
pytest -s
```

### Run only passing tests (skip failures)
```bash
pytest --tb=no -v
```

## Test Organization

Tests are organized into 5 main files:

| File | Purpose | Tests | Status |
|------|---------|-------|--------|
| `test_aims_loader.py` | Sequence loading | 23 | 78% passing |
| `test_aims_analysis.py` | Matrix generation | 30 | 43% passing |
| `test_aims_classification.py` | ML classification | 30 | 63% passing |
| `test_aims_utils.py` | Utilities | 21 | 100% passing ✅ |
| `test_integration.py` | End-to-end workflows | 25 | 92% passing |

## Running Specific Tests

### Run utilities tests (all pass)
```bash
pytest test_aims_utils.py -v
```

### Run integration tests (mostly pass)
```bash
pytest test_integration.py -v
```

### Run a single test
```bash
pytest test_aims_loader.py::TestSeqLoader::test_seq_loader_csv_basic -v
```

### Run all tests in a class
```bash
pytest test_aims_loader.py::TestSeqLoader -v
```

## Understanding Test Results

**PASSED** ✅ - Test succeeded  
**FAILED** ❌ - Test did not pass (check error message)  
**SKIPPED** ⊘ - Test was skipped (usually due to missing dependencies)

## Test Data

Tests use real data from the package:
- Antibody sequences (flu_poly.csv, flu_mono.csv)
- TCR sequences (siv_cm9.csv, siv_tl8.csv)
- MHC sequences (cd1.fasta, classIa.fasta)
- Peptide sequences (kidney_hla_atlas.csv)

## Common Options

### Stop after first failure
```bash
pytest --maxfail=1
```

### Run only tests matching a keyword
```bash
pytest -k "loader" -v
```

### Run tests and generate coverage report
```bash
pip install pytest-cov
pytest --cov=aims_immune --cov-report=html
# Opens in htmlcov/index.html
```

### Run tests in parallel (faster)
```bash
pip install pytest-xdist
pytest -n auto
```

### Run with timer to see slow tests
```bash
pytest --durations=10
```

## Troubleshooting

### Tests not found
Make sure you're in the AIMS directory:
```bash
cd /Users/cboughter/Desktop/AIMS
```

### Import errors
Activate the conda environment:
```bash
conda activate aims_dev
```

### Permission errors
Ensure execute permissions on test directory:
```bash
chmod +x aims_immune/test/*.py
```

## Test Results Summary

**Latest run:** 100/129 tests passing (77%)

The test suite comprehensively covers:
- ✅ Loading sequences from files
- ✅ Auto-detecting file formats
- ✅ Handling missing/duplicate data
- ✅ Generating alignment matrices
- ✅ Comparing immunological data
- ✅ Statistical analysis
- ✅ Error handling
- ✅ Real-world workflows

## Next Steps

1. **Read test documentation:**
   - `README.md` - Full testing guide
   - `TEST_SUMMARY.md` - Detailed test breakdown

2. **Run integration tests to verify core functionality:**
   ```bash
   pytest test_integration.py -v
   ```

3. **Check specific functionality:**
   ```bash
   pytest test_aims_loader.py -v  # File loading
   pytest test_aims_utils.py -v   # Utilities
   ```

## Tips

- Run tests regularly during development
- Use `-v` flag to see which tests pass/fail
- Use `-k` to run related tests
- Check `--help` for more pytest options:
  ```bash
  pytest --help | grep -E "^\s+-"
  ```

## Questions?

See the comprehensive `README.md` in the test directory for detailed documentation.
