# AIMS-Immune Test Suite

This directory contains a comprehensive test suite for the aims_immune package, covering all major modules and workflows.

## Overview

The test suite includes:

- **test_aims_loader.py** - Tests for sequence loading (FASTA, CSV), MSA subsetting, and 3-letter code conversion
- **test_aims_analysis.py** - Tests for matrix generation, sequence analysis, and distance calculations
- **test_aims_classification.py** - Tests for classification pipelines, feature selection, and machine learning models
- **test_aims_utils.py** - Tests for utility functions and helper operations
- **test_integration.py** - End-to-end integration tests for complete workflows
- **conftest.py** - Pytest fixtures for shared test data and utilities

## Running the Tests

### Run all tests
```bash
conda activate aims_dev
pytest
```

### Run specific test file
```bash
pytest test_aims_loader.py
```

### Run specific test class
```bash
pytest test_aims_loader.py::TestSeqLoader
```

### Run specific test
```bash
pytest test_aims_loader.py::TestSeqLoader::test_seq_loader_csv_basic
```

### Run with verbose output
```bash
pytest -v
```

### Run with coverage report
```bash
pip install pytest-cov
pytest --cov=aims_immune --cov-report=html
```

### Run only integration tests
```bash
pytest test_integration.py
```

### Run tests and show print statements
```bash
pytest -s
```

## Test Organization

### Unit Tests

**aims_loader module:**
- Sequence loading (FASTA and CSV formats)
- Auto-detection of file formats
- Duplicate removal
- Index return functionality
- Missing value handling
- Real data loading

**aims_analysis module:**
- Sequence dimension calculation
- TCR matrix generation
- Peptide matrix generation
- MSA matrix generation
- Amino acid frequency analysis
- Distance calculations
- Statistical analysis

**aims_classification module:**
- Matrix generation and manipulation
- Feature selection
- Classification models (LDA, SVM, Random Forest)
- Cross-validation strategies
- Pretrained model application

### Integration Tests

- End-to-end antibody analysis workflows
- TCR comparison pipelines
- MHC sequence analysis
- Peptide data processing
- Cross-dataset analysis
- Error handling and data consistency

## Test Data

Test data is located in `app_data/test_data/` and includes:

- **abs/** - Antibody sequences (polyreactive and monoreactive)
- **tcrs/** - T-cell receptor sequences
- **mhcs/** - MHC sequences (FASTA format)
- **peptides/** - Peptide sequences

Fixtures provide convenient access to these files:

```python
def test_example(abs_data_paths, tcr_data_paths, mhc_data_paths):
    poly_file = abs_data_paths['poly']
    tcr_file = tcr_data_paths['siv_cm9']
    mhc_file = mhc_data_paths['cd1_fasta']
```

## Fixtures

Common fixtures available in conftest.py:

- `test_data_dir` - Path to test data directory
- `abs_data_paths` - Dictionary of antibody file paths
- `tcr_data_paths` - Dictionary of TCR file paths
- `mhc_data_paths` - Dictionary of MHC file paths
- `peptide_data_paths` - Dictionary of peptide file paths
- `sample_dataframe` - Sample DataFrame for testing
- `temp_fasta_file` - Temporary FASTA file for testing
- `temp_csv_file` - Temporary CSV file for testing

## Test Coverage

The test suite covers:

- ✅ Core loading functions
- ✅ Matrix generation algorithms
- ✅ Classification pipelines
- ✅ Data validation and error handling
- ✅ End-to-end workflows
- ✅ Cross-dataset analysis
- ✅ Different molecule types (antibodies, TCRs, MHCs, peptides)

## Adding New Tests

When adding new tests:

1. Create a test class for logical grouping
2. Use descriptive test names starting with `test_`
3. Add docstrings explaining what is being tested
4. Use appropriate fixtures for data access
5. Include both positive and edge case tests

Example:

```python
class TestNewFeature:
    """Tests for new feature"""

    def test_new_feature_basic(self, fixture_name):
        """Test basic functionality"""
        result = function_to_test(fixture_name)
        assert result is not None

    def test_new_feature_edge_case(self):
        """Test edge case"""
        with pytest.raises(ValueError):
            function_to_test(invalid_input)
```

## Dependencies

Required testing dependencies:
- pytest
- numpy
- pandas
- scipy
- scikit-learn
- biopython

Install with:
```bash
pip install pytest
```

## Troubleshooting

### Import errors
Ensure the aims_dev conda environment is activated:
```bash
conda activate aims_dev
```

### Missing test data
Test data should be in `app_data/test_data/`. If missing, run:
```bash
python -m aims_immune.get_tests
```

### Slow tests
Some integration tests may be slow. To skip them:
```bash
pytest -k "not integration"
```

### Memory issues
If running into memory issues with large datasets, run individual test files:
```bash
pytest test_aims_loader.py --maxfail=1
```
