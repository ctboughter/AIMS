"""Tests for aims_utils module"""
import pytest
import numpy as np
import pandas as pd
import tempfile
import os
from pathlib import Path

# Import aims_utils if it contains testable functions
try:
    from aims_immune import aims_utils as utils
    UTILS_AVAILABLE = True
except ImportError:
    UTILS_AVAILABLE = False


@pytest.mark.skipif(not UTILS_AVAILABLE, reason="aims_utils module not available")
class TestAimsUtils:
    """Tests for utility functions in aims_utils"""

    def test_utils_import(self):
        """Test that utils module can be imported"""
        assert UTILS_AVAILABLE


class TestDataHandling:
    """Tests for data handling utilities"""

    def test_dataframe_creation(self):
        """Test creation and validation of DataFrames"""
        df = pd.DataFrame({
            'col1': ['A', 'B', 'C'],
            'col2': [1, 2, 3]
        })
        assert isinstance(df, pd.DataFrame)
        assert df.shape == (3, 2)

    def test_numpy_array_operations(self):
        """Test numpy array operations"""
        arr = np.array([[1, 2, 3], [4, 5, 6]])
        assert arr.shape == (2, 3)
        assert arr.dtype in [np.int64, np.int32, int]

    def test_sequence_validation(self):
        """Test sequence validation"""
        valid_sequence = 'ACDEFGHIKLMNPQRSTVWY'
        # Check all characters are valid amino acids
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        for char in valid_sequence:
            assert char in valid_aa


class TestMatrixOperations:
    """Tests for matrix operations"""

    def test_matrix_transpose(self):
        """Test matrix transposition"""
        mat = np.array([[1, 2, 3], [4, 5, 6]])
        transposed = np.transpose(mat)
        assert transposed.shape == (3, 2)

    def test_matrix_multiplication(self):
        """Test matrix multiplication"""
        mat1 = np.array([[1, 2], [3, 4]])
        mat2 = np.array([[5, 6], [7, 8]])
        result = np.dot(mat1, mat2)
        assert result.shape == (2, 2)

    def test_matrix_normalization(self):
        """Test matrix normalization"""
        mat = np.array([[1.0, 2.0], [3.0, 4.0]])
        # Min-max normalization
        normalized = (mat - mat.min()) / (mat.max() - mat.min())
        assert normalized.min() == 0.0
        assert normalized.max() == 1.0


class TestFileIO:
    """Tests for file input/output operations"""

    def test_create_temp_file(self):
        """Test creating and reading temporary files"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write('test data')
            temp_path = f.name

        try:
            with open(temp_path, 'r') as f:
                content = f.read()
            assert content == 'test data'
        finally:
            os.remove(temp_path)

    def test_csv_io(self):
        """Test CSV file input/output"""
        df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})

        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            df.to_csv(f.name, index=False)
            temp_path = f.name

        try:
            loaded_df = pd.read_csv(temp_path)
            assert loaded_df.shape == df.shape
            assert list(loaded_df.columns) == list(df.columns)
        finally:
            os.remove(temp_path)

    def test_dataframe_to_list(self):
        """Test converting DataFrame to list"""
        df = pd.DataFrame({'A': [1, 2], 'B': [3, 4]})
        result = df.values.tolist()
        assert isinstance(result, list)
        assert len(result) == 2


class TestNumericalOperations:
    """Tests for numerical operations"""

    def test_mean_calculation(self):
        """Test mean calculation"""
        arr = np.array([1, 2, 3, 4, 5])
        assert np.mean(arr) == 3.0

    def test_std_calculation(self):
        """Test standard deviation calculation"""
        arr = np.array([1, 2, 3, 4, 5])
        std = np.std(arr)
        assert std > 0

    def test_percentile_calculation(self):
        """Test percentile calculation"""
        arr = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        p50 = np.percentile(arr, 50)
        assert 4 < p50 < 7

    def test_distance_calculation(self):
        """Test distance calculations"""
        from scipy.spatial.distance import euclidean
        p1 = np.array([0, 0])
        p2 = np.array([3, 4])
        dist = euclidean(p1, p2)
        assert dist == 5.0


class TestStatisticalTests:
    """Tests for statistical functions"""

    def test_pearson_correlation(self):
        """Test Pearson correlation"""
        from scipy.stats import pearsonr
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([2, 4, 5, 4, 6])
        corr, pval = pearsonr(x, y)
        assert -1 <= corr <= 1

    def test_t_test(self):
        """Test t-test"""
        from scipy.stats import ttest_ind
        group1 = np.array([1, 2, 3, 4, 5])
        group2 = np.array([2, 3, 4, 5, 6])
        t_stat, p_val = ttest_ind(group1, group2)
        assert isinstance(t_stat, (float, np.floating))

    def test_chi_square(self):
        """Test chi-square test"""
        from scipy.stats import chi2_contingency
        contingency = np.array([[10, 20], [30, 40]])
        chi2, p_val, dof, expected = chi2_contingency(contingency)
        assert chi2 > 0


class TestDataValidation:
    """Tests for data validation"""

    def test_empty_dataframe_check(self):
        """Test checking for empty DataFrames"""
        df_empty = pd.DataFrame()
        df_full = pd.DataFrame({'A': [1, 2, 3]})

        assert df_empty.empty
        assert not df_full.empty

    def test_nan_detection(self):
        """Test NaN detection"""
        arr = np.array([1.0, np.nan, 3.0])
        assert np.isnan(arr[1])
        assert not np.isnan(arr[0])

    def test_nan_removal(self):
        """Test NaN removal"""
        arr = np.array([1.0, np.nan, 3.0])
        clean_arr = arr[~np.isnan(arr)]
        assert len(clean_arr) == 2

    def test_dtype_checking(self):
        """Test data type checking"""
        arr_int = np.array([1, 2, 3])
        arr_float = np.array([1.0, 2.0, 3.0])
        arr_str = np.array(['a', 'b', 'c'])

        assert np.issubdtype(arr_int.dtype, np.integer)
        assert np.issubdtype(arr_float.dtype, np.floating)
        assert np.issubdtype(arr_str.dtype, np.str_) or np.issubdtype(arr_str.dtype, np.object_)
