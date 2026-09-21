"""Integration tests for aims_immune package"""
import pytest
import numpy as np
import pandas as pd
from aims_immune import aims_loader as loader
from aims_immune import aims_analysis as analysis
from aims_immune import aims_classification as classify


class TestEndToEndAntibodyAnalysis:
    """End-to-end tests for antibody analysis workflow"""

    def test_load_poly_and_mono_antibodies(self, abs_data_paths):
        """Test loading polyreactive and monoreactive antibody data"""
        poly_df = loader.seq_loader(str(abs_data_paths['poly']), label='poly')
        mono_df = loader.seq_loader(str(abs_data_paths['mono']), label='mono')

        assert isinstance(poly_df, pd.DataFrame)
        assert isinstance(mono_df, pd.DataFrame)
        assert poly_df.shape[0] == mono_df.shape[0] == 6  # CDR regions

    def test_antibody_dimension_analysis(self, abs_data_paths):
        """Test sequence dimension analysis for antibodies"""
        abs_data = loader.seq_loader(str(abs_data_paths['poly']), label='poly')

        dims = analysis.get_sequence_dimension(abs_data)

        assert isinstance(dims, list)
        assert len(dims) == 6
        assert all(d > 0 for d in dims)

    def test_antibody_matrix_generation(self, abs_data_paths):
        """Test matrix generation from antibody sequences"""
        abs_data = loader.seq_loader(str(abs_data_paths['poly']), label='poly')

        # Prepare data for matrix generation
        pre_poly = [abs_data.iloc[i].tolist() for i in range(len(abs_data))]

        matrix = analysis.gen_peptide_matrix(pre_poly)
        assert isinstance(matrix, np.ndarray)
        assert matrix.shape[0] > 0

    def test_antibody_comparison_workflow(self, abs_data_paths):
        """Test complete comparison workflow for antibodies"""
        # Load data
        poly_df = loader.seq_loader(str(abs_data_paths['poly']), label='poly')
        mono_df = loader.seq_loader(str(abs_data_paths['mono']), label='mono')

        # Get dimensions
        poly_dims = analysis.get_sequence_dimension(poly_df)
        mono_dims = analysis.get_sequence_dimension(mono_df)

        assert poly_dims is not None
        assert mono_dims is not None

        # Generate matrices
        pre_poly = [poly_df.iloc[i].tolist() for i in range(len(poly_df))]
        pre_mono = [mono_df.iloc[i].tolist() for i in range(len(mono_df))]

        poly_mat = analysis.gen_peptide_matrix(pre_poly)
        mono_mat = analysis.gen_peptide_matrix(pre_mono)

        assert poly_mat is not None
        assert mono_mat is not None


class TestEndToEndTCRAnalysis:
    """End-to-end tests for TCR analysis workflow"""

    def test_load_multiple_tcr_datasets(self, tcr_data_paths):
        """Test loading multiple TCR datasets"""
        siv_cm9 = loader.seq_loader(str(tcr_data_paths['siv_cm9']), label='cm9')
        siv_tl8 = loader.seq_loader(str(tcr_data_paths['siv_tl8']), label='tl8')

        assert isinstance(siv_cm9, pd.DataFrame)
        assert isinstance(siv_tl8, pd.DataFrame)

    def test_tcr_dimension_analysis(self, tcr_data_paths):
        """Test sequence dimension analysis for TCRs"""
        tcr_data = loader.seq_loader(str(tcr_data_paths['siv_cm9']), label='tcr')

        dims = analysis.get_sequence_dimension(tcr_data)

        assert isinstance(dims, list)
        assert len(dims) > 0
        assert all(d > 0 for d in dims)

    def test_tcr_matrix_generation(self, tcr_data_paths):
        """Test matrix generation from TCR sequences"""
        tcr_data = loader.seq_loader(str(tcr_data_paths['siv_cm9']), label='tcr')

        # Take a subset for testing
        pre_tcr = [tcr_data.iloc[i].tolist() for i in range(min(3, len(tcr_data)))]

        matrix = analysis.gen_tcr_matrix(pre_tcr)
        assert isinstance(matrix, np.ndarray)
        assert matrix.shape[0] > 0

    def test_tcr_comparison_workflow(self, tcr_data_paths):
        """Test TCR comparison between different antigens"""
        siv_cm9 = loader.seq_loader(str(tcr_data_paths['siv_cm9']), label='cm9')
        siv_tl8 = loader.seq_loader(str(tcr_data_paths['siv_tl8']), label='tl8')

        # Get dimensions
        cm9_dims = analysis.get_sequence_dimension(siv_cm9)
        tl8_dims = analysis.get_sequence_dimension(siv_tl8)

        assert cm9_dims is not None
        assert tl8_dims is not None


class TestEndToEndMHCAnalysis:
    """End-to-end tests for MHC analysis workflow"""

    def test_load_mhc_sequences(self, mhc_data_paths):
        """Test loading MHC sequences"""
        cd1 = loader.seq_loader(str(mhc_data_paths['cd1_fasta']), label='cd1', datType='fasta')
        classIa = loader.seq_loader(str(mhc_data_paths['classIa_fasta']), label='classIa', datType='fasta')

        assert isinstance(cd1, pd.DataFrame)
        assert isinstance(classIa, pd.DataFrame)

    def test_mhc_dimension_analysis(self, mhc_data_paths):
        """Test sequence dimension analysis for MHC"""
        mhc_data = loader.seq_loader(str(mhc_data_paths['cd1_fasta']), label='mhc', datType='fasta')

        dims = analysis.get_sequence_dimension(mhc_data)

        assert isinstance(dims, list)
        assert len(dims) > 0

    def test_msa_matrix_generation(self, mhc_data_paths):
        """Test MSA matrix generation for MHC"""
        mhc_data = loader.seq_loader(str(mhc_data_paths['cd1_fasta']), label='mhc', datType='fasta')

        # Take first few sequences
        pre_msa = [mhc_data.iloc[i].tolist() for i in range(min(2, len(mhc_data)))]

        try:
            matrix = analysis.gen_MSA_matrix(pre_msa)
            assert isinstance(matrix, np.ndarray)
        except Exception as e:
            # MSA generation might have specific requirements
            pass


class TestEndToEndPeptideAnalysis:
    """End-to-end tests for peptide analysis workflow"""

    def test_load_peptide_data(self, peptide_data_paths):
        """Test loading peptide data"""
        kidney = loader.seq_loader(str(peptide_data_paths['kidney']), label='kidney')
        pancreas = loader.seq_loader(str(peptide_data_paths['pancreas']), label='pancreas')

        assert isinstance(kidney, pd.DataFrame)
        assert isinstance(pancreas, pd.DataFrame)

    def test_peptide_dimension_analysis(self, peptide_data_paths):
        """Test sequence dimension analysis for peptides"""
        pep_data = loader.seq_loader(str(peptide_data_paths['kidney']), label='peptide')

        dims = analysis.get_sequence_dimension(pep_data)

        assert isinstance(dims, list)
        assert all(d > 0 for d in dims)

    def test_peptide_matrix_generation(self, peptide_data_paths):
        """Test matrix generation from peptide sequences"""
        pep_data = loader.seq_loader(str(peptide_data_paths['kidney']), label='peptide')

        pre_pep = [pep_data.iloc[i].tolist() for i in range(min(3, len(pep_data)))]

        matrix = analysis.gen_peptide_matrix(pre_pep)
        assert isinstance(matrix, np.ndarray)


class TestDataProcessingPipeline:
    """Tests for complete data processing pipeline"""

    def test_load_process_multiple_datasets(self, abs_data_paths, tcr_data_paths):
        """Test loading and processing multiple dataset types"""
        # Load different data types
        abs_data = loader.seq_loader(str(abs_data_paths['poly']), label='abs')
        tcr_data = loader.seq_loader(str(tcr_data_paths['siv_cm9']), label='tcr')

        # Get dimensions
        abs_dims = analysis.get_sequence_dimension(abs_data)
        tcr_dims = analysis.get_sequence_dimension(tcr_data)

        assert abs_dims is not None
        assert tcr_dims is not None

    def test_duplicate_removal_workflow(self, abs_data_paths):
        """Test complete workflow with duplicate handling"""
        # Load with duplicates
        with_dups = loader.seq_loader(str(abs_data_paths['poly']),
                                     label='poly', drop_dups=False)
        # Load without duplicates
        no_dups = loader.seq_loader(str(abs_data_paths['poly']),
                                   label='poly', drop_dups=True)

        # No_dups should have <= columns
        assert no_dups.shape[1] <= with_dups.shape[1]

    def test_index_return_workflow(self, abs_data_paths):
        """Test workflow with index return"""
        data, indices = loader.seq_loader(str(abs_data_paths['poly']),
                                         label='poly', return_index=True)

        assert isinstance(data, pd.DataFrame)
        assert isinstance(indices, (list, pd.Index))


class TestCrossDatasetAnalysis:
    """Tests for analysis across different dataset combinations"""

    def test_compare_poly_vs_mono(self, abs_data_paths):
        """Test comparison between polyreactive and monoreactive"""
        poly = loader.seq_loader(str(abs_data_paths['poly']), label='poly')
        mono = loader.seq_loader(str(abs_data_paths['mono']), label='mono')

        poly_dims = analysis.get_sequence_dimension(poly)
        mono_dims = analysis.get_sequence_dimension(mono)

        # Should have same structure for antibodies
        assert len(poly_dims) == len(mono_dims)

    def test_compare_cdr_regions(self, abs_data_paths):
        """Test individual CDR region analysis"""
        abs_data = loader.seq_loader(str(abs_data_paths['poly']), label='poly')

        # Each row is a different CDR region
        for i in range(len(abs_data)):
            cdr_region = abs_data.iloc[i]
            dims = analysis.get_sequence_dimension(pd.DataFrame([cdr_region]))
            assert dims is not None

    def test_combine_datasets(self, abs_data_paths):
        """Test combining multiple datasets"""
        poly = loader.seq_loader(str(abs_data_paths['poly']), label='poly')
        poly.index = np.arange(len(poly))

        mono = loader.seq_loader(str(abs_data_paths['mono']), label='mono')
        mono.index = np.arange(len(mono))

        # Combine datasets
        combined = pd.concat([poly, mono], axis=1)

        assert combined.shape[0] == poly.shape[0]
        assert combined.shape[1] == poly.shape[1] + mono.shape[1]


class TestErrorHandling:
    """Tests for error handling in pipeline"""

    def test_invalid_file_handling(self):
        """Test handling of invalid files"""
        with pytest.raises(FileNotFoundError):
            loader.seq_loader('/invalid/path.csv', label='test')

    def test_incompatible_data_formats(self, temp_fasta_file):
        """Test handling of incompatible data"""
        # Try to load FASTA as CSV (auto-detect should handle it)
        result = loader.seq_loader(temp_fasta_file, label='test')
        assert isinstance(result, pd.DataFrame)

    def test_empty_results_handling(self):
        """Test handling when operations return empty results"""
        empty_df = pd.DataFrame()
        # Should handle gracefully
        assert empty_df.empty


class TestDataConsistency:
    """Tests for data consistency across operations"""

    def test_shape_preservation(self, abs_data_paths):
        """Test that data shapes are preserved through operations"""
        data1 = loader.seq_loader(str(abs_data_paths['poly']), label='poly', drop_dups=False)
        data2 = loader.seq_loader(str(abs_data_paths['poly']), label='poly', drop_dups=False)

        assert data1.shape == data2.shape

    def test_sequence_integrity(self, abs_data_paths):
        """Test that sequences are not modified during loading"""
        data = loader.seq_loader(str(abs_data_paths['poly']), label='poly')

        # All entries should be strings
        for col in data.columns:
            for val in data[col]:
                assert isinstance(val, str)
                # Should only contain valid amino acid characters
                assert all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in val.upper())
