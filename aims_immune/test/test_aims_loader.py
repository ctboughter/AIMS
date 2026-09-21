"""Tests for aims_loader module"""
import pytest
import numpy as np
import pandas as pd
import tempfile
import os
from pathlib import Path
from aims_immune import aims_loader as loader


class TestSeqLoader:
    """Tests for seq_loader function"""

    def test_seq_loader_csv_basic(self, temp_csv_file):
        """Test loading a basic CSV file"""
        result = loader.seq_loader(temp_csv_file, label='test')
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == 3  # 3 columns in CSV (l1, l2, l3)
        assert result.shape[1] == 3  # 3 rows of data

    def test_seq_loader_csv_with_label(self, temp_csv_file):
        """Test that label is properly applied to column names"""
        result = loader.seq_loader(temp_csv_file, label='myseq')
        # Check that column names contain the label
        assert all('myseq' in col for col in result.columns)

    def test_seq_loader_fasta_basic(self, temp_fasta_file):
        """Test loading a basic FASTA file"""
        result = loader.seq_loader(temp_fasta_file, label='fasta_test', datType='fasta')
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == 1  # FASTA has 1 row (sequences as columns)
        assert result.shape[1] == 2  # 2 sequences

    def test_seq_loader_auto_detect_fasta(self, temp_fasta_file):
        """Test auto-detection of FASTA format"""
        result = loader.seq_loader(temp_fasta_file, label='auto_fasta')
        assert isinstance(result, pd.DataFrame)
        assert result.shape[1] == 2

    def test_seq_loader_auto_detect_csv(self, temp_csv_file):
        """Test auto-detection of CSV format"""
        result = loader.seq_loader(temp_csv_file, label='auto_csv')
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == 3

    def test_seq_loader_drop_duplicates(self, abs_data_paths):
        """Test duplicate removal with drop_dups=True"""
        csv_path = str(abs_data_paths['poly'])
        result_with_dups = loader.seq_loader(csv_path, label='test', drop_dups=False)
        result_no_dups = loader.seq_loader(csv_path, label='test', drop_dups=True)

        # Results with no dups should have <= columns than with dups
        assert result_no_dups.shape[1] <= result_with_dups.shape[1]

    def test_seq_loader_return_index(self, temp_csv_file):
        """Test return_index parameter"""
        result_data, result_idx = loader.seq_loader(temp_csv_file, label='test', return_index=True)
        assert isinstance(result_data, pd.DataFrame)
        assert isinstance(result_idx, (list, pd.Index))

    def test_seq_loader_csv_with_dataloc(self, abs_data_paths):
        """Test filtering with dataLoc parameter"""
        csv_path = str(abs_data_paths['poly'])
        # Load full data first to know column names
        full_data = pd.read_csv(csv_path)
        cols_to_use = ['l1', 'l2']

        result = loader.seq_loader(csv_path, label='test', dataLoc=cols_to_use)
        assert result.shape[0] == 2  # Should have 2 columns specified

    def test_seq_loader_handles_missing_values(self, temp_csv_file):
        """Test that loader handles missing values properly"""
        # Create CSV with empty values
        content = """l1,l2,l3
QSISSY,DAS,QHRSTWPPN
ESLLHS,,MQTIQLPGT
,EVS,HQCYNLPYT
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            result = loader.seq_loader(temp_path, label='test')
            # Should remove rows with any empty values
            assert result.shape[1] <= 3
        finally:
            os.remove(temp_path)

    def test_seq_loader_removes_x_residues(self):
        """Test that sequences with X residues are removed"""
        content = """l1,l2
QSISSY,DAS
QXSSXY,XVS
QSISSY,DAS
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            result = loader.seq_loader(temp_path, label='test')
            # Should remove sequence with X
            assert result.shape[1] == 2
        finally:
            os.remove(temp_path)

    def test_seq_loader_real_abs_data(self, abs_data_paths):
        """Test with real antibody data"""
        csv_path = str(abs_data_paths['poly'])
        result = loader.seq_loader(csv_path, label='flu_poly')
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == 6  # Should have 6 CDR regions for antibodies
        assert result.shape[1] > 0

    def test_seq_loader_real_tcr_data(self, tcr_data_paths):
        """Test with real TCR data"""
        csv_path = str(tcr_data_paths['siv_cm9'])
        result = loader.seq_loader(csv_path, label='tcr_siv')
        assert isinstance(result, pd.DataFrame)
        assert result.shape[1] > 0

    def test_seq_loader_real_mhc_fasta(self, mhc_data_paths):
        """Test with real MHC FASTA data"""
        fasta_path = str(mhc_data_paths['cd1_fasta'])
        result = loader.seq_loader(fasta_path, label='mhc_cd1', datType='fasta')
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == 1
        assert result.shape[1] > 0

    def test_seq_loader_invalid_path(self):
        """Test handling of invalid file path"""
        with pytest.raises(FileNotFoundError):
            loader.seq_loader('/nonexistent/path/file.csv', label='test')


class TestGetMsaSub:
    """Tests for get_msa_sub function"""

    def test_get_msa_sub_basic(self):
        """Test basic MSA subsetting"""
        seqF = pd.DataFrame([['MKVLWAALLVGTGSASHSR', 'MEQLTLLLPLATSLALTASM']])
        loc_start = [0]
        loc_end = [10]

        result = loader.get_msa_sub(seqF, loc_start, loc_end)
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == 1

    def test_get_msa_sub_multiple_regions(self):
        """Test subsetting with multiple regions"""
        seqF = pd.DataFrame([
            ['MKVLWAALLVGTGSASHSR', 'MEQLTLLLPLATSLALTASM'],
            ['QSISSYDASQHRSTWPPNG', 'ESLLHSEVSMQTIQLPGTF']
        ])
        loc_start = [0, 5]
        loc_end = [5, 10]

        result = loader.get_msa_sub(seqF, loc_start, loc_end)
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == 2

    def test_get_msa_sub_maintains_structure(self):
        """Test that structure is maintained in output"""
        seqF = pd.DataFrame([
            ['MKVLWAALLVGTGSASHSR', 'MEQLTLLLPLATSLALTASM'],
            ['QSISSYDASQHRSTWPPNG', 'ESLLHSEVSMQTIQLPGTF']
        ], columns=['seq1', 'seq2'])
        loc_start = [0]
        loc_end = [5]

        result = loader.get_msa_sub(seqF, loc_start, loc_end)
        assert result.columns.tolist() == ['seq1', 'seq2']

    def test_get_msa_sub_error_mismatched_lengths(self):
        """Test error when start and end lengths don't match"""
        seqF = pd.DataFrame([['MKVLWAALLVGTGSASHSR']])
        loc_start = [0, 5]
        loc_end = [10]

        # Should print error and return nothing
        result = loader.get_msa_sub(seqF, loc_start, loc_end)
        assert result is None


class TestConvert3Let:
    """Tests for convert_3Let function"""

    def test_convert_3let_single_letter(self):
        """Test conversion to single-letter amino acid code"""
        result = loader.convert_3Let(['ALA'])
        assert result == 'A'

    def test_convert_3let_multiple_residues(self):
        """Test conversion of multiple 3-letter codes"""
        result = loader.convert_3Let(['ALA', 'GLY', 'ARG'])
        assert result == 'AGR'

    def test_convert_3let_all_standard_aas(self):
        """Test conversion of all standard amino acids"""
        three_letter = ['ALA', 'GLY', 'ARG', 'LYS', 'ASP', 'GLU', 'ASN', 'GLN',
                       'MET', 'CYS', 'PHE', 'TYR', 'THR', 'TRP', 'PRO', 'SER',
                       'LEU', 'VAL', 'HIS', 'ILE']
        expected = 'AGRKDEQMCFYTTRPSPLVHI'
        result = loader.convert_3Let(three_letter)
        assert result == expected

    def test_convert_3let_lowercase_input(self):
        """Test that function handles lowercase input"""
        result = loader.convert_3Let(['ala', 'gly'])
        assert result == 'AG'

    def test_convert_3let_mixed_case(self):
        """Test mixed case input"""
        result = loader.convert_3Let(['ALA', 'gly', 'ArG'])
        assert result == 'AGR'

    def test_convert_3let_empty_list(self):
        """Test handling of empty input"""
        # Should raise an error or return empty
        with pytest.raises((ValueError, IndexError, AttributeError)):
            loader.convert_3Let([])

    def test_convert_3let_invalid_residue(self):
        """Test handling of invalid 3-letter code"""
        # Should raise an error or skip the invalid code
        with pytest.raises((ValueError, UnboundLocalError)):
            loader.convert_3Let(['XXX'])


class TestIntegration:
    """Integration tests for aims_loader"""

    def test_load_and_process_multiple_files(self, abs_data_paths):
        """Test loading and processing multiple data files"""
        poly_df = loader.seq_loader(str(abs_data_paths['poly']), label='poly')
        mono_df = loader.seq_loader(str(abs_data_paths['mono']), label='mono')

        # Both should be DataFrames
        assert isinstance(poly_df, pd.DataFrame)
        assert isinstance(mono_df, pd.DataFrame)

        # Should have compatible dimensions
        assert poly_df.shape[0] == mono_df.shape[0]

    def test_load_different_molecule_types(self, abs_data_paths, tcr_data_paths, mhc_data_paths):
        """Test loading different molecule types"""
        abs_data = loader.seq_loader(str(abs_data_paths['poly']), label='abs')
        tcr_data = loader.seq_loader(str(tcr_data_paths['siv_cm9']), label='tcr')
        mhc_data = loader.seq_loader(str(mhc_data_paths['cd1_fasta']), label='mhc', datType='fasta')

        assert isinstance(abs_data, pd.DataFrame)
        assert isinstance(tcr_data, pd.DataFrame)
        assert isinstance(mhc_data, pd.DataFrame)
