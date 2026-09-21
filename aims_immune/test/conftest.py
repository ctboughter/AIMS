import pytest
import os
import numpy as np
import pandas as pd
import tempfile
from pathlib import Path

# Get the path to test data
TEST_DATA_DIR = Path(__file__).parent.parent / 'app_data' / 'test_data'

@pytest.fixture(scope="session")
def test_data_dir():
    """Return path to test data directory"""
    return TEST_DATA_DIR

@pytest.fixture
def abs_data_paths(test_data_dir):
    """Return paths to antibody test data files"""
    return {
        'poly': test_data_dir / 'abs' / 'flu_poly.csv',
        'mono': test_data_dir / 'abs' / 'flu_mono.csv',
        'poly_h3': test_data_dir / 'abs' / 'flu_poly_h3.csv',
        'mono_h3': test_data_dir / 'abs' / 'flu_mono_h3.csv',
    }

@pytest.fixture
def tcr_data_paths(test_data_dir):
    """Return paths to TCR test data files"""
    return {
        'siv_cm9': test_data_dir / 'tcrs' / 'siv_cm9.csv',
        'siv_tl8': test_data_dir / 'tcrs' / 'siv_tl8.csv',
        'vdjdb_seqs': test_data_dir / 'tcrs' / 'vdjdb_seqs_parse.csv',
        'vdjdb_meta': test_data_dir / 'tcrs' / 'vdjdb_meta_parse.csv',
    }

@pytest.fixture
def mhc_data_paths(test_data_dir):
    """Return paths to MHC test data files"""
    return {
        'cd1_fasta': test_data_dir / 'mhcs' / 'cd1.fasta',
        'classIa_fasta': test_data_dir / 'mhcs' / 'classIa.fasta',
        'fish_fasta': test_data_dir / 'mhcs' / 'fish.fasta',
        'cd1_classI_csv': test_data_dir / 'mhcs' / 'ex_cd1_classI.csv',
    }

@pytest.fixture
def peptide_data_paths(test_data_dir):
    """Return paths to peptide test data files"""
    return {
        'kidney': test_data_dir / 'peptides' / 'kidney_hla_atlas.csv',
        'pancreas': test_data_dir / 'peptides' / 'pancreas_hla_atlas.csv',
    }

@pytest.fixture
def sample_dataframe():
    """Create a sample DataFrame for testing"""
    data = {
        'l1': ['QSISSY', 'ESLLHS', 'QDIKNY'],
        'l2': ['DAS', 'EVS', 'HVS'],
        'l3': ['QHRSTWPPN', 'MQTIQLPGT', 'HQCYNLPYT'],
        'h1': ['GGTFSSRA', 'GGIMRRNG', 'GFIFGHFA'],
        'h2': ['IIPIFNTP', 'IIAIFGTP', 'ISGGGLNT'],
        'h3': ['AREMATIFGRMDV', 'VASSGYHLHRETWGY', 'ARFDSSGYNYVRGMVV']
    }
    return pd.DataFrame(data)

@pytest.fixture
def temp_fasta_file():
    """Create a temporary FASTA file for testing"""
    content = """>seq1
MKVLWAALLVGTGSASHSRQVEQAEVGQAAKSAGQSRAELHAHHHHH
>seq2
MEQLTLLLPLATSLALTASMPALAKSAGQSRAELHAHHHHH
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(content)
        temp_path = f.name

    yield temp_path

    # Cleanup
    if os.path.exists(temp_path):
        os.remove(temp_path)

@pytest.fixture
def temp_csv_file():
    """Create a temporary CSV file for testing"""
    content = """l1,l2,l3
QSISSY,DAS,QHRSTWPPN
ESLLHS,EVS,MQTIQLPGT
QDIKNY,HVS,HQCYNLPYT
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(content)
        temp_path = f.name

    yield temp_path

    # Cleanup
    if os.path.exists(temp_path):
        os.remove(temp_path)
