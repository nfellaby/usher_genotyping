import pytest
import src.multi_sequence_alignment as msa
import logging
import os

# Load in reference fasta data
@pytest.fixture
def reference_file():
    "Returns genkbank file in usher genotyping package test_data folder"
    ref_file = "../test_data/A.goose.Guangdong.1.1996_HA.fasta"
    return ref_file

# Load in multi-fasta sequence data
@pytest.fixture
def fasta_file():
    "Returns genbank file in usher genotyping package test_data folder"
    fasta_file = "../test_data/HA_all_sequences.filt.fasta"
    return fasta_file

@pytest.fixture
def reference_index():
    "Returns reference index file"
    index_file = '../test_data/A.goose.Guangdong.1.1996_HA.mni'
    return index_file

@pytest.fixture
def reference_sam():
    " returns SAM file "
    sam_file = '../test_data/HA_all_sequences.filt.sam'
    return sam_file

@pytest.fixture
def output_folder():
    " returns output folder "
    output_folder = './'
    return output_folder

@pytest.fixture
def aln_file():
    " returns output folder "
    aln_fn = '../test_data/HA_all_sequences.filt.sam.aln'
    return aln_fn

# Unit tests
def test_check_depency_installed(capsys):
    " Test whether a dependecy is installed, bash should always be installed"
    msa.check_dependency_installed('bash')
    captured = capsys.readouterr()
    assert captured.out == f"bash found in 'which' path. Looks to be installed.\n"

def test_index(reference_file):
    # Warning when overwriting existing files
    # assert caplog.record_tuples == [("root", logging.WARNING, "root:logger.py:50 ======= Index file for reference already exists. This files will be over-written.")]
    reference_index = msa.generate_reference_index(reference_file)
    assert reference_index == 'A.goose.Guangdong.1.1996_HA.mni'
    assert os.path.exists('A.goose.Guangdong.1.1996_HA.mni')

def test_generate_alignment(fasta_file, reference_index):
    " Check that alignment is generated "
    sam_filename = msa.generate_alignment(fasta_file, reference_index)
    assert sam_filename == 'HA_all_sequences.filt.sam'
    assert os.path.exists('HA_all_sequences.filt.sam')

def test_sam_to_aln(reference_sam, reference_file, output_folder):
    " Check expected vcf file is generated "
    msa.sam_to_aln(reference_sam, reference_file, output_folder)
    assert os.path.isfile(os.path.join(output_folder, 'HA_all_sequences.filt.sam.aln'))

def test_faToVCF(aln_file):
    ''' Check VCF file generated '''
    msa.faToVCF(aln_file)
    assert os.path.isfile('HA_all_sequences.filt.sam.vcf')
