import pytest
import os
from Bio import SeqIO
import src.parse_genbank_ref as pgr


# Load in Genbank data
@pytest.fixture
def genkbank_file():
    "Returns genkbank file in usher genotyping package test_data folder"
    gb_file = "../test_data/A.goose.Guangdong.1.1996_HA.gb"
    return gb_file

# Load in Genbank data
@pytest.fixture
def fasta_file():
    "Returns genbank file in usher genotyping package test_data folder"
    fasta_file = "../test_data/HA_all_sequences.filt.fasta"
    return fasta_file

@pytest.fixture
def bad_genbank_file():
    " Returns a genbank file missing essential parameters "
    bad_gb_file = "../test_data/bad_genbank_data.gb"
    return bad_gb_file

@pytest.fixture
def genbank_data(genkbank_file):
    '''
    Take genbank filename and generate genkbank object
    :params: file name as string
    :return: BioSeq Genbank Object
    '''
    with open(genkbank_file, "r") as file:
        gb_data = SeqIO.read(file, 'gb')
        return gb_data
  

# Unit Testing - genbank object attributes
def test_wrong_file_type(fasta_file):
    """ Test exception is raised when the wrong file type is provided """
    with pytest.raises(ValueError):
        pgr.parse_genbank_file(fasta_file)

def test_bad_genbank_file(bad_genbank_file):
     """ Test exception is raised when a corrupt file is provided """
     with pytest.raises(ValueError):
         pgr.parse_genbank_file(bad_genbank_file)


def test_check_correct_genbank_attributes(genbank_data, capsys):
    '''
    Assert expected attributes are present in Genbank data
    :params: genkbank object
    :return: asserts output print statement (does not fail)
    '''
    pgr.check_genbank_attributes(genbank_data)
    captured = capsys.readouterr()
    assert captured.out == 'Genbank File appears as expected\n'    

def test_gbmunge_installed(genkbank_file, capsys):
    '''
    Test behaviour for when gbmunge is installed
    '''
    pgr.run_gbmunge(genkbank_file)
    captured = capsys.readouterr()
    # Check gbmunge runs without errors
    assert captured.err == '', "Error running gbmunge"
    # Check fasta files are as expected
    assert captured.out == "Fasta and metadata generated from Genbank file.\n", "outputf ilse are not as expected."   