import pytest
import src.annotate_phylogenetic_tree as utga
import logging
import os

# Load in reference vcf data
@pytest.fixture
def vcf_file():
    " Return test VCF file"
    vcf_fn = '../test_data/HA_all_sequences.filt.sam.vcf'
    return vcf_fn

def test_build_tree(vcf_file):
    ' Test Usher generates phylogentic tree '
    utga.build_tree(vcf_file)
    expected_output_file = 'HA_all_sequences.filt.sam.pb'
    assert os.path.isfile(expected_output_file), 'expected Usher phylogenetic tree file not found.'

def test_annotate_tree():
    ' Test annotate tree file is generated '

    pass
