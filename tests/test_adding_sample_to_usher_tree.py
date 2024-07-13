import pytest
import src.adding_sample_to_usher_tree as astut
import logging
import os
import pandas
import logging
import src.annotate_phylogenetic_tree as apt


#################################################
#### Fixtures
#################################################

@pytest.fixture
def project_dir():
    " Return root directory of project "
    script_fp = os.path.dirname(__file__)
    project_fp = '/'.join(script_fp.split('/')[:-2])
    return project_fp

@pytest.fixture
def output_folder(project_dir):
    " Return output folder"
    pb_fn = os.path.join(str(project_dir), 'test_out')
    return pb_fn

# Load in reference pb data
@pytest.fixture
def pb_file(project_dir):
    " Return test PB tree file"
    pb_fn = os.path.join(str(project_dir), 'test_data/HA_all_sequences.filt.pb')
    return pb_fn

@pytest.fixture
def query_seqs(project_dir):
    " Return query multi-fasta fasta file "
    fasta_fn = os.path.join(str(project_dir), 'test_data/sequences_to_add_test_set.fasta')
    return fasta_fn

@pytest.fixture
def inappropriate_query_seqs():
    " Return a poor quality multi-fasta file (COVID not AI genome) "
    fasta_fn = os.path.joni(str(project_dir), 'test_data/inappropriate_sequences_to_add.fasta')
    return fasta_fn

@pytest.fixture
def distant_query_seqs():
    ' Return a phylogenetically distant multi-fasta file (H1N2 not H5N1) '
    fasta_fn = os.path.join(str(project_dir), 'test_data/H1N2_sequences_to_add.fasta')
    return fasta_fn

# Load in reference fasta data
@pytest.fixture
def reference_file(project_dir):
    "Returns genkbank file in usher genotyping package test_data folder"
    ref_file = os.path.join(str(project_dir), 'test_data/A.goose.Guangdong.1.1996_HA.fasta')
    return ref_file

# Test vcf file
@pytest.fixture
def vcf_file(project_dir):
    ' Returns vcf file containing samples to add to pb tree file '
    vcf_file = os.path.join(str(project_dir), 'test_data/sequences_to_add.vcf')
    return vcf_file

# Test PB tree file with additional samples
@pytest.fixture
def expanded_pb_file(project_dir):
    ' Returns a PB tree file that has had samples added to it '
    expanded_pb_fp = os.path.join(str(project_dir), 'test_data/HA_all_sequences.filt.updated.pb')
    return expanded_pb_fp

@pytest.fixture
def perfect_meta_file(project_dir):
    " Return tsv file with sample name and clade, no header."
    meta_fn = os.path.join(project_dir, 'test_data/20230723_test_genotypes.perfect.csv')
    return meta_fn

@pytest.fixture
def annotated_updated_pb_file(project_dir):
    ' Return tree with additional samples and updated clade annotations '
    annotated_updated_fn = os.path.join(project_dir, 'test_data/HA_all_sequences.filt.updated.annotated.pb')
    return annotated_updated_fn


#################################################
#### Unit Tests
#################################################

def test_generate_query_vcf(reference_file, output_folder, query_seqs):
    ' Test vcf file is generated for query sequences '

    #TODO: not all sequences get mapped from fasta file?
    vcf = astut.generate_query_vcf(reference_file, output_folder, 2, query_seqs)
    assert vcf

def test_add_vcf_to_usher_tree(pb_file, vcf_file, output_folder, caplog):
    ' Test that additional sample(s) get added to Protobuf Tree using a known vcf file'
    caplog.set_level(logging.INFO)

    # Create folder to save summary data out to
    add_samples_output_folder = str('.'.join(os.path.basename(pb_file).split('.')[:-1])+str('.additional_samples'))
    add_samples_output_folder = os.path.join(output_folder, add_samples_output_folder)
    if not os.path.isdir(add_samples_output_folder):
        os.mkdir(add_samples_output_folder)

    astut.add_vcf_to_tree(pb_file, vcf_file, 2, output_folder)

    # Assert expected samples are generated
    for outfile in ['uncondensed-final-tree.nh', 'placement_stats.tsv', 'mutation-paths.txt']:
        file_path = os.path.join(add_samples_output_folder, outfile)
        assert os.path.isfile(file_path)

    # Create name for new protobuf tree
    updated_pb_tree_fn = str('.'.join(os.path.basename(pb_file).split('.')[:-1])+('.updated.pb'))
    pb_samples_fp = os.path.join(output_folder, updated_pb_tree_fn)
    assert pb_samples_fp


def test_add_test_sample_set_to_usher_tree(reference_file, output_folder, query_seqs, pb_file, caplog):
    ' Test adding a mixture of good and bad quality samples to Usher tree '

    # Create folder to save summary data out to
    add_samples_output_folder = str('.'.join(os.path.basename(pb_file).split('.')[:-1])+str('.additional_samples'))
    add_samples_output_folder = os.path.join(output_folder, add_samples_output_folder)
    if not os.path.isdir(add_samples_output_folder):
        os.mkdir(add_samples_output_folder)
    
    vcf = astut.generate_query_vcf(reference_file, output_folder, 2, query_seqs)
    astut.add_vcf_to_tree(pb_file, vcf, 2, output_folder)

    # Assert expected samples are generated
    for outfile in ['uncondensed-final-tree.nh', 'placement_stats.tsv', 'mutation-paths.txt']:
        file_path = os.path.join(add_samples_output_folder, outfile)
        assert os.path.isfile(file_path)

    # Create name for new protobuf tree
    updated_pb_tree_fn = str('.'.join(os.path.basename(pb_file).split('.')[:-1])+('.updated.pb'))
    pb_samples_fp = os.path.join(output_folder, updated_pb_tree_fn)
    assert pb_samples_fp

    # TODO: assert poor quality samples?

def test_annotate_tree(expanded_pb_file, perfect_meta_file, output_folder):
    ' test annotated tree file is generated with log files '
    apt.annotate_tree(expanded_pb_file, perfect_meta_file, output_folder, 2)
    # Test output file has been generated
    pb_annotated_fn = str('.'.join(os.path.basename(expanded_pb_file).split('.')[:-1])+str('.annotated.pb'))
    pb_annotated_fp = os.path.join(output_folder, pb_annotated_fn)
    assert os.path.isfile(pb_annotated_fp) # test tree file generated

    # Log file
    log_fn = str('.'.join(os.path.basename(expanded_pb_file).split('.')[:-1])+str('.annotated.log'))
    log_fp = os.path.join(output_folder, log_fn)
    assert os.path.isfile(log_fp) # check log file has been generated

    # TODO: probably need to check if clade output file is generated


def test_write_summary_output_files(annotated_updated_pb_file, output_folder):
    ' Test if summary outputs write as expected '
    astut.write_summary_output_files(annotated_updated_pb_file, output_folder)
    
    summary_output_folder = str('.'.join(os.path.basename(annotated_updated_pb_file).split('.')[:-1])+str('.summary_stats'))
    summary_output_folder = os.path.join(output_folder, summary_output_folder)

    for out_file in ['samples.tsv','mutations.tsv','clades.tsv','aberrant.tsv']:
        file_path = os.path.join(summary_output_folder, str('samples.tsv'))
        assert os.path.isfile(file_path)