import pytest
import src.annotate_phylogenetic_tree as apt
import logging
import os
import pandas


#################################################
#### Fixtures
#################################################

@pytest.fixture
def project_dir():
    " Return root directory of project "
    script_fp = os.path.dirname(__file__)
    project_fp = '/'.join(script_fp.split('/')[:-2])
    return project_fp

# Load in reference pb data
@pytest.fixture
def pb_file(project_dir):
    " Return test VCF file"
    pb_fn = os.path.join(str(project_dir), 'test_data/HA_all_sequences.filt.pb')
    return pb_fn

# Output folder for tests (probably should be a temporary construct)
@pytest.fixture
def output_folder(project_dir):
    " Return output folder"
    pb_fn = os.path.join(str(project_dir), 'test_out')
    return pb_fn

@pytest.fixture
def perfect_meta_file(project_dir):
    " Return tsv file with sample name and clade, no header."
    meta_fn = os.path.join(project_dir, 'test_data/20230723_test_genotypes.perfect.csv')
    return meta_fn

@pytest.fixture
def good_meta_file(project_dir):
    " Return tsv file with sample name and clade, no header."
    meta_fn = os.path.join(project_dir, 'test_data/20230723_test_genotypes.good.csv')
    return meta_fn

@pytest.fixture
def poor_meta_file(project_dir):
    " Return tsv file with sample name and clade, no header."
    meta_fn = os.path.join(project_dir, 'test_data/20230723_test_genotypes.poor.csv')
    return meta_fn

@pytest.fixture
def fail_meta_file(project_dir):
    " Return tsv file with sample name and clade, no header."
    meta_fn = os.path.join(project_dir, 'test_data/20230723_test_genotypes.fail.csv')
    return meta_fn

#################################################
#### Unit Tests
#################################################


def test_get_pb_sample_ids(pb_file):
    ' Test that PB Tree dataframes have been generated and populated'
    pb_samples_df, pb_lookup_df = apt.get_pb_sample_ids(pb_file)
    assert pb_samples_df.empty == False
    assert pb_lookup_df.empty == False


def test_parse_metadata_file(perfect_meta_file):
    ' Check metdata dataframes have been populated '
    meta_df, meta_lookup_df = apt.parse_metadata_file(perfect_meta_file)
    assert meta_df.empty == False
    assert  meta_lookup_df.empty == False


# This seems verbose?
def test_cross_reference_lookup_perfect(pb_file, perfect_meta_file):
    ' Test that if all samples are linked there is appropriate output '
    pb_samples_df, pb_lookup_df = apt.get_pb_sample_ids(pb_file)
    meta_df, meta_lookup_df = apt.parse_metadata_file(perfect_meta_file)

    summary_df = apt.cross_reference_lookup_dfs(pb_lookup_df, meta_lookup_df)
    assert summary_df.empty == False
    assert summary_df['Percentage'][0] == 100 # assert that there is an exact match

def test_reviewing_match_results(output_folder, caplog):
    ' Check the summary output file is correct '
    caplog.set_level(logging.INFO)
    
    # Perfect results
    summary_dict = {'Summary': ['Total Samples', 'Tree and Genotyped', 'Tree only', 'Genotype Only'], 
                    'Count':  [100, 100, 0, 0], 
                    'Percentage':[float(100), float(100), float(0), float(0)] }
    perfect_summary_df = pandas.DataFrame(summary_dict)
    apt.reviewing_matched_results(perfect_summary_df, output_folder)

    # Test output file is generated
    output_fp = os.path.join(output_folder, 'tree_genotyping_sample_matching_summary.tsv')
    assert os.path.isfile(output_fp)

    # Good results
    summary_dict = {'Summary': ['Total Samples', 'Tree and Genotyped', 'Tree only', 'Genotype Only'], 
                    'Count':  [100, 51, 49, 0], 
                    'Percentage':[float(100), float(51), float(49), float(0)] }
    good_summary_df = pandas.DataFrame(summary_dict)
    apt.reviewing_matched_results(good_summary_df, output_folder)
    
    # Poor Results
    summary_dict = {'Summary': ['Total Samples', 'Tree and Genotyped', 'Tree only', 'Genotype Only'], 
                    'Count':  [100, 49, 51, 0], 
                    'Percentage':[float(100), float(49), float(51), float(0)] }
    poor_summary_df = pandas.DataFrame(summary_dict)
    apt.reviewing_matched_results(poor_summary_df, output_folder)
    
    # Test logging
    assert 'All samples in the tree have been matched in genotyping tsv' in caplog.records[0].msg
    assert 'Total percentage of matched samples: 51.0%' in caplog.records[1].msg
    assert 'Less than 50% of samples have been matched between tree and gentoyping file. (49.0%)' in caplog.records[2].msg
    

def test_cross_reference_lookup_fail(pb_file, fail_meta_file, caplog):
    ' Test that if all samples are linked there is appropriate output '
    caplog.set_level(logging.CRITICAL)
    pb_samples_df, pb_lookup_df = apt.get_pb_sample_ids(pb_file)
    meta_df, meta_lookup_df = apt.parse_metadata_file(fail_meta_file)

    summary_df = apt.cross_reference_lookup_dfs(pb_lookup_df, meta_lookup_df)
    assert 'CRITICAL' in caplog.text # make sure a critical log is raised


def test_annotate_tree(pb_file, perfect_meta_file, output_folder):
    ' test annotated tree file is generated with log files '
    apt.annotate_tree(pb_file, perfect_meta_file, output_folder, 2)
    # /home/phe.gov.uk/nicholas.ellaby/Documents/git_repos/dev_code/usher_genotyping/test_out
    # Test output file has been generated
    pb_annotated_fn = str('.'.join(os.path.basename(pb_file).split('.')[:-1])+('.annotated.pb'))
    pb_annotated_fp = os.path.join(output_folder, pb_annotated_fn)
    assert os.path.isfile(pb_annotated_fp) # test tree file generated

    # Log file
    log_fn = str('.'.join(os.path.basename(pb_file).split('.')[:-1])+('.annotated.log'))
    log_fp = os.path.join(output_folder, log_fn)
    assert os.path.isfile(log_fp) # check log file has been generated