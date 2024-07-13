from typing import Optional
from typing import Sequence
import argparse
from usher_genotyping.src.logger import log
import os
import subprocess
import pandas
from usher_genotyping.src.multi_sequence_alignment import check_dependency_installed
import sys

"""
Annotate Protobuff Phylogenetic Tree with Genotypes and SNP mutations
"""

def check_output_file(file_name):
    '''
    check output files exist
    :params: file name to check if present
    :return:
    '''
    if os.path.isfile(file_name):
        log("info", f"======= {file_name} generated.")
    else:
        log("critical", f"======= {file_name} has not been generated.")


def get_pb_sample_ids(pb_fn):
    '''
    Get sample IDs from PB tree file
    :params: protobuf file
    :return:    dataframe with PB file data ['sample','parsimony','parent_id']
                dataframe for sample lookup ['sample','tree (boolean)']
    '''

    # Set up directories and file name
    pb_samples_fn = str('.'.join(os.path.basename(pb_fn).split('.')[:-1])+('.pb_samples.txt'))
    working_dir = os.path.dirname(pb_fn)
    pb_samples_fp = os.path.join(working_dir, pb_samples_fn)
    # get txt file containing all sample names
    pb_sample_extract_process = subprocess.Popen(['matUtils', 'summary', '-i', pb_fn, 
                                '--samples', pb_samples_fn, '-d', working_dir])    
    pb_sample_extract_process.wait()

    # check file is generated
    check_output_file(pb_samples_fp)

    # Import file and check if file has contents
    pb_samples_df = pandas.read_csv(pb_samples_fp, sep='\t')

    # Create look-up table
    pb_lookup_df = pb_samples_df.loc[:, ['sample']]
    pb_lookup_df['Tree'] = True

    if not pb_samples_df.empty:
        return pb_samples_df, pb_lookup_df
    else:
        log("critical", f"======= {pb_samples_fp} has no samples present. Check PB tree file was generated as expected")
    

def parse_metadata_file(meta_fn):
    '''
    Read in metadata tsv
    :params: metadata tsv containing genotype and sample id, no header
    :return:    pandas df ['genotype','sample']
                lookup df ['sample', 'genotype (bool)']
    '''
    meta_df =pandas.read_csv(meta_fn, sep='\t', header=None)
    meta_df.columns = ['genotype','sample']

    # Create lookup table for metadata
    meta_lookup_df = meta_df.loc[:, ['sample']]
    meta_lookup_df['Genotype'] = True

    if not meta_df.empty:
        return meta_df, meta_lookup_df
    else:
        log("critical", f"======= {meta_fn} failed to be parsed successfully. Check formatting and context")


def sample_percentage(dataframe, total):
    return (100 / total) * dataframe.index.size

def cross_reference_lookup_dfs(pb_lookup_df, meta_lookup_df):
    '''
    Cross reference the sample ids in pb tree and metadata, report matches discrepencies
    :params:    pb samples data pandas dataframe ['sample', 'parsimony', 'parent_id'], 
                genotype metadata pandas dataframe ['genotype', 'sample']
    :return: tsv report on matches, none-matches, fail if no samples are matched
    '''
    merged_df = pb_lookup_df.merge(meta_lookup_df, on='sample', how='outer')
    merged_df = merged_df.fillna(False)

    # Dfs for samples with both Tree and Genotype, only tree, only genotype
    matched_df = merged_df.loc[(merged_df['Tree']  == True) & (merged_df['Genotype']  == True)]
    tree_only_df = merged_df.loc[(merged_df['Tree']  == True) & (merged_df['Genotype']  == False)]
    genotype_only_df = merged_df.loc[(merged_df['Tree']  == False) & (merged_df['Genotype']  == True)]

    if matched_df.empty == True: # CRITICAL error if no matches are found
        log('critical', f"====== No matches identified between tree and genoptyping sample ids. Please check files.")
        sys.exit() #TODO: change this to do automatic annotation of clades
    else:
        # Counts
        matched_count = matched_df.index.size
        tree_only_count = tree_only_df.index.size
        genotype_only_count = genotype_only_df.index.size

        # Percentage
        matched_perc = sample_percentage(matched_df, merged_df.index.size)
        tree_only_perc = sample_percentage(tree_only_df, merged_df.index.size)
        genotype_only_perc = sample_percentage(genotype_only_df, merged_df.index.size)

        # Generate summary output
        summary_dict = {'Summary': ['Total Samples', 'Tree and Genotyped', 'Tree only', 'Genotype Only'], 
                    'Count':  [merged_df.index.size, matched_count, tree_only_count, genotype_only_count], 
                    'Percentage':[float(100), matched_perc, tree_only_perc, genotype_only_perc] }
        summary_df = pandas.DataFrame(summary_dict)


        return summary_df

def reviewing_matched_results(summary_df, output):
    '''
    Review summary dataframe and report appropraitely
    :params: pandas dataframe with useful results
    :return: appropriate logs
    '''
    # Logging results
    if summary_df['Percentage'][1] == float(100):
        log('info', f"====== All samples in the tree have been matched in genotyping tsv")
    elif summary_df['Percentage'][1] <= float(50):
        log('warning', f"====== Less than 50% of samples have been matched between tree and gentoyping file. ({summary_df['Percentage'][1]}%)")
    else:
        log('info', f"====== Total percentage of matched samples: {summary_df['Percentage'][1]}%")

    output_fp = os.path.join(output, 'tree_genotyping_sample_matching_summary.tsv')
    summary_df.to_csv(output_fp, sep='\t', index=False)


def annotate_tree(pb_tree_fn, metadata, output, threads):
    '''
    Annotate the phylogenetic tree generated by usher with Clade/Genotyping information
    :param: usher phylogenetic tree, metadata tsv
    :return: annotate phylogenetic tree and log file
    '''
    # craete annotations sub-directory
    annotations_dir = os.path.join(output, str('annotations'))
    if not os.path.isdir(annotations_dir):
        os.mkdir(annotations_dir)

    # Create filename for annotated pb file
    pb_annotated_fn = str('.'.join(os.path.basename(pb_tree_fn).split('.')[:-1])+str('.annotated.pb'))
    pb_annotated_fp = os.path.join(output, pb_annotated_fn)

    # Details file
    details_fp = os.path.join(annotations_dir, 
                              str('.'.join(os.path.basename(pb_tree_fn).split('.')[:-1])
                                  +str('.annotation_details.tsv')))
    
    
    # Mutations file
    mutations_fp = os.path.join(annotations_dir, 
                              str('.'.join(os.path.basename(pb_tree_fn).split('.')[:-1])
                                  +str('.mutation_details.tsv')))


    # Annotate Phylogenetic Tree with clade information
    annotation_process = subprocess.Popen(['matUtils', 'annotate', 
                                           '-i', pb_tree_fn, # input mutation anntotated tree
                                           '-c', metadata, # sample clade names 
                                           '-o', pb_annotated_fp, # Path to processed mutation annotated tree
                                           '-D', details_fp, # details of nodes considered for each clade root
                                           '-u', mutations_fp, # listing each clade and the mutations found in at least 
                                           '-T', str(threads)] # threads
                                        ) 
    annotation_process.wait()

    for file_name in [pb_annotated_fp, pb_annotated_fp, mutations_fp]:
        check_output_file(file_name)

    return annotations_dir

def main(argv: Optional[Sequence[str]] = None) -> int:
    # read in arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--pb_file', '-p', required=True, help='PB Tree input file.')
    parser.add_argument('--metadata', '-m', required=True, help='Clade metadata tsv file.')
    parser.add_argument('--output', '-o', required=False, help='Output directory to save files to, if non-specified uses PB Tree file dir.')
    parser.add_argument('--threads', '-T', required=False, default=2, help='VCF input file.')

    args = parser.parse_args(argv)

    # handle output dir
    if not args.output:
        output = os.path.dirname(args.pb_file)
    else:
        output = args.output

    # Check matUtils
    check_dependency_installed('matUtils')

    # Get sample ids from PB file
    pb_samples_df, pb_lookup_df = get_pb_sample_ids(args.pb_file)

    # Parse genotype metadata
    meta_df, meta_lookup_df = parse_metadata_file(args.metadata)

    # Cross reference samples in tree with genotyped samples
    matched_summary_df = cross_reference_lookup_dfs(pb_lookup_df, meta_lookup_df)

    # Review cross referening results:
    reviewing_matched_results(matched_summary_df, output)

    # Annotate Tree with clade annotations
    annotations_dir = annotate_tree(args.pb_file, args.metadata, output, args.threads)


if __name__ == "__main__":
  exit(main())