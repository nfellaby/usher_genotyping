import argparse
import sys
import os
import logging
import glob
from usher_genotyping.src.logger import log, setup_logger
import usher_genotyping.src.adding_sample_to_usher_tree as astut
import usher_genotyping.src.annotate_phylogenetic_tree as apt
import usher_genotyping.src.summarise_phylogenetic_annotations as spa
try:
    from .version import version as __version__
    from .version import version_tuple
except ImportError:
    __version__ = "0.01"
    version_tuple = (0, 0, "unknown version")



# Main script to run each sub-script
"""
The "main" function in this file provides the command line entry point for the usher genotyping package
"""

def args_parser():
    """
    Command line argument parser
    """
    parser = argparse.ArgumentParser(description="Process for Usher Genotyping.")
    parser.add_argument('--pb_file', '-p', required=True, help='PB Tree input file.')
    parser.add_argument('--sequences', '-s', required=True, help='FASTA sequences used to add to phylogenetic tree.')
    parser.add_argument('--reference', '-r', required=True, help='Reference FASTA used to make vcf.')
    parser.add_argument('--genotypes', '-g', required=True, help='Clade metadata tsv file.')
    parser.add_argument('--output', '-o', required=True, help='Updated PB Tree')
    parser.add_argument('--threads', '-T', required=False, default=2, help='VCF input file.')
    args=parser.parse_args()
    return args

def add_samples_to_usher_tree(reference, output, threads, sequences, pb_file):
    '''
    Add Samples to existing Usher Protobuf Tree
    :params: 
    :return:
    '''
    # Generate query vcf file
    vcf = astut.generate_query_vcf(reference, output, threads, sequences)
    # Add vcf samples to existing Usher tree
    updated_tree_fp = astut.add_vcf_to_tree(pb_file, vcf, threads, output)
    return updated_tree_fp

def write_summary_output_files(annotations_dir, output_folder):
    ''''
    Write summary output files using the annotated Protobuf tree file
    :params: protobuf tree filepath, output folder
    :return:
    '''
    # Identify and pull in the annotated pb file
    annotated_pb_file = glob.glob(output_folder+str('/*annotated.pb'))[0]

    # write all possible summary output files
    summary_samples_tsv, summary_output_folder = astut.write_summary_output_files(annotated_pb_file, output_folder)

    # Calculate uncertaintly for samples added to tree
    astut.calculate_sample_uncertainty(annotated_pb_file, summary_samples_tsv, output_folder)
    
    return summary_output_folder

def summarise_annotations(assigned_clades, summary_stats_folder, output):
    '''
    Summarise the inferred genotypes and phylogenetic placements
    :params: 
        assigned_clades:str filepath
        summary_stats_folder:str filepath
        output:str - output directory
    :return:
    '''
    # Read in files
    assigned_clade_df = spa.import_original_clade_assignments(assigned_clades)
    sample_clades_df, parsimonious_placements_df, clades_summary_df = spa.import_matutils_summary_data(summary_stats_folder)

    # Review sample placements in phylogenetic tree via parsimonious information
    pp_iqr_df, pp_iqr_values = spa.parsimonious_placements_review(parsimonious_placements_df)

    # Review sample placements in phylogenetic tree via neighbourhood information
    ns_iqr_df, ns_iqr_values = spa.neighbourhood_size_review(parsimonious_placements_df)

    # Annotate quality groups to samples
    summary_df = spa.grouping_samples_phylo_assignment_scores(parsimonious_placements_df, pp_iqr_values, ns_iqr_values)

    # Add inferred and annotated sample clades to final table
    summary_df_merge = spa.add_inferred_sample_clades(summary_df, sample_clades_df, assigned_clade_df)

    # Annotate summary of inferred vs assigned clades
    summary_df_merge_annotated  = spa.summarise_clade_assignments(summary_df_merge)
    
    # # Write to file
    output_filename = os.path.join(output, str('annotated_summary_stats.csv'))
    spa.write_to_csv(summary_df_merge_annotated, output_filename)


def main():
    """ 
    main function
    
    """

    # create a global logger
    #  setup_logger()
    try:
        args = args_parser()
        # add samples to usher tree
        updated_tree_fp = add_samples_to_usher_tree(args.reference, args.output, args.threads, args.sequences, args.pb_file)
        # Annotate Tree with Genotypes
        annotations_dir = apt.annotate_tree(updated_tree_fp, args.genotypes, args.output, args.threads)
        # Write summary files
        summary_output_folder = write_summary_output_files(annotations_dir, args.output)
        # summarise analysis
        summarise_annotations(args.genotypes, summary_output_folder, args.output)

    except:
        # log( 'critical', '')
        raise

if __name__ == "__main__":
    exit(main())