import argparse
import sys
import os
import logging
import glob
from usher_genotyping.src.logger import log, setup_logger
import usher_genotyping.src.parse_genbank_ref as pgr
import usher_genotyping.src.parse_gff_ref as pgff
import usher_genotyping.src.multi_sequence_alignment as msa
import usher_genotyping.src.build_phylogenetic_tree as bpt
import usher_genotyping.src.annotate_phylogenetic_tree as apt
import usher_genotyping.src.summarise_phylogenetic_annotations as spa
import usher_genotyping.src.adding_sample_to_usher_tree as astut
import subprocess

try:
    from .version import version as __version__
    from .version import version_tuple
except ImportError:
    __version__ = "0.01"
    version_tuple = (0, 0, "unknown version")



# Main script to run each sub-script
"""
The main function of this file is to handle the phylogenetic processing of sequences using usher and matutuils
"""

def args_parser():
    """
    Command line argument parser
    """
    parser = argparse.ArgumentParser(description="Process for Usher Genotyping.")
    parser.add_argument('--sequences', '-s', required=True, help='Contextual FASTA sequences used to make phylogenetic tree.')
    ref_group = parser.add_mutually_exclusive_group(required=True)
    ref_group.add_argument('--genbank', '-r',  help='Reference Genbank used to make phylogenetic tree.')
    ref_group.add_argument('--fasta', '-f',  help='Reference FASTA used to make phylogenetic tree.')
    parser.add_argument('--gff',  help='Reference GFF used to make phylogenetic tree. This is required to infer mutations.')
    parser.add_argument('--genotypes', '-g', required=True, help='Assigned Genotypes (tsv). First column should be genotype, second column sample id. No header.')
    parser.add_argument('--output', '-o', required=True, help='Output directory to save vcf')
    parser.add_argument('--threads', '-T', required=False, default=2, help='VCF input file.')
    parser.add_argument('--overwrite', '-w', required=False, default=str('no'), 
                        help='If files are already present, overwrite files. Otherwise, skip step.',
                        choices=['yes','no'])
    args=parser.parse_args()
    return args

def parse_genbank(reference:str, output:str):
    '''
    Process reference genbank into fasta sequence
    :params: reference file path, output folder filepath
    :return: output fasta file path
    '''
    gb_object = pgr.parse_genbank_file(reference)
    gb_object = pgr.check_genbank_attributes(gb_object)
    pgr.run_gbmunge(reference, output)
    # output file
    ref_name = '.'.join(os.path.basename(reference).split('.')[:-1])
    ref_fasta_name = ref_name+str('.fasta')
    ref_fasta_file = os.path.join(output, ref_fasta_name)
    return ref_fasta_file

def parse_gff(reference:str, output:str):
    '''
    Extract mutation information from gff files
    :params: gff reference file path, output folder filepath
    :return: mutation file (GTF) filepath
    '''
    #TODO: parse gff into GTF for mutation inference
    pass


def run_analysis_test(output_fp, overwrite):
    '''
    Check if output files or directories are present, and handle the 'overwrite' argument.
    :params: output file path to check (str), overwrite argument ['yes','no']
    :return: True/False
    '''
    if os.path.isfile(output_fp):
        log("info", f"======= Output file/dir already present: {output_fp}")
        if overwrite == 'yes':
            log("warning", f"======= Overwrite argument set to: {overwrite}. Re-creating file/dir...")
            return True
        else:
            log("info", f"======= Overwrite argument set to: {overwrite}. Using original file/dir.")
            return False
    else:
        log("info", f"======= File/dir not present: {output_fp}. Creating reference index...")
        return True

def alignment(sequences:str, reference:str, output:str, threads:int, overwrite:str):
    '''
    Use minimap to map sequences to reference to generated alignment
    :params: sequences filepath, reference filepath
    :return: alignment filepath
    '''
    # Create reference file
    reference_index = str('.'.join(os.path.basename(reference).split('.')[:-1])+('.mni'))
    reference_index = os.path.join(output, reference_index)
    run_ref_index_analysis = run_analysis_test(reference_index, overwrite)
    if run_ref_index_analysis == True:
        reference_index = msa.generate_reference_index(reference, output, threads)

    # Create SAM file
    sam_filename = str('.'.join(os.path.basename(sequences).split('.')[:-1])+('.sam'))
    sam_file_path = os.path.join(os.path.dirname(reference_index), sam_filename)
    run_sam_analysis = run_analysis_test(sam_file_path, overwrite)
    if run_sam_analysis == True:
        sam_file_path = msa.generate_alignment(sequences, reference_index, threads)

    # Conver SAM to alignment file
    out_msa_fn =   str('.'.join(os.path.basename(sam_file_path).split('.')[:-1])+('.aln'))
    aln_filename =  os.path.join(os.path.dirname(sam_file_path), out_msa_fn)
    run_aln_analysis = run_analysis_test(aln_filename, overwrite)
    if run_aln_analysis == True:
        aln_filename = msa.sam_to_aln(sam_file_path, reference)
     
    # Convert Alignment file to VCF
    vcf_fh =  str('.'.join(os.path.basename(aln_filename).split('.')[:-1])+('.vcf'))
    vcf_file_path = os.path.join(os.path.dirname(aln_filename), vcf_fh)
    run_vcf_analysis = run_analysis_test(vcf_file_path, overwrite)
    if run_vcf_analysis:
        msa.faToVCF(aln_filename)
    return vcf_file_path

def build_phylogenetic_tree(vcf_file:str, threads:int, output:str, overwrite:str):
    '''
    Use Usher to generate protobuf phylogenetic tree
    :params: vcf filepath
    :return: protobuf filepath
    '''
    # create pb_file_path for checking if file is already present and skipping if present
    pb_fn = '.'.join(os.path.basename(vcf_file).split('.vcf')[:-1])+str('.pb')
    working_dir = os.path.join(os.path.dirname(vcf_file), 'usher_out')

    if not os.path.isdir(working_dir):
        os.mkdir(working_dir)
    pb_file_path = os.path.join(working_dir, pb_fn)

    # Check Usher
    bpt.check_dependency_installed('usher')

    # Build tree from VCF file
    run_pb_analysis = run_analysis_test(pb_file_path, overwrite)
    if run_pb_analysis == True:
        sam_file_path = bpt.build_tree(vcf_file, threads)

    return pb_file_path

def annotate_tree(pb_file, metadata, output, threads, overwrite):
    '''
    Use matUTils to annotate protobuf tree file
    :params: protobuf tree filepath
    :return: annotated protobuf tree filepath
    '''
    # Check matUtils
    apt.check_dependency_installed('matUtils')
    # Get sample ids from PB file
    pb_samples_df, pb_lookup_df = apt.get_pb_sample_ids(pb_file)
    # Parse genotype metadata
    meta_df, meta_lookup_df = apt.parse_metadata_file(metadata)
    # Cross reference samples in tree with genotyped samples
    matched_summary_df = apt.cross_reference_lookup_dfs(pb_lookup_df, meta_lookup_df)
    # Review cross referening results:
    apt.reviewing_matched_results(matched_summary_df, output)

    # Annotate Tree with clade annotations
    annotations_dir = os.path.join(output, str('annotations'))
    if not os.path.isdir(annotations_dir):
        os.mkdir(annotations_dir)

    pb_annotated_fn = str('.'.join(os.path.basename(pb_file).split('.')[:-1])+str('.annotated.pb'))
    pb_annotated_fp = os.path.join(output, pb_annotated_fn)
    
    run_annotation_analysis = run_analysis_test(pb_annotated_fp, overwrite)
    if run_annotation_analysis == True:
        annotations_dir = apt.annotate_tree(pb_file, metadata, output, threads)

    annotated_pb_file = glob.glob(output+str('/*annotated.pb'))[0]

    return annotations_dir, annotated_pb_file

def write_summary_output_files(annotations_dir, annotated_pb_file, overwrite):
    ''''
    Write summary output files using the annotated Protobuf tree file
    :params: protobuf tree filepath, output folder
    :return:
    '''
    # Identify and pull in the annotated pb file
    

    # Expected output folders/files
    
    # Check if expected summary files are present
    expected_files = ['samples.tsv','sample_clades.tsv','mutations.tsv','clades.tsv','aberrant.tsv']
    rerun_analysis = [] # list of True/False statements for each expected file
    for expected_file in expected_files: # cycle through expected output files
        expected_file_fp = os.path.join(annotations_dir, expected_file)
        rerun_analysis.append(run_analysis_test(expected_file, overwrite))
    if True in rerun_analysis:
        summary_samples_tsv = astut.write_summary_output_files(annotated_pb_file, annotations_dir)

    # Calculate uncertaintly for samples added to tree
    equally_parsimonious_placements_tsv = os.path.join(annotations_dir, 'equally_parsimonious_placements.tsv')
    run_uncertainty_analysis = run_analysis_test(equally_parsimonious_placements_tsv, overwrite)
    if run_uncertainty_analysis == True:
        astut.calculate_sample_uncertainty(annotated_pb_file, summary_samples_tsv)

def summarise_annotations(assigned_clades, summary_stats_folder, output, overwrite):
    '''
    Summarise the inferred genotypes and phylogenetic placements
    :params: 
        assigned_clades:str filepath
        summary_stats_folder:str filepath
        output:str - output directory
    :return:
    '''
    output_filename = os.path.join(output, str('annotated_summary_stats.csv'))
    run_sum_annotations_analysis = run_analysis_test(output_filename, overwrite)
    if run_sum_annotations_analysis == True:
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
        spa.write_to_csv(summary_df_merge_annotated, output_filename)
        
def visualise_tree(output_dir):
    '''
    Use R script with tree file and metadata information to generate plots
    :params: folder containing tree file, folder containing metadata file
    :return: 
    '''
    # Find appropriate files
    metadata_fp = os.path.join(output_dir, str('annotated_summary_stats.csv')) # metadata file
    tree_fp = os.path.join(output_dir, str('../usher_out/final-tree.nh'))

    # Read in R script (since it isn't a python module)
    main_script_location = os.path.dirname(os.path.realpath(__file__))
    r_script_location = os.path.join(main_script_location, 'usher_genotyping/src/tree_visualisation.r')
    if not os.path.isfile(r_script_location):
        log('warning', 'R script file not found in expected location: {r_script_location}. Cannot generate tree figures.')
    else:
        
        visualise_tree_process = subprocess.Popen(
            ["Rscript", '--vanilla', r_script_location, tree_fp, metadata_fp, output_dir]
        )
        print(' '.join(visualise_tree_process.args))
        visualise_tree_process.wait()


def main():
    """ 
    main function
    """
    setup_logger()
    try:
        args = args_parser()

        # Make output directory
        if not os.path.isdir(args.output):
            os.mkdir(args.output)

        # Process Reference Genbank File
        if args.genbank != None:
            ref_fasta_file = parse_genbank(args.genbank, args.output)
        else:
            ref_fasta_file = args.fasta
        # mapping
        sequences_vcf_file = alignment(args.sequences, ref_fasta_file, args.output, args.threads, args.overwrite)
        # Build phylogenetic Tree
        pb_file_path = build_phylogenetic_tree(sequences_vcf_file, args.threads, args.output, args.overwrite)
        # Annotate Tree with Genotypes
        annotations_dir, annotated_pb_file = annotate_tree(pb_file_path, args.genotypes, args.output, args.threads, args.overwrite)
  
        # Write summary files
        write_summary_output_files(annotations_dir, annotated_pb_file, args.overwrite)
        # summarise analysis
        summarise_annotations(args.genotypes, annotations_dir, annotations_dir, args.overwrite)
        # Visualise Tree
        visualise_tree(annotations_dir) #TODO: This doesn't work automatically

    except:
        # log( 'critical', '')
        raise

if __name__ == "__main__":
    exit(main())
