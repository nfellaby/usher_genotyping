import subprocess
from usher_genotyping.src.logger import log
import multi_sequence_alignment
import os
import argparse
from typing import Optional
from typing import Sequence
import annotate_phylogenetic_tree as apt 
import pandas

'''
Adding additional Samples to Usher Tree
https://usher-wiki.readthedocs.io/en/latest/UShER.html#placing-new-samples

1. Generate VCF for new samples
2. Add to existing tre
3. Summarise data

'''

def check_dependency_installed(dependency):
    '''
    check if a program required is installed
    '''
    check_process = subprocess.run(["which", dependency ], 
                        capture_output=True, text=True)
    
    if dependency in check_process.stdout:
        print(f"{dependency} found in 'which' path. Looks to be installed.")
        log("info", f"======= {dependency} installed")
    else:
        log("critical", f"======= {dependency} not found in path. Check {dependency} is installed")
        exit()

def faToVCF(out_msa_path, reference, output):
    '''
    Convert alignment file to VCF
    :params: alignment file
    :return: vcf file
    '''
    vcf_fh =  str('.'.join(os.path.basename(out_msa_path).split('.')[:-1])+('.vcf'))
    vcf_file_path = os.path.join(os.path.dirname(out_msa_path), vcf_fh)
    check_dependency_installed('faToVcf')

    # Get reference fasta header
    reference_content = open(reference, 'r')
    for line in reference_content:
        if line.startswith('>'):
            ref_header = line.split('>')[1]

    vcf_process = subprocess.Popen([
        'faToVcf', 
        out_msa_path, 
        vcf_file_path
        ])
    vcf_process.wait()

    # check output files are generated
    if os.path.isfile(vcf_file_path):
        log('info',f'======= VCF Alignment file generated: {vcf_file_path}')
        return vcf_file_path
    else:
        log('critical', f"======= VCF alignment file not generated. Please check input files.")
        exit()

def generate_query_vcf(reference, output, threads, query_sequences):
    '''
    Generate VCF file for Query sequences
    :params: reference fasta, query fasta sequecnes
    :return: vcf file path
    '''
    reference_index = multi_sequence_alignment.generate_reference_index(reference, output, threads)
    sam_file_path = multi_sequence_alignment.generate_alignment(query_sequences, reference_index, str(threads))
    aln_filename = multi_sequence_alignment.sam_to_aln(sam_file_path, reference)
    vcf = faToVCF(aln_filename, reference, output)

    # clean up directory
    # Create folder to save intermediate files to
    int_output_folder = os.path.join(output, str('intermediate_files'))
    if not os.path.isdir(int_output_folder):
        os.mkdir(int_output_folder)

    for int_file in [reference_index, sam_file_path, aln_filename]:
        new_name = os.path.join(int_output_folder, os.path.basename(int_file))
        os.replace(int_file, new_name)
    
    return vcf

def add_vcf_to_tree(pb_tree_fn, vcf_fn, threads, output_folder):
    '''
    Add samples described in vcf to PB tree
    :params: vcf file and protobuf tree
    :return: updated protobuf tree with additional samples

    usher -i global_assignments.pb -v test/new_samples.vcf -u -d output/
    '''
    # Create folder to save summary data out to
    add_samples_output_folder = str('.'.join(os.path.basename(pb_tree_fn).split('.')[:-1])+str('.additional_samples'))
    add_samples_output_folder = os.path.join(output_folder, add_samples_output_folder)
    if not os.path.isdir(add_samples_output_folder):
        os.mkdir(add_samples_output_folder)

    # Create name for new protobuf tree
    updated_pb_tree_fn = str('.'.join(os.path.basename(pb_tree_fn).split('.')[:-1])+('.updated.pb'))
    pb_samples_fp = os.path.join(output_folder, updated_pb_tree_fn)

    add_samples_prcess = subprocess.Popen(
        ['usher', 
         '-i', pb_tree_fn,
         '-v', vcf_fn,
         '-u',
         '-o', pb_samples_fp,
         '-T', str(threads),
         '-d', add_samples_output_folder])
    
    add_samples_prcess.wait()

    if not os.path.isfile(pb_samples_fp):
        log('critical', f"Updated tree file not found. Check output files and logs.")
        exit()
    else:
        return pb_samples_fp

def write_summary_output_files(pb_tree_fn, output_folder):
    '''
    Generate summary output files
    :params: Protobuf tree
    :return: summary output files
    '''

    # Run matUtils summary command
    write_summary_output_files_proc = subprocess.Popen(
        ['matUtils', 
         'summary', 
         '-i', pb_tree_fn,
         '-A', output_folder,
         '-C', str('sample_clades.tsv'),
         '-d', output_folder])
        
    write_summary_output_files_proc.wait()
    summary_samples_tsv = os.path.join(output_folder, 'samples.tsv')
    return summary_samples_tsv

def calculate_sample_uncertainty(pb_file, summary_samples_tsv):
    '''
    Calculate the uncertainty, these metricss support contact tracing and reliable identification of the origin of a newly placed sample.
    :params: protobuf tree file, samples summary file, output folder
    :return: uncertainty metrics in a tsv file
    '''

    # Extract names from summary_stats/samples.tsv
    samples_df = pandas.read_csv(summary_samples_tsv, sep='\t')
    samples_only_df = samples_df[['sample']].dropna()

    samples_only_df_fn = os.path.join(os.path.dirname(summary_samples_tsv), str('sample_list.txt'))
    samples_only_df.to_csv(samples_only_df_fn, header=False, index=False)
    
    # Record placements output file
    record_placements_fn = os.path.join(os.path.dirname(summary_samples_tsv), str('record_placements.tsv'))

    # Find EPP
    epp_fn = os.path.join(os.path.dirname(summary_samples_tsv), str('equally_parsimonious_placements.tsv'))

    # Calculate uncertainty for all samples in tree
    calculate_uncertainty_proc = subprocess.Popen(
        ['matUtils',
         'uncertainty',
         '-i', pb_file,
         '-s', samples_only_df_fn,
         '-e', epp_fn,
         '-o', record_placements_fn])
    calculate_uncertainty_proc.wait()

def main(argv: Optional[Sequence[str]] = None) -> int:
    '''
    Main function for running script directly
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--pb_file', '-p', required=True, help='PB Tree input file.')
    parser.add_argument('--sequences', '-s', required=True, help='FASTA sequences used to add to phylogenetic tree.')
    parser.add_argument('--reference', '-r', required=True, help='Reference FASTA used to make vcf.')
    parser.add_argument('--metadata', '-m', required=True, help='Clade metadata tsv file.')
    parser.add_argument('--output', '-o', required=True, help='Updated PB Tree')
    parser.add_argument('--threads', '-T', required=False, default=2, help='VCF input file.')
    args = parser.parse_args(argv)

    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    # Check Usher is installed
    check_dependency_installed('usher')
    
    # Generate query vcf file
    vcf = generate_query_vcf(args.reference, args.output, args.threads, args.sequences)
    
    # Add vcf samples to existing Usher tree
    updated_tree_fp = add_vcf_to_tree(args.pb_file, vcf, args.threads, args.output)

    # Annotate Tree with clades
    apt.annotate_tree(updated_tree_fp, args.metadata, args.output, args.threads)

    # Write all possible summary output files to a specific directory
    summary_samples_tsv, summary_output_folder = write_summary_output_files(args.pb_file, args.output)

    # Calculate uncertaintly for samples added to tree
    calculate_sample_uncertainty(args.pb_file, summary_samples_tsv, args.output)

if __name__ == "__main__":
    exit(main())
