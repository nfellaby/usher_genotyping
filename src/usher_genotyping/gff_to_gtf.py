from typing import Optional
from typing import Sequence
import argparse
from usher_genotyping.src.logger import log
import os
import subprocess

"""
Create GTF file from GFF file
:params: GFF file
:return: GTF file
"""


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
        print(f"{dependency} not found in 'which' path. Check it has been installed.")

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

def aa_annotation(pb_tree, ref_fasta, gtf, output):
    '''
    Convert GFF to GTF file
    :params: GFF file path
    :return: GTF file path.
    '''
    check_dependency_installed('agat_convert_sp_gff2gtf.pl')
    aa_fn = str('.'.join(os.path.basename(pb_tree).split('.')[:-1])+('.aa_annotations.tsv'))
    aa_fp = os.path.dirname(output, aa_fn)
    annotate_aa_process = subprocess.Popen(
        [ 'matUtils summary',
            '--translate', aa_fp,
            '-i', pb_tree,
            '-r', 'NC_007362_HA_ref.gtf',
            '-f', ref_fasta
        ])
    annotate_aa_process.wait()
    check_output_file(aa_fp)
    return aa_fp

def main(argv: Optional[Sequence[str]] = None) -> int:
    '''
    Main function for running script directly
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--pb_tree', '-p', required=True, help='Protobuf Tree file.')
    parser.add_argument('--gtf', '-g', required=True, help='GFT File.')
    parser.add_argument('--reference', '-r', required=True, help='Reference fasta file.')
    parser.add_argument('--output', '-o', required=True, help='Output folder to save GTF file')
    args = parser.parse_args(argv)


    aa_fp = aa_annotation(args.pb_tree, args.ref_fasta, args.gtf, args.output)


if __name__ == '__main__':
    exit(main())