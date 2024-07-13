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


def gff2gtf(gff_fp):
    '''
    Convert GFF to GTF file
    :params: GFF file path
    :return: GTF file path.
    '''
    check_dependency_installed('agat_convert_sp_gff2gtf.pl')
    gtf_fn = str('.'.join(os.path.basename(gff_fp).split('.')[:-1])+('.gtf'))
    gtf_fp = os.path.join(os.path.dirname(gtf_fp), gtf_fn)
    gtf_convert_process = subprocess.Popen(['agat_convert_sp_gff2gtf.pl', '--gff', gff_fp, '-o', gtf_fp])
    gtf_convert_process.wait()

def main(argv: Optional[Sequence[str]] = None) -> int:
    '''
    Main function for running script directly
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff', '-g', required=True, help='GFF file to convert into GTF file.')
    parser.add_argument('--output', '-o', required=True, help='Output folder to save GTF file')
    args = parser.parse_args(argv)

    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    reference_index = generate_reference_index(args.reference, args.output, args.threads)
    sam_file_path = generate_alignment(args.sequences, reference_index, args.threads)
    aln_filename = sam_to_aln(sam_file_path, args.reference)
    faToVCF(aln_filename)

if __name__ == '__main__':
    exit(main())