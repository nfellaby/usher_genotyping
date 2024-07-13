from typing import Optional
from typing import Sequence
import argparse
from usher_genotyping.src.logger import log
import os
import subprocess
import pandas
from usher_genotyping.src.multi_sequence_alignment import check_dependency_installed

"""
Build Phylogenetic Tree with Usher

Build tree (if one doesn't exist)
- input: "{sample}.vcf", "seed.nwk" is an "empty" tree file (just contains "();")
- output: {sample}.pb

`usher -i {config[input_pb]} -v {input} -o {output}`

Annotate Tree with known clades/genotyping information

matUtils annotate -i HA_all_sequences.filt.pb -o HA_all_sequences.filt.clade_annotated.pb -c 20230723_ha_genotypes.v2.csv

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

def build_tree(vcf_file, threads):
    '''
    Use Usher to generate phylogenetic tree
    :params: vcf file
    :returns: phylognetic tree file
    '''
    # pb filename
    pb_fn = '.'.join(os.path.basename(vcf_file).split('.vcf')[:-1])+str('.pb')
    working_dir = os.path.join(os.path.dirname(vcf_file), 'usher_out')
    if not os.path.isdir(working_dir):
        os.mkdir(working_dir)
    pb_file_path = os.path.join(working_dir, pb_fn)

    # Sort out empty seed file for Usher
    seed_nwk = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data/seed.nwk')
    if os.path.isfile(seed_nwk):
        log('info', f"======= Using {seed_nwk} as blank template to generate phylogenetic tree")
    else:
        log('warning', f"======= Expected blank template file not found where it was expected. Generating new one.")
        seed_nwk = './seed.nwk' 
        with open('./seed.nwk','w') as seed_nwk:
            pass
        log('info', f"======= {seed_nwk} has been created.")

    # Run Usher
    usher_process = subprocess.Popen(["usher", "-t", seed_nwk, "-v", vcf_file, '-o', pb_file_path, '-d', working_dir, '-T', str(threads)])
    usher_process.wait()
    
    # Check output file has been generated
    check_output_file(pb_file_path)
    return pb_file_path


def main(argv: Optional[Sequence[str]] = None) -> int:
    # read in arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_file', '-v', required=True, help='VCF input file.')
    parser.add_argument('--threads', '-T', required=False, default=2, help='VCF input file.')
    args = parser.parse_args(argv)

    # Check Usher
    check_dependency_installed('usher')

    # Build tree from VCF file
    pb_out = build_tree(args.vcf_file, args.threads)

    # TODO: Write all possible summary output files to a specific directory


if __name__ == "__main__":
  exit(main())