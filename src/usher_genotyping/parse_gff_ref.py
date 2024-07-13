import argparse
import sys
import os
from usher_genotyping.src.logger import log
import glob
import logging
from typing import Optional
from typing import Sequence
import gffutils
# from usher_genotyping.src.logger import log
import os
import subprocess


'''
Parsing Genbank reference files used in Usher Tree generation
Functions Required:
1. Read in GFF file
2. Generate fasta and tsv file outputs
3. Generate GTF mutation annotate file
'''

def parse_gff_file(file_name):
    '''
    Take genbank filename and generate genkbank object
    :params: file name as string
    :return: BioSeq Genbank Object
    '''
    try:
        # Run Usher
        usher_process = subprocess.Popen(["agat_convert_sp_gff2gtf.pl", "-t", seed_nwk, "-v", vcf_file, '-o', pb_file_path, '-d', working_dir, '-T', str(threads)])
        usher_process.wait()
        # log('log', f"======= File {file_name} read in as GenBank File")
        # return gb_object
    except:
        # log('critical', f"======= File {file_name} is not a GenBank File or the Genbank file is corrupted.")
        raise 

def main(argv: Optional[Sequence[str]] = None) -> int:
    '''
    Main function for running script directly
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff_file', '-g', required=True, help='Genkbank File for reference used to make phylogenetic tree.')
    parser.add_argument('--output','-o', required=True, help='Output directory for reference files.')
    args = parser.parse_args(argv)

    if args.gff_file == '':
        logging('critical', f"======= requires genbank file to be specified")

    gff_object = parse_gff_file(args.gff_file)
    # gb_object = check_genbank_attributes(gb_object)
    # run_gbmunge(args.genbank. args.output)

if __name__ == '__main__':
    exit(main())
