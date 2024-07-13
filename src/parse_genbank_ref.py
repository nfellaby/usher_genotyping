import argparse
from typing import Optional
from typing import Sequence
from usher_genotyping.src.logger import log
from Bio import SeqIO
import os
import subprocess

'''
Parsing Genbank reference files used in Usher Tree generation
Functions Required:
1. Read in Genbank file
2. Run invoke and run gbmunge
3. Generate fasta and tsv file outputs
4. Generate GTF2 mutation annotate file
'''

def parse_genbank_file(file_name):
    '''
    Take genbank filename and generate genkbank object
    :params: file name as string
    :return: BioSeq Genbank Object
    '''
    try:
        gb_object = SeqIO.read(open(file_name), 'gb')
        log('log', f"======= File {file_name} read in as GenBank File")
        return gb_object
    except:
        log('critical', f"======= File {file_name} is not a GenBank File or the Genbank file is corrupted.")
        raise 
    
def check_genbank_attributes(gb_object):
    '''
    Check GenBank object has the expected attributes
    :params: genbank object
    '''
    # Critical Attributes
    # try except loop? or sys.exit() for failing appropriately?
    # reverse logic?
    if hasattr(gb_object, 'seq') is False:
        log('critical', f"======= Genbank file does not have sequence data attribute")
        raise
    if hasattr(gb_object, 'id') is False:
        log('critical', f"======= Genbank file does not have sequence id attribute")
        raise
    if hasattr(gb_object, 'annotations') is False:
        log('critical', f"======= Genbank file  does not have annotations attribute")
        raise

    log('info', f"Genbank object read. Has all expected attributes")
    print('Genbank File appears as expected')

def run_gbmunge(genkbank_file, output_dir):
    ''''
    Running gbmunge to generate fasta and metadata csv files from genbank object
    :params: genbank file name, environment name (which should be pathogen-protobuf)
    :return:
    '''
    output_tsv = str('.'.join(os.path.basename(genkbank_file).split('.')[:-1])+('.metadata.tsv'))
    output_tsv = os.path.join(output_dir, output_tsv)
    
    output_fasta = str('.'.join(os.path.basename(genkbank_file).split('.')[:-1])+('.fasta'))
    output_fasta = os.path.join(output_dir, output_fasta)

    # check if output files already exist
    if (os.path.isfile(output_tsv) and os.path.isfile(output_fasta)):
        log('warning', f"======= Fasta and Metadata output files already exist. These files will be over-written.")
    else:
        log('info', f"======= Genbank Fasta and Metadata output files will be written as:{output_fasta}\n{output_tsv}\n")

    results = subprocess.run(["gbmunge", "-i", genkbank_file, "-f", output_fasta, "-o", output_tsv ], 
                            capture_output=True, text=True)
    
    # Check if gbmunge runs (not output) or failed 'output'
    if results.stderr != '':
        log('critical', f"Error running gbmunge. Please check error message and try again:\n{results.stderr}")
        raise

    # Check output files are present
    if not (os.path.isfile(output_tsv) and os.path.isfile(output_fasta)):
        log('critical', f"======= Expected output metadata ({output_tsv}) and fasta ({output_fasta}) file not generated from Genbank.")
        raise
    else:
        log('info', f'======= Fasta and metadata files generated from Genbank file.')
        print(f"Fasta and metadata generated from Genbank file.")


def main(argv: Optional[Sequence[str]] = None) -> int:
    '''
    Main function for running script directly
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--genbank', '-g', required=True, help='Genkbank File for reference used to make phylogenetic tree.')
    parser.add_argument('--output','-o', required=True, help='Output directory for reference files.')
    args = parser.parse_args(argv)

    if args.genbank == '':
        logging('critical', f"======= requires genbank file to be specified")

    gb_object = parse_genbank_file(args.genbank)
    gb_object = check_genbank_attributes(gb_object)
    run_gbmunge(args.genbank. args.output)

if __name__ == '__main__':
    exit(main())
