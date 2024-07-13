from typing import Optional
from typing import Sequence
import argparse
from usher_genotyping.src.logger import log
import os
import subprocess
from usher_genotyping.src.sam_to_aln import aln_to_fasta


"""
Mutliple sequence alignment
- input: "{sample}.fasta"
- output: "{sample}.sam.vcf"

- Create index:
    1. `minimap2 -t 2 -d ./A.goose.Guangdong.1.1996_HA.mni A.goose.Guangdong.1.1996_HA.fasta`
- Alignment:
    `minimap2 -t 2 --score-N=0 --secondary=no --sam-hit-only -a -o ./HA_all_sequences.filt.sam A.goose.Guangdong.1.1996_HA.mni HA_all_sequences.filt.no_ref.fasta`
- Convert SAM to FASTA (using bbmap)
    `python ../fasta_to_vcf.py -a HA_all_sequences.filt.sam -r A.goose.Guangdong.1.1996_HA.fasta -o ./`

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

def generate_reference_index(reference_file_name, output, threads):
    """
    Use minimap2 to generate a n index
    :params: reference fasta file path
    :return: reference index file path
    """
    reference_index = str('.'.join(os.path.basename(reference_file_name).split('.')[:-1])+('.mni'))
    ref_index_fn = os.path.join(output, reference_index)

    # check if output files already exist
    if os.path.isfile(ref_index_fn):
        log('warning', f"======= Index file for reference already exists. This files will be over-written.")
    else:
        log('info', f"======= Generating Index file:{ref_index_fn}")

    # check minimap2 is installed
    check_dependency_installed('minimap2')

    # run minimap with reference file
    minimap_process = subprocess.Popen(["minimap2", "-t", str(threads), "-d", ref_index_fn, reference_file_name])
    minimap_process.wait()

    #TODO: handle the stdout stderr better for this subprocess.
    # check output files exist
    if os.path.isfile(ref_index_fn):
        log("info", f"======= Minimap2 Index generated: {ref_index_fn}")
        return ref_index_fn
    else:
        log("critical", f"Minimap2 failed to produce expected output file {ref_index_fn}")
    

def generate_alignment(sequences, reference_index, threads):
    """
    Generate alignment from multi-fasta sequence file and reference index
    :param: multi-fasta sequence filename
    :return: index reference file name
    """
    sam_filename = str('.'.join(os.path.basename(sequences).split('.')[:-1])+('.sam'))
    sam_file_path = os.path.join(os.path.dirname(reference_index), sam_filename)
    alignment_process = subprocess.Popen([
        "minimap2", 
        "-t", str(threads), 
        "--score-N=0", 
        "--secondary=no", 
        "--sam-hit-only", 
        "-a", 
        "-o", sam_file_path, 
        reference_index, 
        sequences
    ])
    alignment_process.wait()

    if os.path.isfile(sam_file_path):
        log("info", f"======= SAM Alignment file generated: {sam_file_path}")
        return sam_file_path
    else:
        log('critical', f"======= SAM alignment file not generated. Please check log file:\n{alignment_process.stderr}")

def sam_to_aln(out_aln_path, ref_fasta_fn):
    '''
    Use the sam_to_aln.py script to generate an alignment file from the sam file
    # python ../fasta_to_vcf.py -a HA_all_sequences.filt.sam -r A.goose.Guangdong.1.1996_HA.fasta -o ./
    :params: sam file name, reference fasta
    :return: alignment file
    '''
    out_msa_fn =   str('.'.join(os.path.basename(out_aln_path).split('.')[:-1])+('.aln'))
    out_msa_path =  os.path.join(os.path.dirname(out_aln_path), out_msa_fn)
    omit_ref = True
    buffer_size = int(1048576)
    aln_to_fasta(out_aln_path, out_msa_path, ref_fasta_fn, bufsize=buffer_size)

    # check output files are generated
    if os.path.isfile(out_msa_path):
        log('info',f'======= Alignment file generated: {out_msa_path}')
    else:
        log('critical', f"======= alignment file not generated. Please check input files.")
    
    return out_msa_path

def faToVCF(out_msa_path):
    '''
    Convert alignment file to VCF
    :params: alignment file
    :return: vcf file
    '''
    vcf_fh =  str('.'.join(os.path.basename(out_msa_path).split('.')[:-1])+('.vcf'))
    vcf_file_path = os.path.join(os.path.dirname(out_msa_path), vcf_fh)
    check_dependency_installed('faToVcf')
    vcf_process = subprocess.Popen(['faToVcf', out_msa_path, vcf_file_path])
    vcf_process.wait()
    
    # check output files are generated
    if os.path.isfile(out_msa_path):
        log('info',f'======= VCF Alignment file generated: {out_msa_path}')
    else:
        log('critical', f"======= VCF alignment file not generated. Please check input files.")
    
def main(argv: Optional[Sequence[str]] = None) -> int:
    '''
    Main function for running script directly
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--sequences', '-s', required=True, help='Contextual FASTA sequences used to make phylogenetic tree.')
    parser.add_argument('--reference', '-r', required=True, help='Reference FASTA used to make phylogenetic tree.')
    parser.add_argument('--output', '-o', required=True, help='Output directory to save vcf')
    parser.add_argument('--threads', '-T', required=False, default=2, help='VCF input file.')
    args = parser.parse_args(argv)

    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    reference_index = generate_reference_index(args.reference, args.output, args.threads)
    sam_file_path = generate_alignment(args.sequences, reference_index, args.threads)
    aln_filename = sam_to_aln(sam_file_path, args.reference)
    faToVCF(aln_filename)

if __name__ == '__main__':
    exit(main())