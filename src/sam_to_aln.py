#!/usr/bin/env python

"""Convert fasta files to VCF files.

This script converts a fasta file and a given fasta reference file to
a file in VCF. A reference genome can be provided in fasta format as
input. If no reference is given, the first sequence in the fasta file
will be used as reference.

"""

import argparse
from multiprocessing import cpu_count
import sys
import os
from hashlib import md5
from shutil import copy, move
from os.path import abspath, expanduser, isdir, isfile, split



DEFAULT_BUFSIZE = 1048576 # 1 MB #8192 # 8 KB
DEFAULT_ALIGNER = 'minimap2'
DEFAULT_THREADS = cpu_count()
VERSION = 0.1
CIGAR_LETTERS = {'M','D','I','S','H','=','X'}


# count the number of IDs in a FASTA file
def count_IDs_fasta(fn, bufsize=DEFAULT_BUFSIZE):
    return sum(l.startswith('>') for l in open(fn, buffering=bufsize))

# parse a CIGAR string
def parse_cigar(s):
    out = list(); ind = len(s)-1
    while ind >= 0:
        let = s[ind]; ind -= 1; num = ''
        while s[ind] not in CIGAR_LETTERS:
            num += s[ind]; ind -= 1
        out.append((let, int(num[::-1])))
    return out[::-1]

# parse user args
def parse_args():
     # use argparse to parse user arguments
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-a', '--out_aln_path', required=True, type=str, help="Aligned sequences (SAM format)")
    parser.add_argument('-r', '--reference', required=True, type=str, help="Reference")
    parser.add_argument('-o', '--output', required=True, type=str, help="Output Directory")
    parser.add_argument('-t', '--threads', required=False, type=int, default=DEFAULT_THREADS, help="Number of Threads")
    parser.add_argument('--omit_ref', action="store_true", help="Omit reference sequence from output alignment")
    parser.add_argument('-b', '--buffer_size', required=False, type=int, default=DEFAULT_BUFSIZE, help="File Stream Buffer Size (bytes)")
    args = parser.parse_args()

    # check user args for validity
    if args.threads < 1:
        print("ERROR: Number of threads must be positive", file=sys.stderr); exit(1)
    if args.buffer_size < 1:
        print("ERROR: Output buffer size must be positive", file=sys.stderr); exit(1)
    args.output = abspath(expanduser(args.output))
    if not os.path.isdir(args.output):
        os.path.mkdir(args.output)
    if isfile(args.reference):
        if count_IDs_fasta(args.reference, bufsize=args.buffer_size) != 1:
            print("ERROR: Reference file (%s) must have exactly 1 sequence in the FASTA format" % args.reference, file=sys.stderr); exit(1)

    args.ref_genome_path = abspath(expanduser(args.reference))

    return args



# convert alignment (SAM/PAF) to FASTA
def aln_to_fasta(out_aln_path, out_msa_path, ref_genome_path, bufsize=DEFAULT_BUFSIZE, omit_ref=False):
    if out_aln_path.lower().endswith('.sam'):
        aln_type = 's' # SAM
    elif out_aln_path.lower().endswith('.paf'):
        aln_type = 'p' # PAF
    else:
        print("ERROR: Invalid alignment extension: %s" % out_aln_path, file=sys.stderr); exit(1)
    msa = open(out_msa_path, 'w', buffering=bufsize); ref_seq = list()
    for line in open(ref_genome_path):
        if len(line) == 0:
            continue
        if line[0] != '>':
            ref_seq.append(line.strip())
        elif not omit_ref:
            msa.write(line)
    if not omit_ref:
        for l in ref_seq:
            msa.write(l)
        msa.write('\n')
    ref_seq_len = sum(len(l) for l in ref_seq)
    num_output_IDs = 0
    for l in open(out_aln_path):
        if l == '\n' or l[0] == '@':
            continue
        parts = l.split('\t')
        if aln_type == 's': # SAM
            flags = int(parts[1])
            if flags != 0 and flags != 16:
                continue
            ID = parts[0].strip(); num_output_IDs += 1
            ref_ind = int(parts[3])-1
            seq = parts[9].strip()
            cigar = parts[5].strip()
        elif aln_type == 'p': # PAF
            raise RuntimeError("PAF alignments are not yet supported")
        edits = parse_cigar(cigar)
        msa.write(">%s\n" % ID)
        if ref_ind > 0:
            msa.write('-'*ref_ind) # write gaps before alignment
        ind = 0; seq_len = ref_ind
        for e, e_len in edits:
            if e == 'M' or e == '=' or e == 'X': # (mis)match)
                msa.write(seq[ind:ind+e_len])
                ind += e_len; seq_len += e_len
            elif e == 'D':                       # deletion (gap in query)
                msa.write('-'*e_len)
                seq_len += e_len
            elif e == 'I':                       # insertion (gap in reference; ignore)
                ind += e_len
            elif e == 'S' or e == 'H':           # starting/ending segment of query not in reference (i.e., span of insertions; ignore)
                ind += e_len
        if seq_len < ref_seq_len:
            msa.write('-'*(ref_seq_len-seq_len)) # write gaps after alignment
        msa.write('\n')
    msa.close()
    return num_output_IDs


# main content
def main():
    # parse user args and prepare run
    args = parse_args()
    print(args)

    # align viral genomes against reference
    out_aln_path = args.out_aln_path

    # convert alignment (SAM/PAF) to MSA FASTA
    out_msa_path = '%s/%s.aln' % (args.output, args.out_aln_path.split('/')[-1])

    num_output_IDs = aln_to_fasta(out_aln_path, out_msa_path, args.ref_genome_path, omit_ref=args.omit_ref, bufsize=args.buffer_size)


# run tool
if __name__ == "__main__":
    main()
