# Testing:
- To run the CLI, from the top directory: `python3 -m usher_genotyping.command_line`
- To run the create_tree command, from the top directory: `python3 -m usher_genotyping.command_line create_tree`


# Current Questions
- Do we need an seperate scripts for querying additional sequences against a known pb file with clades already annotated to them

# Requirements

## Starting requirements: Software
	- gbmunge: generate fasta and tsv from Genbank file (for reference)
	- minimap: generate alignment
	- fasta_to_vcf.py (custom script): generate fasta alignment from SAM file (needs renaming)
	- faToVcf: create vcf from aligned fastas
	- Usher: Generate tree
	- matUtils: Infer clade/genotype

## Starting input files:
	- Reference (genbank)
	- sequences (fasta)
	- tsv with known genotypes


## Pipeline Staging
0. Create test directory:
- reference
- fasta sequences
- tsv with genotypes
- samples to add to reference tree * - not added yet

1. Generating Reference Usher Tree and Clade Annotations:
- Reference (gb)
- sequences (fasta)
- tsv with known genotypes

- Testing requirements:
	- Check pre-required software is available
	- Genbank is as expected
	- Sequences are in fasta format, and more than 1 sequence
	- tsv genotypes matches fasta sequences in dataset

2. Adding additional sequences to existing Usher Tree and estimating Clades
- Usher Tree with clade annotations
- Query sequence(s)

- Testing requirements:
	- Check pre-required software is available
	- Check Usher tree has clade annotations
	- Check query sequences are in the appropriate format
	- Sanity check estimated clades

3. Generating mutations for PB tree, requries GTF file
	- Convert GB file to GTF file
	- Add GTF file to PB tree

- Testing: `pytest -vv -rP usher_genotyping/tests/test_annotate_phylogenetic_tree.py`

