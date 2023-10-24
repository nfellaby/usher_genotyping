# Designing a Usher Pipeline for annotation of genotype

## Starting requirements: Software
	- gbmunge: generate fasta and tsv from Genbank file (for reference)
	- minimap: generate alignment
	- fasta_to_vcf.py (custom script): generate fasta alignment from SAM file (needs renaming)
	- faToVcf: create vcf from aligned fastas
	- Usher: Generate tree
	- matUtils: Infer clade/genotype
	- agat: generate GTF amino acid mutation files (optional)

## Starting input files:
	- Reference (genbank)
	- sequences (fasta)
	- tsv with known genotypes

## Usage
If starting with no tree:
- Provide a reference GenBank file, fasta sequences to create tree from, and existing genotypes to annotate against the tree
`python main.py -r reference.gb -o output_directory -s sequences.fasta -g genotypes.csv`
- If you want to add samples to existing tree, provide the tree file (in protobuf format), reference fasta, sequences to add, genotype information, output folder
`python add_samples_to_tree.py -p tree.pb -r reference.fasta -o output_dir -s sequences_to_add.fasta -g genotypes.csv`

## To Do
- Clean up extra files formatting
- Clean up logging
- Handling no genotype information
- Add snp-dists component?
