# Orthopneumovirus Comparison

### save.py

This script downloads and saves human and bovine RSV files from NCBI genbank. Output is in gbk format

### split.py

The input genbank files for all sequences are split by gene (for nucleotides and amino acids). Outputs are in fasta format, and separate for each gene.

### graph.py

A graph comparing distances between human RSV subtypes is , bovine diversity and human and bovine RSV is constructed. 
The inputs are phylogenetic trees for each gene (separate for amino acids and nucleotides) in treefile format.
The distances between each clade of interest is calculated by summing up branch lengths.
The output is a graph png showing these distances.
