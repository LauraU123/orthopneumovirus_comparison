# Orthopneumovirus Comparison

This workflow compares the distances between human and bovine RSV for each gene. 
It downloads RSV-A, RSV-B and bovine sequences from NCBI GenBank and outputs a graph with
distances between clades in each group. 

### Running the Workflow

To run the workflow:

```
snakemake --cores all
```

### Workflow steps

* GenBank files of the sequences of interest are downloaded and saved from GenBank

* Nucleotide sequences and amino acid sequences are split by gene

* Each gene is aligned separately using MSA (separate for amino acids and nucleotides)

* Phylogenetic trees are constructed for each gene for amino acids and nucleotides

* The distances between each clade of interest are calculated and a graph is constructed (distances are based on summing branch lengths between the clades of interest)


### Results

The output of this workflow is a graph in pdf format comparing human RSV subtypes, human and bovine RSV and bovine diversity. 
 ![RSV](https://github.com/LauraU123/genome_coverage/blob/master/example_graphs/a_coverageGraph.png)

