GENE = ['NS1', 'NS2', 'P', 'N', 'G', 'F', 'M', 'L']
IDS = ['KJ627301.1','MG431251.1', 'ON152648.1','NC_038272.1', 'NC_001989.1', 'OM328114.1', 'OM328115.1', 'MT861050.1', 'MG947594.1', 'AF295543.1', 'AF092942.1','MT107528']


rule all:
    input:
        graph = "results/graph.pdf",
        ids = expand('data/downloaded_data/{ids}_data.gbk', ids=IDS)

rule save:
    message:
    	"""
    	saving files from Genbank
    	"""
    params:
        id_ = lambda w: w.ids,
    output:
        file = 'data/downloaded_data/{ids}_data.gbk'
    shell:
        """
        python code/save.py \
        --file {output.file} \
        --ids {params.id_}
        """
  
rule sequences:
    message:
        """
        dividing into amino acid sequences and nucleotide sequences for each gene
        """
    input:
        files = expand('data/downloaded_data/{ids}_data.gbk', ids=IDS)
    output: 
        aa = "results/aa_and_nuc/{gene}_aminoacids.fasta",
        nuc = "results/aa_and_nuc/{gene}_nucleotides.fasta",
    params:
        gene = lambda w: w.gene
    shell:
        """
        python code/saving_files.py \
        --gene {params.gene}
        """

rule align:
    input:
        aminoacids = rules.sequences.output.aa,
        nucleotides = rules.sequences.output.nuc
    output:
        alignedaa = "results/alignedfiles/{gene}_alignedaa.fasta",
        alignednuc = "results/alignedfiles/{gene}_alignednuc.fasta"
    shell:
        """
        augur align \
        --sequences {input.aminoacids} \
        --output {output.alignedaa}

        augur align \
        --sequences {input.nucleotides} \
        --output {output.alignednuc}
        """

rule treesaa:
    input:
        aa = rules.align.output.alignedaa
    output:
        treefile = "results/alignedfiles/{gene}_alignedaa.fasta.treefile",
        tmpfiles = temp(f"results/alignedfiles/{{gene}}_alignedaa.fasta.{x}" for x in ["ckp.gz", "bionj", "mldist","log"])
    shell:
        """
        iqtree -s {input.aa} -m LG
        """

rule treesnuc:
    input:
        nuc = rules.align.output.alignednuc
    output:
        treefile = "results/alignedfiles/{gene}_alignednuc.fasta.treefile",
        tmpfiles = temp(f"results/alignedfiles/{{gene}}_alignednuc.fasta.{x}" for x in ["ckp.gz", "bionj", "mldist","log", "model.gz"])
    shell:
        """
        iqtree -s {input.nuc}
        """


def get_trees_nuc():
    all_trees = ""
    genes = ['NS1', 'NS2', 'P', 'N', 'G', 'F', 'M', 'L']
    for gene_ in genes:
        all_trees+=f"results/alignedfiles/{gene_}_alignednuc.fasta.treefile+"
    return(all_trees)

def get_trees_aa():
    all_trees = ""
    genes = ['NS1', 'NS2', 'P', 'N', 'G', 'F', 'M', 'L']
    for gene_ in genes: 
        all_trees+=f"results/alignedfiles/{gene_}_alignedaa.fasta.treefile+"
    print(all_trees)
    return(all_trees)

rule graph:
    input:
        treesaa = expand(rules.treesnuc.output.treefile, gene=GENE),
        treesnuc = expand(rules.treesaa.output.treefile, gene=GENE)
    output:
        tsv = "results/dataframe.tsv",
        graph = "results/graph.pdf"
    params:
        treesaa = get_trees_nuc(),
        treesnuc = get_trees_aa()
    shell:
        """
        python code/graph.py \
        --aatrees {params.treesaa} \
        --nttrees {params.treesnuc} \
        --dataframe {output.tsv} \
        --output {output.graph}
        """




