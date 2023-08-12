import argparse
from Bio import Phylo
from os import path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="graph a comparison of human and bovine mutation distances",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--aatrees', required=True, type=str, help="aa tree file")
    parser.add_argument('--nttrees', required=True, type=str, help="nt tree file")
    parser.add_argument('--output', required=True, type=str, help="graph file")
    parser.add_argument('--dataframe', required=True, type=str, help="dataframe file")
    args = parser.parse_args()

    listofgenes = ['NS1', 'NS2', 'N','P', 'M', 'G', 'F', 'L' ]
    files_aa = args.aatrees
    files_aa = files_aa.split("+") 
    files_aa.pop(-1)
    files_nuc = args.nttrees
    files_nuc = files_nuc.split("+")
    files_nuc.pop(-1)

    print("files", files_aa, files_nuc)

    dictionary_aa, dictionaryhrsv_aa, dictionarybrsv_aa, dictionarybrsv_nuc, dictionary_nuc, dictionaryhrsv_nuc  = (dict() for i in range(6))
    distances_aa, distances_hrsv_aa, distances_brsv_aa, distances_nuc, distances_hrsv_nuc, distances_brsv_nuc = ([] for i in range(6))

    for file in files_aa:
        print(file)
        Tree= Phylo.read(file, 'newick')
        Tree.root.clades
        Tree.root_at_midpoint()
        Tree.ladderize
        hrsv_anc = Tree.common_ancestor('MT107528.1', 'KJ627301.1', 'ON152648.1')
        brsv_anc = Tree.common_ancestor('MG947594.1', 'MT861050.1', 'OM328114.1', 'AF092942.1')
        distances_aa.append(Tree.distance(hrsv_anc, brsv_anc))
        print(distances_aa)

        rsvA_anc = hrsv_anc.common_ancestor('KJ627301.1', 'ON152648.1')
        rsvB_anc = hrsv_anc.common_ancestor('MT107528.1','MG431251.1')
        distances_hrsv_aa.append(hrsv_anc.distance(rsvA_anc, rsvB_anc))
        distances_brsv_aa.append(brsv_anc.total_branch_length() - brsv_anc.branch_length)

    for file in files_nuc:
        Tree = Phylo.read(file, 'newick')
        Tree.root.clades
        Tree.root_at_midpoint()
        Tree.ladderize
        print(Tree)
        hrsv_anc = Tree.common_ancestor('MT107528.1', 'KJ627301.1', 'ON152648.1')
        print(hrsv_anc)
        brsv_anc = Tree.common_ancestor('MG947594.1', 'MT861050.1', 'OM328114.1', 'AF092942.1')
        distances_nuc.append(Tree.distance(hrsv_anc, brsv_anc))

        rsvA_anc = hrsv_anc.common_ancestor('KJ627301.1', 'ON152648.1')
        rsvB_anc = hrsv_anc.common_ancestor('MT107528.1','MG431251.1')
        distances_hrsv_nuc.append(hrsv_anc.distance(rsvA_anc, rsvB_anc))
        distances_brsv_nuc.append(brsv_anc.total_branch_length() - brsv_anc.branch_length)

    for i, j in zip(listofgenes, distances_aa): 
            print(j)
            dictionary_aa[i] = j
    for i, j in zip(listofgenes, distances_hrsv_aa): 
            print(j)
            dictionaryhrsv_aa[i] = j
    for i, j in zip(listofgenes, distances_brsv_aa): 
            print(j)
            dictionarybrsv_aa[i] = j
    for i, j in zip(listofgenes, distances_nuc): 
            print(j)
            dictionary_nuc[i] = j
    for i, j in zip(listofgenes, distances_hrsv_nuc): 
            print(j)
            dictionaryhrsv_nuc[i] = j
    for i, j in zip(listofgenes, distances_brsv_nuc): 
            dictionarybrsv_nuc[i] = j

    df = pd.DataFrame([dictionary_aa, dictionaryhrsv_aa, dictionarybrsv_aa, dictionary_nuc, dictionaryhrsv_nuc, dictionarybrsv_nuc])
    print(df)
    df.to_csv(args.dataframe, sep='\t')

    tsv_file = open(args.dataframe)
    read_tsv = pd.read_csv(tsv_file, delimiter='\t').T
    tsv_dict =read_tsv.to_dict()
    dictionaries = []

    for i, j in tsv_dict.items():
        del j['Unnamed: 0']
        dictionaries.append(j)

    labels = ["human-bovine aa", "rsv-a/b aa", "bovine diversity aa", "human-bovine nuc", "rsv-a/b nuc", "bovine diversity nuc"]
    for c, (l, i) in enumerate(zip(labels, dictionaries)):
        if "human-bovine" in l: color_ = "violet"
        if "rsv-a/b" in l: color_ = "blue"
        if "bovine diversity" in l: color_ = "darkgreen"
        if "aa" in l:  
            marker_ = "o"
            linestyle_ = "-"
        if "nuc" in l: 
            marker_ = "x"
            linestyle_ = "--"
        keys = i.keys()
        values = i.values()
        plt.plot(np.arange(len(keys))+c*0.1, values, label=l, linestyle= linestyle_, marker=marker_, color=color_)
        plt.xticks(range(len(keys)), keys)
        plt.xlabel('Gene')
        plt.ylabel('Distance')
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(args.output, bbox_extra_artists=(lgd,), format='pdf')