import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="graph a comparison of human and bovine mutation distances",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--gene', required=True, type=str, help="gene of interest")
    args = parser.parse_args()

    id_list = ['KJ627301.1',' MG431251.1', 'ON152648.1','NC_038272.1', 'NC_001989.1', 'OM328114.1', 'OM328115.1', 'MT861050.1', 'MG947594.1', 'AF295543.1', 'AF092942.1','MT107528']
    sequencesaa, sequencesnuc = ([] for i in range(2))
    for i in id_list:
        data = SeqIO.read(f'/home/laura/bovine_rsv_cleaned_up/downloaded_data/{i}_data.gbk', 'genbank')

        for feature in data.features:
            if feature.type == 'CDS':
                if feature.qualifiers['gene'][0] == args.gene:
                    #nucl sequence
                    loc_nuc = feature.location
                    start = list(loc_nuc)[0]
                    end = list(loc_nuc)[-1]
                    sequencenuc = Seq(data.seq[start:end])
                    recordnuc = SeqRecord(sequencenuc, id =data.id, description= data.description)
                    # aa sequence
                    aaasequence = (feature.qualifiers['translation'][0])
                    sequenceaa = Seq(aaasequence)
                    print(aaasequence)
                    recordaa = SeqRecord(sequenceaa, id=data.id, description= data.description)

                    sequencesnuc.append(recordnuc)
                    sequencesaa.append(recordaa)

    SeqIO.write(sequencesnuc, f'results/aa_and_nuc/{args.gene}_nucleotides.fasta', 'fasta')
    SeqIO.write(sequencesaa, f'results/aa_and_nuc/{args.gene}_aminoacids.fasta','fasta')