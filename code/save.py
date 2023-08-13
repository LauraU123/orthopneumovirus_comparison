from Bio import Entrez, SeqIO
import argparse

"""This script downloads gbk files from NCBI genbank given the id and output file"""

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="downloading genbank files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ids', required=True, type=str, help="output file")
    parser.add_argument('--file', required=True, type=str, help="output file")
    args = parser.parse_args()

    Entrez.email = 'laura.urbanska@stud.unibas.ch'
    handle = Entrez.efetch(db="nuccore", id=args.ids, rettype="gb", retmode="text")
    genome = SeqIO.parse( handle, "genbank")
    SeqIO.write(genome, args.file, 'genbank')
        
