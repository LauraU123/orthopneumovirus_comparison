from Bio import Entrez, SeqIO
import argparse

"""THIS SCRIPT DOWNLOADS GBK FILES FROM THE FOLLOWING ID LIST"""


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
        
