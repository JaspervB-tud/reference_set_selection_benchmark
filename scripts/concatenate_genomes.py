from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()

    fasta_files = [file for file in os.listdir(args.genomes) if file.endswith(".fa")] #gather all fasta files
    genomes = []
    for file in fasta_files:
        concatenated_sequence = Seq("".join(str(record.seq) for record in SeqIO.parse(f"{args.genomes}/{file}", "fasta"))) #concatenate current fasta file
        genomes.append(SeqRecord(Seq(concatenated_sequence), id=file, description=""))
    SeqIO.write(genomes, f"{args.output}/all_genomes.fasta", "fasta")

if __name__ == "__main__":
    main()
