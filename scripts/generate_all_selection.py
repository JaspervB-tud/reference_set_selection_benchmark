from Bio import SeqIO
import os
import argparse

def find_genomes_per_species(parent_folder_path):
	'''
	This function will search "parent_folder_path" and for every subfolder find all files ending in .fa, storing them in a dictionary with species as keys.
	Note that this function assumes that the subfolders of the parent folder are formatted as "species_X" where X refers to the taxid of the species.
	Input:
		- parent_folder_path (str): path to the folder that contains subfolders for all species
	Output:
		- files_per_species (dict[str] -> list[str]): dict with all genome file names for all species (keys)
		- accessions_per_species (dict[str] -> list[str]): dict with all sequence ids for all genomes for all species (keys)
	'''
	files_per_species = {}
	accessions_per_species = {}
	for folder in os.listdir(parent_folder_path):
		subfolder_path = os.path.join(parent_folder_path, folder)
		if os.path.isdir(subfolder_path): #only consider folders
			cur_species = folder.split("_")[-1]
			files_per_species[cur_species] = []
			accessions_per_species[cur_species] = []
			for file in os.listdir(subfolder_path): #iterate over .fa files and retrieve all sequence ids for accessionid2taxid
				if file.endswith(".fa"):
					files_per_species[cur_species].append(file)
					accessions_per_species[cur_species] += [record.id for record in SeqIO.parse(os.path.join(subfolder_path, file), "fasta")]

	return files_per_species, accessions_per_species

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--genomes", type=str, required=True)
	parser.add_argument("--output", type=str, required=True) #this should be the folder to store the output
	args = parser.parse_args()

	files_per_species, accessions_per_species = find_genomes_per_species(args.genomes)
	# Write accession2taxid file to be used by profilers
	with open(f"{args.output}/nucl_gb.accession2taxid", "w") as f_out:
		f_out.write("accession\t accession.version\ttaxid\tgi\n")
		for species in accessions_per_species: #iterate over all species, storing the sequence ids and corresponding taxids
			for seq_id in accessions_per_species[species]:
				#Note that in general, sequences in NCBI have an id structured as id.version, so we split it
				cur_id = seq_id.split(".")
				if len(cur_id) == 1:
					cur_id = seq_id
				else:
					cur_id = ".".join(cur_id[:-1])
				print(f"{cur_id}\t{seq_id}\t{species}\t0")
				f_out.write(f"{cur_id}\t{seq_id}\t{species}\t0\n")
	# Write file for "all" genomes
	with open(f"{args.output}/all_selection.tsv", "w") as f_out:
		for species in files_per_species:
			for file in files_per_species[species]:
				"""
				Here we generate a file called "all_selection.txt" which will contain all genome names organized per species.
				This file will have the following structure on every line:
					$SPECIES \t $FILENAME \t +/-
				The final column will be a '+' if a selection was made for the species, and a '-' otherwise
				"""
				f_out.write(f"{species}\t{file}\t+\n")
