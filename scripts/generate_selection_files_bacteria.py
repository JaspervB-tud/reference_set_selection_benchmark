from Bio import SeqIO
import os
import argparse
import random

# Scripts for gathering genomes
def read_standard_output(selection_path):
	"""
	This covers selections made by: hierarchical, centroid, ggrasp
	"""
	selection = []
	with open(selection_path, "r") as f_in:
		for line in f_in:
			selection.append(line.strip())
	return selection

def read_gclust_output(selection_path):
	selection = []
	with open(selection_path, "r") as f_in:
		for line in f_in:
			line = line.strip().split()
			if line[-1] == "*": #cluster representative
				selection.append(line[2][1:-3]) #strips the leading ">" and trailing "..."
	return selection

def read_meshclust_output(selection_path):
	selection = []
	with open(selection_path, "r") as f_in:
		for line in f_in:
			line = line.strip()
			if len(line) > 0:
				line = line.split("\t")
				if line[-1] == "C": #cluster representative
					selection.append(line[1][1:])
	return selection

def read_vsearch_output(selection_path):
	selection = [record.id for record in SeqIO.parse(selection_path, "fasta")] #vsearch outputs fasta file
	return selection


def find_genomes_per_species(parent_folder_path):
	'''
	This function will search "parent_folder_path" and for every subfolder find all files ending in .fa, storing them in a dictionary with species as keys.
	Note that this function assumes that the subfolders of the parent folder are formatted as "species_X" where X refers to the taxid of the species.
	Input:
		- parent_folder_path (str): path to the folder that contains subfolders for all species
	Output:
		- files_per_species (dict[str] -> list[str]): dict with all genome file names for all species (keys)
	'''
	files_per_species = {}
	for folder in os.listdir(parent_folder_path):
		subfolder_path = os.path.join(parent_folder_path, folder)
		if os.path.isdir(subfolder_path): #only consider folders
			cur_species = folder.split("_")[-1]
			files_per_species[cur_species] = []
			for file in os.listdir(subfolder_path): #iterate over .fa files and retrieve all sequence ids for accessionid2taxid
				if file.endswith(".fa"):
					files_per_species[cur_species].append(file)

	return files_per_species

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--genomes", type=str, required=True)
	parser.add_argument("--selection", type=str, required=True) #path to folder that contains the selection
	parser.add_argument("--filename", type=str, required=True) #naming convention for selection (e.g. ggrasp, or complete-linkage_0.01)
	parser.add_argument("--output", type=str, required=True) #this should be the folder to store the output
	args = parser.parse_args()

	files_per_species = find_genomes_per_species(args.genomes) #these are used to determine which species only have 1 genome

	selection_per_species = {}
	for species in files_per_species:
		cur_folder_path = f"{args.selection}/{species}"
		selection_per_species[species] = []
		if len(files_per_species[species]) <= 1: #up to 1 genome for the species -> automatically select
			selection_per_species[species] = [(files_per_species[species][0], "-")] #we set the "selected" value (+/-) to - which is used when comparing reference sets
		else:
			# Check if output file exists, otherwise randomly select
			file_path = f"{cur_folder_path}/{args.filename}"
			print(file_path)
			if os.path.exists(file_path):
				# Check which method
				if "gclust" in args.selection:
					cur_selection = read_gclust_output(file_path)
				elif "meshclust" in args.selection:
					cur_selection = read_meshclust_output(file_path)
				elif "vsearch" in args.selection:
					cur_selection = read_vsearch_output(file_path)
				else:
					cur_selection = read_standard_output(file_path)
				# Store selection
				for selected in cur_selection:
					selection_per_species[species].append( (selected, "+") )
			else: #file does not exist, add 
				selection_per_species[species] = [(random.choice(files_per_species[species]), "-")]
	with open(f"{args.output}/{args.filename}", "w") as f_out:
		for species in selection_per_species:
			for file, val in selection_per_species[species]:
				f_out.write(f"{species}\t{file}\t{val}\n")
	print(selection_per_species)



	'''
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
	'''

if __name__ == "__main__":
	main()