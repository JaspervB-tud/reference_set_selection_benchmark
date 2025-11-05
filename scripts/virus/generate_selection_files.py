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

def read_vlq_output(selection_path, genomes_path):
	"""
	Note that this works differently compared to the other functions since the output of VLQ
	is significantly different (i.e. it produces a fasta file with all targets, and a metadata
	file with metadata for the selected targets). The output will be the same, however
	Input: 
		- selection_path (str): path to the selection folder
		- genomes_path (str): path to the genomes folder
	Output:
		- selection (dict[str]->str): dictionary mapping lineages to included genomes
	"""
	lin2ids = {} #maps lineage to sequence ids
	with open(f"{selection_path}/metadata.tsv", "r") as f_in:
		next(f_in) #skip header
		for line in f_in:
			line = line.strip().split("\t")
			cur_id = line[0]
			cur_lineage = line[11]
			if cur_lineage not in lin2ids:
				lin2ids[cur_lineage] = set()
			lin2ids[cur_lineage].add(cur_id)

	# Now parse original genome files and determine which sequence ids correspond to which genomes
	id2seq = {} #maps sequence id to genome filename
	selection = {lineage: [] for lineage in lin2ids}
	for folder in os.listdir(genomes_path):
		subfolder_path = os.path.join(genomes_path, folder)
		if os.path.isdir(subfolder_path): #only consider folders
			cur_lineage = folder.split("_")[-1]
			for file in os.listdir(subfolder_path): #iterate over .fa files and retrieve all sequence ids
				if file.endswith(".fa"):
					cur_ids = [record.id for record in SeqIO.parse(os.path.join(subfolder_path, file), "fasta")]
					for cur_id in lin2ids[cur_lineage]:
						if cur_id in cur_ids:
							if len([f for f in os.listdir(subfolder_path) if f.endswith(".fa")]) > 1: #this is for comparing against other reference sets
								selection[cur_lineage].append((file, "+"))
							else:
								selection[cur_lineage].append((file, "-")) 
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


	selection_per_species = {}
	if not "vlq" in args.filename:
		files_per_species = find_genomes_per_species(args.genomes) #these are used to determine which species only have 1 genome
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
				else: #file does not exist, add randomly
					selection_per_species[species] = [(random.choice(files_per_species[species]), "-")]
	else:
		selection_per_species = read_vlq_output(args.selection, args.genomes)
		
	with open(f"{args.output}/{args.filename}.tsv", "w") as f_out:
		for species in selection_per_species:
			for file, val in selection_per_species[species]:
				f_out.write(f"{species}\t{file}\t{val}\n")
	print(selection_per_species)
	

if __name__ == "__main__":
	main()