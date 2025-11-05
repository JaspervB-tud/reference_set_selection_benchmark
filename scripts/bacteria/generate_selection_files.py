from Bio import SeqIO
import os
import argparse
import random
from TaxTree import *

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

def find_species_per_taxon(path, taxtree):
	taxa = ["genus", "family", "order"]
	selected_taxa = {
		"genus": 1301, # Streptococcus
		"family": 543, # Enterobacteriaceae
		"order": 91347 # Enterobacterales
	}
	species_per_taxon = {taxon: set() for taxon in taxa}
	for taxon in taxa:
		with open(f"{path}/selected_{taxon}_species.tsv", "r") as f_in:
			for line in f_in:
				taxid = int(line.strip())
				try:
					cur_taxon = taxtree.nodes[taxid].taxonomy[taxon]
					if cur_taxon == selected_taxa[taxon]:
						species_per_taxon[taxon].add(taxid)
				except:
					continue
	return species_per_taxon
			

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--genomes", type=str, required=True)
	parser.add_argument("--selection", type=str, required=True) #path to folder that contains the selection
	parser.add_argument("--filename", type=str, required=True) #naming convention for selection (e.g. ggrasp, or complete-linkage_0.01)
	parser.add_argument("--output", type=str, required=True) #this should be the folder to store the output
	args = parser.parse_args()

	tax_tree = read_taxtree('root/taxdmp/nodes.dmp')
	tax_tree = Tree(tax_tree)

    print("Populating taxonomic tree")
    for node in tax_tree.nodes:
        tax_tree.find_lineage(tax_tree.nodes[node].taxid)
    taxid2name, taxid2syn = read_taxnames('root/taxdmp/names.dmp')
    print("Done populating taxonomic tree")

	files_per_species = find_genomes_per_species(args.genomes)
	species_per_taxon = find_species_per_taxon("root", tax_tree)

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
				else:
					cur_selection = read_standard_output(file_path)
				# Store selection
				for selected in cur_selection:
					selection_per_species[species].append( (selected, "+") )
			else: #file does not exist, add 
				selection_per_species[species] = [(random.choice(files_per_species[species]), "-")]
	# Write output
	with open(f"{args.output}/{args.filename}.tsv", "w") as f_out:
		for species in selection_per_species:
			for file, val in selection_per_species[species]:
				f_out.write(f"{species}\t{file}\t{val}\n")
	# Write output per experiment
	with open(f"{args.output}/genus_experiments/{args.filename}.tsv", "w") as f_out:
		for species in selection_per_species:
			if int(species) in species_per_taxon["genus"]:
				for file, val in selection_per_species[species]:
					f_out.write(f"{species}\t{file}\t{val}\n")
	with open(f"{args.output}/family_experiments/{args.filename}.tsv", "w") as f_out:
		for species in selection_per_species:
			if int(species) in species_per_taxon["family"]:
				for file, val in selection_per_species[species]:
					f_out.write(f"{species}\t{file}\t{val}\n")
	with open(f"{args.output}/order_experiments/{args.filename}.tsv", "w") as f_out:
		for species in selection_per_species:
			if int(species) in species_per_taxon["order"]:
				for file, val in selection_per_species[species]:
					f_out.write(f"{species}\t{file}\t{val}\n")

if __name__ == "__main__":
	main()