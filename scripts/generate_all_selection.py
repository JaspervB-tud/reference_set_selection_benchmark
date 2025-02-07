from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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



if __name__ == "__main__":
	main()

"""
import os
import glob
import gzip
import tqdm
import math
import argparse
import time
from multiprocessing import Pool
from Bio import SeqIO

def read_accessions_in_genome(data):
	'''
	Intention for this is to be multiprocessing -> data consists of [filename, kingdom, taxid].
	This function parses the file corresponding to "filename" in the kingdom "kingdom" and in species "taxid", and
	returns a tuple containing the accession ids in the file.
	'''
	filename, kingdom, taxid = data
	sequences = []
	with gzip.open(f"/tudelft.net/staff-umbrella/refsetbenchmark/species/{kingdom}/{taxid}/{filename}.fna.gz", "rt") as f_in:
		records = SeqIO.parse(f_in, "fasta")
		for record in records:
			sequences.append(record)

	return sequences

def read_reference_set(method):
	file_path = f"/tudelft.net/staff-umbrella/refsetbenchmark/reference_sets/_genomes_included_per_method/{method}_corrected.txt"

	accessions_per_species_per_kingdom = {"2157": {}, "2": {}, "10239": {}}
	tot_accessions = 0
	with open(file_path, "r") as f_in:
		for line in f_in:
			kingdom, species, accession, _ = line.strip().split("\t")
			if species not in accessions_per_species_per_kingdom[kingdom]:
				accessions_per_species_per_kingdom[kingdom][species] = []
			accessions_per_species_per_kingdom[kingdom][species].append(accession)
			tot_accessions += 1

	print(f"Reference set {method} has {tot_accessions} total accessions (files)")
	return accessions_per_species_per_kingdom


def main():
	batch_size = 1000
	num_cores = 8

	parser = argparse.ArgumentParser()
	parser.add_argument("--method", type=str, required=True)
	args = parser.parse_args()
	
	methods = [
		"centroid",
		"ggrasp",
		"single-linkage-0.95",
		"single-linkage-0.97",
		"single-linkage-0.99",
		"complete-linkage-0.95",
		"complete-linkage-0.97",
		"complete-linkage-0.99",
		"meshclust-0.95",
		"meshclust-0.97",
		"meshclust-0.99",
		"gclust-0.95",
		"gclust-0.97",
		"gclust-0.99",
		"all"
	]

	for method in [args.method]:
		print(f"Processing method \"{method}\"")
		accessions_per_species_per_kingdom = read_reference_set(method)
		all_accessions = []

		kingdom = "2"
		for species in accessions_per_species_per_kingdom[kingdom]:
			for accession in accessions_per_species_per_kingdom[kingdom][species]:
				all_accessions.append((accession, kingdom, species)) #filename, kingdom, taxid

		# Start writing out to "combined" fasta file that includes all sequences
		print(f"Starting to write sequence to /tudelft.net/staff-umbrella/refsetbenchmark/reference_sets/_genomes_included_per_method/{method}_all_sequences_bacteria.fna.gz")
		print(f"Working in batches of size {batch_size} ({math.ceil(len(all_accessions) / batch_size)} batches)")
		
		if method == "all":
			output_location = f"/tudelft.net/staff-umbrella/refsetopt/{method}_all_sequences_bacteria.fna.gz"
		else:
			output_location = f"/tudelft.net/staff-umbrella/refsetbenchmark/reference_sets/_genomes_included_per_method/{method}_all_sequences_bacteria.fna.gz"
		start_time = time.time()
		with gzip.open(output_location, "wt") as f_out:
			if num_cores > 1:
				j = 0
				for i in range(0, len(all_accessions), batch_size):
					batch = all_accessions[i:i+batch_size]
					batch_results = []
					with Pool(processes=num_cores) as pool:
						batch_results = pool.map(read_accessions_in_genome, batch)
					for sequences in batch_results:
						for sequence in sequences:
							SeqIO.write(sequence, f_out, "fasta")
					j += 1
					print(f"{j} batch(es) done ({(time.time() - start_time)/(j)}s per batch)")
			else:
				count = 0
				batch_count = 0
				for accession in all_accessions:
					cur_sequences = read_accessions_in_genome(accession)
					for sequence in cur_sequences:
						SeqIO.write(sequence, f_out, "fasta")
					count += 1
					if count == batch_size:
						print(f"{batch_count+1} batch(es) done")
						batch_count += 1
						count = 0	

if __name__ == "__main__":
	main()
"""