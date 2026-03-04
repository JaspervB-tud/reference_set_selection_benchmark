from TaxTree import *
import argparse
import copy
import os
from Bio import SeqIO
import random

def find_genomes(parent_folder_path, tax_tree):
    '''
    This function will parse "parent_folder_path", and searches every subfolder to find all files ending in .fa, storing them in a dictionary with species/strains as keys.
    NOTE: to comply with the proposed folder structure, this function assumes that the subfolders of the parent folder are formatted
    as "species_X" or "strain_X" where X refers to the taxid of the species/strain.
    Input:
        - parent_folder_path (str): path to the folder that contains subfolders for all species/strains
        - tax_tree (Tree): taxonomic tree containing all taxa of interest, used for parsing the lineage of each taxon and retrieving the corresponding taxid
    Output:
        - files_per_taxon (dict[str] -> list[str]): dict with all genome file names for all species/strains (keys)
        - accessions_per_taxon (dict[str] -> list[str]): dict with all sequence ids for all genomes for all species/strains (keys)
        - genomes (dict[str] -> list[dict]): dict with all genome information for all species/strains (keys). The value is a list of dicts, where each dict corresponds to a genome and has the following structure:
    '''
    files_per_taxon = {}
    accessions_per_taxon = {}
    genomes = {}
    for folder in os.listdir(parent_folder_path):
        subfolder_path = os.path.join(parent_folder_path, folder)
        genomes[folder] = []
        if os.path.isdir(subfolder_path): #only consider folders
            cur_taxon = int(folder.split("_")[-1])
            cur_strain = tax_tree.nodes[cur_taxon].taxonomy["strain"]
            cur_genus = tax_tree.nodes[cur_taxon].taxonomy["genus"]
            cur_family = tax_tree.nodes[cur_taxon].taxonomy["family"]
            cur_order = tax_tree.nodes[cur_taxon].taxonomy["order"]
            files_per_taxon[cur_taxon] = []
            accessions_per_taxon[cur_taxon] = []
            for file in os.listdir(subfolder_path): #iterate over .fa files and retrieve all sequence ids for accessionid2taxid
                if file.endswith(".fa"):
                    files_per_taxon[cur_taxon].append(file)
                    accessions_per_taxon[cur_taxon] += [record.id for record in SeqIO.parse(os.path.join(subfolder_path, file), "fasta")]
                    genome_info = {}
                    genome_info["accession"] = file.split(".fa")[0] #strip any version and .fa extension
                    records = SeqIO.parse(os.path.join(subfolder_path, file), "fasta")
                    genome_info["records"] = []
                    for record in records:
                        record_copy = copy.deepcopy(record)
                        cur_id = record_copy.id
                        record_copy.id = cur_id
                        record_copy.name = cur_id
                        record_copy.description = ""
                        genome_info["records"].append(record_copy)
                    genome_info["strain"] = cur_strain
                    genome_info["genus"] = cur_genus
                    genome_info["family"] = cur_family
                    genome_info["order"] = cur_order
                    genomes[folder].append(genome_info)

    return files_per_taxon, accessions_per_taxon, genomes

def read_medoid(path, filename):
    selection = []
    with open(f"{path}/{filename}", "r") as f_in:
        for line in f_in:
            selection.append(line.strip())

    return selection

def read_hierarchical(path, filename):
    selection = []
    with open(f"{path}/{filename}", "r") as f_in:
        for line in f_in:
            selection.append(line.strip())

    return selection

def read_ggrasp(path, filename):
    selection = []
    with open(f"{path}/{filename}", "r") as f_in:
        for line in f_in:
            selection.append(line.strip())

    return selection

def read_gclust(path, filename):
    selection = []
    with open(f"{path}/{filename}", "r") as f_in:
        for line in f_in:
            line = line.strip().split()
            if line[-1] == "*":
                selection.append(line[2][1:-3]) #strips the leading ">" and trailing "..."

    return selection

def read_meshclust(path, filename):
    selection = []
    with open(f"{path}/{filename}", "r") as f_in:
        for line in f_in:
            line = line.strip()
            if len(line) > 0:
                line = line.split("\t")
                if line[-1] == "C":
                    selection.append(line[1][1:]) #strips the leading ">"

    return selection

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes_folder", type=str, required=True)
    parser.add_argument("--selections_folder", type=str, required=True)
    parser.add_argument("--taxdmp_folder", type=str) #this should be a taxdmp from NCBI
    parser.add_argument("--taxdmp_output_folder", type=str) #this should be the folder to store the custom taxdmp files
    args = parser.parse_args()

    # Taxa of interest
    taxa_of_interest = {
        "strain": 562, # Escherichia coli (species)
        "genus": 1301, # Streptococcus
        "family": 543, # Enterobacteriaceae
        "order": 91347 # Enterobacterales
    }

    # Gather taxonomic information
    tax_tree = read_taxtree(f"{args.taxdmp_folder}/nodes.dmp")
    tax_tree = Tree(tax_tree)

    # Populate taxonomic tree
    for node in tax_tree.nodes:
        tax_tree.find_lineage(tax_tree.nodes[node].taxid)
    taxid2name, _ = read_taxnames(f"{args.taxdmp_folder}/names.dmp")

    # Fetch all genome files and sequence ids for all species/strains
    files_per_taxon, accessions_per_taxon, genomes = find_genomes(args.genomes_folder, tax_tree)

    # Process strain data
    E_COLI_TAXID = 562 #strain experiments were based on E. coli
    lineage = tax_tree.find_lineage(E_COLI_TAXID)
    taxa = [l[0] for l in lineage]
    strains = [taxon for taxon in files_per_taxon]

    # Create custom taxonomy for strains
    taxdmp_lines = [] #this will contain all the nodes.dmp lines for the custom taxonomy
    with open(f"{args.taxdmp_folder}/nodes.dmp", "r") as f_in:
        for line in f_in:
            line = line.strip().split("\t|\t")
            taxid = int(line[0])
            if taxid in taxa:
                taxdmp_lines.append(line)
            elif taxid in strains:
                line[1] = "562" #set parent taxid to E. coli so that all strains are a direct child of E. coli
                taxdmp_lines.append("\t|\t".join(line) + "\t|\n")
    # Write custom taxonomy to file
    os.makedirs(args.taxdmp_output_folder, exist_ok=True)
    with open(f"{args.taxdmp_output_folder}/nodes.dmp", "w") as f_out:
        for line in taxdmp_lines:
            f_out.write("\t|\t".join(line) + "\t|\n")

    # Create custom taxonomy names file
    taxnames_lines = [] #this will contain all the names.dmp lines for the custom taxonomy
    with open(f"{args.taxdmp_folder}/names.dmp", "r") as f_in:
        for line in f_in:
            line = line.strip().split("\t|\t")
            taxid = int(line[0])
            if taxid in taxa or taxid in strains:
                taxnames_lines.append(line)
    # Write custom taxonomy names to file
    with open(f"{args.taxdmp_output_folder}/names.dmp", "w") as f_out:
        for line in taxnames_lines:
            f_out.write("\t|\t".join(line) + "\t|\n")

    # Create custom accession2taxid file
    accession2taxid_lines = [] #this will contain all the lines for the custom accession2taxid file
    accession2taxid_centrifuge_lines = [] #this will contain all the lines for the custom accession2taxid file in centrifuge format
    accession2taxid = {}
    for taxon in files_per_taxon:
        for genome_info in genomes[taxon]:
            try: #set taxon to strain for E. coli, original taxon otherwise
                if int(taxon) == 562 and genome_info["strain"] is not None:
                    taxon = genome_info["strain"]
            except:
                pass
            accession2taxid[genome_info["accession"]] = taxon
            for record in genome_info["records"]:
                accession2taxid_lines.append(f"{'.'.join(record.id.split('.')[:-1])}\t{record.id}\t{taxon}\t0")
                accession2taxid_centrifuge_lines.append(f"{record.id}\t{taxon}")
    # Write custom accession2taxid files
    with open(f"{args.taxdmp_output_folder}/full.accession2taxid", "w") as f_out:
        for line in accession2taxid_lines:
            f_out.write(line + "\n")
    with open(f"{args.taxdmp_output_folder}/partial.accession2taxid", "w") as f_out:
        for line in accession2taxid_centrifuge_lines:
            f_out.write(line + "\n")

    # Process reference sets
    for experiment in ["strain", "genus", "family", "order"]:
        os.makedirs(f"{args.selections_folder}/{experiment}_experiments", exist_ok=True)
        taxon_of_interest = taxa_of_interest[experiment]

        # All reference set
        with open(f"{args.selections_folder}/{experiment}_experiments/all.tsv", "w") as f_out:
            for taxon in files_per_taxon:
                if tax_tree.nodes[int(taxon)].taxonomy[experiment] == taxon_of_interest: #only consider species/strains that are part of the experiment
                    for file in files_per_taxon[taxon]:
                        f_out.write(f"{taxon}\t{file}\t+\n")

        # Medoid reference set
        with open(f"{args.selections_folder}/{experiment}_experiments/medoid.tsv", "w") as f_out:
            for taxon in files_per_taxon:
                if tax_tree.nodes[int(taxon)].taxonomy[experiment] == taxon_of_interest: #only consider species/strains that are part of the experiment
                    # Try-except block to handle cases where medoid selection was not successful for a taxon, in which case we randomly select a genome from the taxon
                    try:
                        selection = read_medoid(f"{args.selections_folder}/medoid/{taxon}", "medoid")
                        success = "+"
                    except:
                        selection = [random.choice(files_per_taxon[taxon])]
                        success = "-"
                    for selected in selection:
                        f_out.write(f"{taxon}\t{selected}\t{success}\n")

        # Hierarchical clustering reference sets
        if experiment == "strain":
            thresholds = ["0.95", "0.99", "0.999"]
        else:
            thresholds = ["0.95", "0.97", "0.99"]
        for linkage in ["single", "complete"]:
            for threshold in thresholds:
                with open(f"{args.selections_folder}/{experiment}_experiments/hierarchical_{linkage}_{threshold}.tsv", "w") as f_out:
                    for taxon in files_per_taxon:
                        if tax_tree.nodes[int(taxon)].taxonomy[experiment] == taxon_of_interest: #only consider species/strains that are part of the experiment
                            # Try-except block to handle cases where hierarchical clustering selection was not successful for a taxon, in which case we randomly select a genome from the taxon
                            try:
                                selection = read_hierarchical(f"{args.selections_folder}/hierarchical/{taxon}", f"{linkage}-linkage_{threshold}")
                                success = "+"
                            except:
                                selection = [random.choice(files_per_taxon[taxon])]
                                success = "-"
                            for selected in selection:
                                f_out.write(f"{taxon}\t{selected}\t{success}\n")

        # GGRaSP reference set
        with open(f"{args.selections_folder}/{experiment}_experiments/ggrasp.tsv", "w") as f_out:
            for taxon in files_per_taxon:
                if tax_tree.nodes[int(taxon)].taxonomy[experiment] == taxon_of_interest: #only consider species/strains that are part of the experiment
                    # Try-except block to handle cases where GGRaSP selection was not successful for a taxon, in which case we randomly select a genome from the taxon
                    try:
                        selection = read_ggrasp(f"{args.selections_folder}/ggrasp/{taxon}", "ggrasp")
                        success = "+"
                    except:
                        selection = [random.choice(files_per_taxon[taxon])]
                        success = "-"
                    for selected in selection:
                        f_out.write(f"{taxon}\t{selected}\t{success}\n")

        # MeShClust reference sets
        if experiment == "strain":
            thresholds = ["0.95", "0.99"]
        else:
            thresholds = ["0.95", "0.97", "0.99"]
        for threshold in thresholds:
            with open(f"{args.selections_folder}/{experiment}_experiments/meshclust_{threshold}.tsv", "w") as f_out:
                for taxon in files_per_taxon:
                    if tax_tree.nodes[int(taxon)].taxonomy[experiment] == taxon_of_interest: #only consider species/strains that are part of the experiment
                        # Try-except block to handle cases where MeShClust selection was not successful for a taxon, in which case we randomly select a genome from the taxon
                        try:
                            selection = read_meshclust(f"{args.selections_folder}/meshclust/{threshold}/{taxon}", f"meshclust_{threshold}")
                            success = "+"
                        except:
                            selection = [random.choice(files_per_taxon[taxon])]
                            success = "-"
                        for selected in selection:
                            f_out.write(f"{taxon}\t{selected}\t{success}\n")

        # Gclust reference sets
        if experiment == "strain":
            thresholds = ["0.95", "0.99", "0.999"]
        else:
            thresholds = ["0.95", "0.97", "0.99"]
        for threshold in thresholds:
            with open(f"{args.selections_folder}/{experiment}_experiments/gclust_{threshold}.tsv", "w") as f_out:
                for taxon in files_per_taxon:
                    if tax_tree.nodes[int(taxon)].taxonomy[experiment] == taxon_of_interest: #only consider species/strains that are part of the experiment
                        # Try-except block to handle cases where Gclust selection was not successful for a taxon, in which case we randomly select a genome from the taxon
                        try:
                            selection = read_gclust(f"{args.selections_folder}/gclust/{threshold}/{taxon}", f"gclust_{threshold}.clusters")
                            success = "+"
                        except:
                            selection = [random.choice(files_per_taxon[taxon])]
                            success = "-"
                        for selected in selection:
                            f_out.write(f"{taxon}\t{selected}\t{success}\n")                

if __name__ == "__main__":
    main()