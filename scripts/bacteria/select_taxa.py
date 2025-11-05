from TaxTree import *
import numpy as np

def calculate_average_distance(species):
    # Calculate distance matrix
    dist_matrix = np.genfromtxt(f"root/genomes/{species}/converted_matrix.mat", delimiter="\t", skip_header=1, dtype=np.float64)
    # Retrieve upper triangular submatrix to calculate average distance -> use as criterion for selection
    dist_matrix = dist_matrix[np.triu_indices_from(dist_matrix, k=1)] #k=1 offset from diagonal
    return np.mean(dist_matrix)


def main():
    # Build taxonomic tree and populate it
    tax_tree = read_taxtree('root/taxdmp/nodes.dmp')
    tax_tree = Tree(tax_tree)

    print("Populating taxonomic tree")
    for node in tax_tree.nodes:
        tax_tree.find_lineage(tax_tree.nodes[node].taxid)
    taxid2name, taxid2syn = read_taxnames('root/taxdmp/names.dmp')
    print("Done populating taxonomic tree")

    # Determine number of genomes per species
    num_per_species = {}
    with open("root/reference_sets/all.tsv", "r") as f_in:
        next(f_in)
        for line in f_in:
            species, accession, _ = line.strip().split("\t")
            species = int(species)
            if species not in num_per_species:
                num_per_species[species] = 0
            num_per_species[species] += 1

    # Select species with at least 100 genomes and calculate average distance
    COUNT_THRESHOLD=100
    S = []
    for taxid in num_per_species:
        count = num_per_species[taxid]
        try:
            if count >= COUNT_THRESHOLD:
                species = tax_tree.nodes[taxid].taxonomy["species"]
                genus = tax_tree.nodes[taxid].taxonomy["genus"]
                family = tax_tree.nodes[taxid].taxonomy["family"]
                order = tax_tree.nodes[taxid].taxonomy["order"]
                avg_dist = calculate_average_distance(species)
                S.append((taxid, count, avg_dist, (species, taxid2name[species]), (genus, taxid2name[genus]), (family, taxid2name[family]), (order, taxid2name[order])))
        except:
            continue
    # Sort by average distance (descending)
    S.sort(key = lambda x: x[2], reverse=True)
    for s in S:
        print(s)
    print(len(S))

    # Now select 1 order -> 1 family -> 1 genus
    selection_per_order = {}
    selected_order = None
    selection_per_family = {}
    selected_family = None
    selection_per_genus = {}
    selected_genus = None
    selection_per_species = {}
    selected_species = None
    for s in S:
        if s[6][0] not in selection_per_order:
            selection_per_order[s[6][0]] = []
        if not selected_order:
            selection_per_order[s[6][0]].append(s)
            if len(selection_per_order[s[6][0]]) >= 5:
                selected_order = s[6]

        if s[5][0] not in selection_per_family:
            selection_per_family[s[5][0]] = []
        if not selected_family:
            selection_per_family[s[5][0]].append(s)
            if len(selection_per_family[s[5][0]]) >= 5:
                selected_family = s[5]

        if s[4][0] not in selection_per_genus:
            selection_per_genus[s[4][0]] = []
        if not selected_genus:
            selection_per_genus[s[4][0]].append(s)
            if len(selection_per_genus[s[4][0]]) >= 5:
                selected_genus = s[4]

    print(selected_order, selected_family, selected_genus)
    print("Order:", selected_order)
    for s in selection_per_order[selected_order[0]]:
        print("species:", s[0])
    #print(selection_per_order[selected_order[0]])
    print("Family:", selected_family)
    for s in selection_per_family[selected_family[0]]:
        print("species:", s[0])
    #print(selection_per_family[selected_family[0]])
    print("Genus:", selected_genus)
    for s in selection_per_genus[selected_genus[0]]:
        print("species:", s[0])
    #print(selection_per_genus[selected_genus[0]])

    # Create files with all species belonging to selected taxa
    with open("root/selected_genus_species.txt", "w") as f_out:
        for taxid in num_per_species:
            try:
                genus = tax_tree.nodes[taxid].taxonomy["genus"]
                if genus == selected_genus[0]:
                    f_out.write(f"{taxid}\n")
            except:
                continue
    with open("root/selected_family_species.txt", "w") as f_out:
        for taxid in num_per_species:
            try:
                family = tax_tree.nodes[taxid].taxonomy["family"]
                if family == selected_family[0]:
                    f_out.write(f"{taxid}\n")
            except:
                continue
    with open("root/selected_order_species.txt", "w") as f_out:
        for taxid in num_per_species:
            try:
                order = tax_tree.nodes[taxid].taxonomy["order"]
                if order == selected_order[0]:
                    f_out.write(f"{taxid}\n")
            except:
                continue

if __name__ == "__main__":
    main()