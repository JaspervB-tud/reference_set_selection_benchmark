import json

class Node():
    def __init__(self, taxid, parent, rank):
        self.taxid = taxid
        self.parent = parent
        self.rank = rank
        self.taxonomy = {'strain' : None, 'subspecies' : None, 'species' : None, 'genus' : None, 'family' : None, 'order' : None, 'class' : None, 'phylum' : None, 'superkingdom' : None}
        self.full_path = None
        
    def __eq__(self, other):
        if type(other) == Node:
            return self.taxid == other.taxid
        return False
    
    def set_taxonomy(self, taxonomy):
        for entry in taxonomy:
            if entry[1] in self.taxonomy:
                self.taxonomy[entry[1]] = entry[0]

class Tree():
    def __init__(self, nodes):
        self.nodes = nodes
        
    def find_lineage(self, taxid):
        taxid = int(taxid)
        cur_node = self.nodes[taxid]
        path = [(cur_node.taxid, cur_node.rank)]
        while cur_node.taxid != 1:
            cur_node = self.nodes[cur_node.parent]
            path.append((cur_node.taxid, cur_node.rank))
        self.nodes[taxid].set_taxonomy(path)
        return path
    
    def find_lineages(self, taxids):
        paths = {}
        total_not_found = 0
        for taxid in taxids:
            try:
                paths[taxid] = self.find_lineage(taxid)
            except:
                paths[taxid] = 'Not found'
                total_not_found += 1
        print('A total of', total_not_found, 'taxids were not found')
        return paths
    
    def find_lineages_from_accessionid(self, accessions2taxids):
        paths = {}
        not_found = []
        total_not_found = 0
        for accession in accessions2taxids:
            cur_taxid = accessions2taxids[accession]
            try:
                paths[accession] = self.find_lineage(cur_taxid)
            except:
                paths[accession] = (None, None)
                not_found.append(accession)
                total_not_found += 1
        print('A total of', total_not_found, 'accessions did not have a taxid that is in taxdump')
        return paths, not_found
    
    def find_ids_per_species(self, taxids):
        not_found = set()
        _ = self.find_lineages(taxids)
        ids_per_species = {}
        for taxid in taxids:
            try:
                cur_node = self.nodes[int(taxid)]
                cur_species = cur_node.taxonomy['species']
                if cur_species:
                    if cur_species not in ids_per_species:
                        ids_per_species[cur_species] = []
                    ids_per_species[cur_species].append(taxid)
            except:
                not_found.add(taxid)
        return ids_per_species, not_found

def read_taxtree(path):
    taxtree = {}
    with open(path, 'r') as f:
        for line in f:
            line = line.split('|')
            cur_taxid = int(line[0].strip())
            cur_parent = int(line[1].strip())
            cur_rank = line[2].strip()
            cur_node = Node(cur_taxid, cur_parent, cur_rank)
            taxtree[cur_taxid] = cur_node
    return taxtree

def read_taxnames(path):
    names = {}
    synonyms = {}
    with open(path, 'r') as f:
        for line in f:
            if 'scientific name' in line:
                line = line.split('|')
                cur_taxid = int(line[0].strip())
                names[cur_taxid] = line[1].strip()
                if cur_taxid not in synonyms:
                    synonyms[cur_taxid] = []
                synonyms[cur_taxid].append(line[1].strip())
            elif "synonym" in line:
                line = line.split("|")
                cur_taxid = int(line[0].strip())
                if cur_taxid not in synonyms:
                    synonyms[cur_taxid] = []
                synonyms[cur_taxid].append(line[1].strip())
    return names, synonyms
    
def read_accession2taxid(path):
    accession2taxid = {}
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                line = line.split('\t')
                accession2taxid[line[0].strip()] = int(line[1])
    return accession2taxid

def read_accession2path(path):
    accession2path = {}
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                line = line.split('\t')
                accession2path[line[0].strip()] = [eval(l) for l in line[1:]]
    return accession2path

def read_taxid2accession(path):
    taxid2accession = {}
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                line = line.split('\t')
                taxid2accession[int(line[0])] = line[1].strip()
    return taxid2accession

def read_taxid2name(path):
    taxid2name = {}
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                line = line.split('*|_|*')
                taxid2name[int(line[0])] = line[1].strip()
    return taxid2name

def main():
    print("Hello from TaxTree.py")
    
if __name__ == '__main__':
    main()
    
    
    
    
