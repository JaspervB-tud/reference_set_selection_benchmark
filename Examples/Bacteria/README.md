The goal of selection is to find a subset of genomes that collectively represent all of the available genomes for a species. We have done this with the following approaches:

- Centroid selection based on MASH distances
- Hierarchical clustering based on MASH distances (as in dRep)
- GGRaSP based on MASH distances
- Gclust
- MeShClust

After selecting reference genomes, we assess the performance of taxonomic profilers based on the selected reference genomes. In our work we have done this for:

- Kraken2 + Bracken
- Centrifuge
- BWA + DUDes

In what follows we will provide the steps taken to generate the results in our manuscript.

## Genome selection
### MASH-based selection
We first generate MASH sketches, which are used to estimate the (dis)similarity between genomes. For the example provided here, all genomes are decompressed and stored in the "genomes_per_species" folder, where every subfolder is labeled according to the (fake) taxid of the species and contains the genomes of the species.

Assuming that the user is in this folder, we would now run the following two lines for every species (we only show the lines for running it for species "0"):

```bash
mash sketch -s 1000 -S 123456 -k 21 -p 2 -o genomes_per_species/0/mash_sketches.msh genomes_per_species/0/*.fa
mash triangle -p 2 genomes_per_species/0/mash_sketches.msh > genomes_per_species/0/mash_distances.dist
```

After doing so, we obtain a MASH sketch for every genome (located in the `mash_sketches.msh` files), as well as a lower triangular distance matrix between all genomes of a species.

#### Hierarchical clustering selection
To obtain the hierarchical clustering selection, we can now run `run_hierarchical.py`, storing the results in the "selections" folder. In the manuscript, we ran the script with both single-linkage and complete-linkage linkages, and thresholds equal to 0.01, 0.03 and 0.05 which correspond to similarities of 99%, 97% and 95% respectively. The following snippet shows how this was done (for species 0) with all three thresholds:

```bash
python ../../Bacteria/run_hierarchical.py --matrix genomes_per_species/0/mash_distances.dist --threshold 0.01 --output selections/0
python ../../Bacteria/run_hierarchical.py --matrix genomes_per_species/0/mash_distances.dist --threshold 0.03 --output selections/0
python ../../Bacteria/run_hierarchical.py --matrix genomes_per_species/0/mash_distances.dist --threshold 0.05 --output selections/0
```