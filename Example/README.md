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