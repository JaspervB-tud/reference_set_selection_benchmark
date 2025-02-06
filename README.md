# Benchmarking the impact of reference genome selection on taxonomic profiling accuracy (this is currently JUST the abstract)
Motivation: Over the past decades, genome databases have expanded exponentially, often incorporating
highly similar genomes at the same taxonomic level. This redundancy can hinder taxonomic classification,
leading to difficulties distinguishing between closely related sequences and increasing computational
demands. While some novel taxonomic classification tools address this redundancy by filtering the input
genome database, there is limited work exploring the impact of different sequence dereplication methods
across taxonomic classification tools.

Results: We assess the effect of genome selection approaches on taxonomic classification using both bacterial
and viral datasets. Our results demonstrate that careful selection of reference genomes can improve classification
accuracy. We also show that using prior knowledge about a metagenomic sample, such as sampling
location, can significantly improve classification accuracy. Finally, we find that using a redundancy-filtered
genome database generally reduces the computational resources required, with minimal loss in classification
accuracy.

## Data
The accession ids of all genomes used for our manuscript are available on [Zenodo](https://doi.org/10.5281/zenodo.14727633). The bacterial assemblies can be downloaded from NCBI (e.g. through NCBI's entrez), and the SARS-CoV-2 genomes as well as corresponding metadata can be downloaded from GISAID using the provided accession ids.

## Organization
In our experiments, we ran every tool on a "per-taxon basis", meaning that we made a selection for every taxon (e.g. species or SARS-CoV-2 lineage) independent of the selections for other species. When using the scripts provided in this github we assume the following general folder structure (Species would be replaced by Lineage in case of SARS-CoV-2):

```
root
├── genomes
│   ├── species_1
|   │   ├── accession_1.fa
|   │   └── accession_2.fa
│   └── species_2
|   │   └── accession_3.fa
├── selections
│   ├── method_1
|   │   ├── threshold_1
|   │   |   ├── species_1
|   │   |   └── species_2
│   ├── method_2
|   │   ├── threshold_1
|   │   |   ├── species_1
|   │   |   └── species_2
|   │   ├── threshold_2
|   │   |   ├── species_1
|   │   |   └── species_2
```

## Genome selection
The first step in our experiments is to select reference genomes, to be used by taxonomic profilers. In our experiments we used the following approaches:
- MASH (v2.3) + picking a centroid per cluster
- MASH (v2.3) + picking based on hierarchical clustering (single and complete linkage)
- MASH (v2.3) + GGRaSP (v1.3)
- MeShClust (Identity v1.2, MeShClust v2.0)
- Gclust (v1.0)

In addition, for the SARS-CoV-2 experiments we also used VSEARCH (v2.28.1) and VLQ's selection methodology.

### MASH
We first generate MASH sketches for all genomes, which are then used to estimate (dis)similarity between them. From hereon we assume that the user is in the `Genomes` directory as shown in the folder structure in the Organization section. Note that we run these commands FOR EVERY SPECIES/LINEAGE, in order to ensure that there is a selection for every species/lineage. If a species/lineage only has a single available assembly/genome, then the selection is omitted, and the single available genome is always used.
```bash
mash sketch -s 1000 -S 123456 -k 21 -o genomes_per_species/0/mash_sketches.msh genomes_per_species/0/*.fa
mash triangle genomes_per_species/0/mash_sketches.msh > genomes_per_species/0/mash_distances.dist
```