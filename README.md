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
In our experiments, we ran every tool on a "per-taxon basis", meaning that we made a selection for every taxon (e.g. species or SARS-CoV-2 lineage) independent of the selections for other species. When using the scripts provided in this github we assume the following general folder structure for the genomes:

```
Root
├── Species 1
│   ├── Accession 1
│   └── Accession 2
└── Species 2
    ├── Accession 3
```

The "Bacteria" and "SARS-CoV-2" folders contain all scripts that were used to run the selection tools, as well as subfolders with the scripts for running the profiler(s). In most cases, the selection tools required a single multi-fasta file containing all of the genomes. When a genome had multiple segments or parts they were concatenated in the order in which they appeared in the corresponding fasta file downloaded from NCBI.


## Hierarchical clustering and centroid selections
Both the hierarchical clustering (based on MASH distances), and the centroid selections were obtained by running a simple Python script. For hierarchical clustering, this script loaded the MASH distance matrix for a species, and afterwards used `scipy`'s clustering functionality for hierarchical clustering. Then, for every cluster, the centroid was selected as its representative.

The centroid selection followed a similar selection strategy, but instead of clustering first, immediately returns the centroid of all genomes for a species/lineage.
