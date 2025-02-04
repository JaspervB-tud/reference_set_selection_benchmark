# Benchmarking the impact of reference genome selection on taxonomic profiling accuracy
This repository contains the scripts that were used to generate the results in the manuscript.

## General
The accession ids of all genomes used in this study are available on [Zenodo](https://doi.org/10.5281/zenodo.14727633). In our experiments, we ran every tool on a "per-taxon basis", meaning that we made a selection for every taxon (e.g. species or SARS-CoV-2 lineage) independent of the selections for other species.

The "Bacteria" and "SARS-CoV-2" folders contain all scripts that were used to run the selection tools, as well as subfolders with the scripts for running the profiler(s). In most cases, the selection tools required a single multi-fasta file containing all of the genomes. When a genome had multiple segments or parts they were concatenated in the order in which they appeared in the corresponding fasta file downloaded from NCBI.


## Hierarchical clustering and centroid selections
Both the hierarchical clustering (based on MASH distances), and the centroid selections were obtained by running a simple Python script. For hierarchical clustering, this script loaded the MASH distance matrix for a species, and afterwards used `scipy`'s clustering functionality for hierarchical clustering. Then, for every cluster, the centroid was selected as its representative.

The centroid selection followed a similar selection strategy, but instead of clustering first, immediately returns the centroid of all genomes for a species/lineage.
