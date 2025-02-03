# Benchmarking the impact of reference genome selection on taxonomic profiling accuracy
This repository contains the scripts that were used to generate the results in the manuscript.

## General
The accession ids of all genomes used in this study are available on [Zenodo](https://doi.org/10.5281/zenodo.14727633). In our experiments, we ran every tool on a "per-species basis", meaning that we made a selection for every species independent of the selections for other species.

The "Bacteria" and "SARS-CoV-2" folders contain all scripts that were used to run the selection tools, as well as subfolders with the scripts for running the profiler(s). In most cases, the selection tools required a single multi-fasta file containing all of the genomes. When a genome had multiple segments or parts they were concatenated in the order in which they appeared in the corresponding fasta file downloaded from NCBI.


## Hierarchical clustering and centroid selections
Both the hierarchical clustering (based on MASH distances), and the centroid selections were obtained by running a simple Python script. This script loaded the distance matrix for a species, and afterwards ran the following lines of code for hierarchical clustering:

```python
Z = linkage(squareform(distance_matrix), method=args.linkage)
clusters = fcluster(Z, t=threshold, criterion="distance")
```

Then, for every `cluster` in `clusters`, we obtain the centroid as a representative for the corresponding cluster. This is also used for the centroid selection, except in that case we simply provided the full set of genomes rather than every cluster.
