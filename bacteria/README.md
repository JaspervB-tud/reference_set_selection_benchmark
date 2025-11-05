# Workflow for Bacteria
This document details all of the steps required to reproduce the bacterial results in the manuscript.

## Data
All of the bacterial genomes were downloaded on October 9, 2024 from [NCBI](https://www.ncbi.nlm.nih.gov), and the taxonomy data was downloaded on September 24, 2024 from NCBI. The accession ids for all genomic sequences as well as simulated reads used throughout this study are available on [Zenodo](https://zenodo.org/records/17524569).

## Organization
In our experiments, we ran every tool on a "per-species" basis, leading to a selection for every bacterial species that is independent of the selections for other species. Throughout this workflow we assume the following general folder structure:
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
│   ├── method_2_threshold_1
│   ├── method_2_threshold_2
│   ├── genus_experiments
│   |   ├── method_1
│   |   ├── method_2_threshold_1
│   |   ├── method_2_threshold_2
│   ├── family_experiments
│   └── order_experiments
├── reference_sets
│   ├── method_1
│   ├── method_2_threshold_1
│   ├── method_2_threshold_2
│   ├── genus_experiments
│   |   ├── method_1
│   |   ├── method_2-threshold_1
|   |   └── method_2-threshold_2
|   ├── family_experiments
|   └── order_experiments
├── indexes
│   ├── genus_experiments
│   |   ├── profiler_1
│   |   |   ├── method_1
│   |   |   ├── method_2-threshold_1
|   |   |   └── method_2-threshold_2
|   |   ├── profiler_2
|   |   └── profiler_3
│   ├── family_experiments
|   └── order_experiments
├── samples
│   ├── genus_experiments
│   |   ├── sample_1
│   |   |   ├── sample_1_1.fq
|   |   |   └── sample_1_2.fq
│   ├── family_experiments
|   └── order_experiments
├── estimations
│   ├── genus_experiments
│   |   ├── sample_1
│   |   |   ├── method_1_predictions
|   |   |   ├── method_2-threshold_1_predictions
|   |   |   └── method_2-threshold_2_predictions
│   ├── family_experiments
|   └── order_experiments
├── taxdmp
│   ├── names.dmp
│   └── nodes.dmp
└── scripts
```

## Data collection
For our bacterial experiments, we used all available "Complete Genome" status from NCBI Assembly (now NCBI Genomes) downloaded on October 9th, 2024. In this, we only included species that had at least one such assembly for our analysis. Throughout this, we assume that the individual genomes are located in the `root/genomes/${species}` folders for all corresponding species.

## Genome selection
The first step in our experiments is to perform reference genome selection. For the bacterial experiments we used the following approaches:
- MASH (v2.3) + picking a centroid per lineage
- MASH (v2.3) + picking based on hierarchical clustering (single & complete linkage)
- MASH (v2.3) + GGRaSP (v1.3)
- MeShClust (Identity v1.2, MeShClust v2.0)
- Gclust (v1.0)
Unfortunately, due to scaling problems, we were unable to run VSEARCH for these experiments.

### MASH
We first generated MASH sketches for all genomes, which are used to estimate (dis)similarities. For all species with at least 2 genomes and for all experiments (genus, family, order), we now run the following commands:
```bash
mash sketch -p 16 -s 1000 -S 123456 -k 21 -o root/genomes/${species}/mash_sketches.msh root/genomes/${species}/*.fa
mash triangle -p 16 root/genomes/${species}/mash_sketches.msh > root/genomes/${experiment}/mash_distances.dist
```
This creates both a MASH sketch for every genome (stored in the `mash_sketches.msh` files) using a sketch size of 1,000 and a k-mer size of 21, as well as a lower triangular distance matrix capturing the dissimilarity between all genomes within a lineage (stored in the `mash_distances.dist` files).

#### Distance matrix conversion
Since GGRaSP requires a full distance matrix instead of a lower triangular matrix, we used the `scripts/convert_matrix.py` Python script for converting the distance matrices as follows:
```bash
python -u scripts/convert_matrix.py --matrix root/genomes/${species}/mash_distances.dist --output root/genomes/${species}
```
This produces a `converted_matrix.mat` file with the full distance matrix that can be removed after running GGRaSP

### Centroid selection
To obtain the centroid genome for every species, we used the `scripts/run_centroid.py` Python script as follows:
```bash
python -u scripts/run_centroid.py --matrix root/genomes/${species}/mash_distances.dist --output root/selections/centroid/${species}
```
Afterwards, the selection can be found in `root/selections/centroid/${species}` for every and bacterial species.

### Hierarchical clustering selection (dRep-derived approach)
Unlike the viral experiments, we use fixed similarity thresholds of 95%, 97% and 99% for the bacterial experiments. To obtain the hierarchical clustering-based selections, we run the `scripts/bacteria/run_hierarchical.py` Python script as follows:
```bash
python -u scripts/bacteria/run_hierarchical.py --matrix root/genomes/${species}/mash_distances.dist --threshold 0.01 --output root/selections/hierarchical/${species} # 0.01 -> 99% similarity
python -u scripts/bacteria/run_hierarchical.py --matrix root/genomes/${species}/mash_distances.dist --threshold 0.03 --output root/selections/hierarchical/${species} # 0.03 -> 97% similarity
python -u scripts/bacteria/run_hierarchical.py --matrix root/genomes/${species}/mash_distances.dist --threshold 0.05 --output root/selections/hierarchical/${species} # 0.05 -> 95% similarity
```
This produces selections that can be found in `root/selections/hierarchical/${species}/{single-linkage, complete-linkage}_${threshold}`.

### GGRaSP selection
After running the `scripts/convert_matrix.py` script, GGRaSP can be run by calling:
```bash
Rscript scripts/run_ggrasp.R root/genomes/${species}/converted_matrix.mat root/selections/ggrasp/${species}
```
This uses the `scripts/run_ggrasp.R` R script which runs GGRaSP using the automatic threshold selection strategy, and the resulting selections can be found in `root/selections/ggrasp/${species}` for every location type and species. \
**NOTE**: In our experiments, GGRaSP failed to select reference genomes for several species, often for unclear reasons. 

### MeShClust, Gclust and VSEARCH selection
MeShClust, Gclust and VSEARCH all require input to be in the form of a single multi-fasta file containing all genomes to select from. For this we used the `scripts/concatenate_genomes.py` Python script as follows:
```bash
python -u scripts/concatenate_genomes.py --genomes root/genomes/${species} --output root/genomes/${species}
```
This concatenates genomes consisting of multiple compartments (e.g. chromosomes) and creates a single multi-fasta file called `all_genomes.fasta` for every genome in a species. \
**NOTE**: this script assumes that all genomes of a species are stored in the same input directy, and that they end with the `.fa` extension.

#### MeShClust selection
We ran MeShClust with thresholds of 95%, 97% and 99% as follows:
```bash
meshclust -d root/genomes/${species}/all_genomes.fasta -o root/selections/meshclust/0.95/${species}/meshclust_0.95 -t 0.95 -c 64
meshclust -d root/genomes/${species}/all_genomes.fasta -o root/selections/meshclust/0.97/${species}/meshclust_0.97 -t 0.97 -c 64
meshclust -d root/genomes/${species}/all_genomes.fasta -o root/selections/meshclust/0.99/${species}/meshclust_0.99 -t 0.99 -c 64
```
This creates a file called `meshclust_${threshold}` in the `root/selections/meshclust/${threshold}/${species}` folder containing the clusters as well as assigned cluster representatives.

#### Gclust selection
To run Gclust, all genomes in the `all_genomes.fasta` file have to be sorted according to length (longest to shortest). For this, we used the perl script provided with Gclust:
```bash
perl sortgenome.pl --genomes-file root/genomes/${species}/all_genomes.fasta --sortedgenomes-file root/genomes/${species}/all_genomes_sorted.fasta
```
Afterwards, we ran Gclust with similarity thresholds of 95%, 97% and 99% with parameters as described in the [paper](https://academic.oup.com/gpb/article/17/5/496/7229746):
```bash
gclust -threads 64 -minlen 41 -both -chunk 400 -ext 1 -sparse 4 -nuc -loadall -rebuild -memiden 95 root/genomes/${species}/all_genomes_sorted.fasta > root/selections/gclust/0.95/${species}/gclust_0.95
gclust -threads 64 -minlen 41 -both -chunk 400 -ext 1 -sparse 4 -nuc -loadall -rebuild -memiden 97 root/genomes/${species}/all_genomes_sorted.fasta > root/selections/gclust/0.97/${species}/gclust_0.97
gclust -threads 64 -minlen 41 -both -chunk 400 -ext 1 -sparse 4 -nuc -loadall -rebuild -memiden 99 root/genomes/${species}/all_genomes_sorted.fasta > root/selections/gclust/0.99/${species}/gclust_0.99
```

# Post-processing
After running the tools, the selected reference genomes are all located in the `root/selection` folder for all species. From here, we proceed with three steps. First, we generate index files that store all genomes (as their filenames) and their corresponding species. This can be done by running the `generate_all_selection.py` Python script in the `scripts` folder:
```bash
python scripts/generate_all_selection.py --genomes root/genomes --output root/reference_sets --a2t
```
This creates two files. The first file is a tab-delimited file called `all.tsv` in the `root/reference_sets` folder containing all of the filenames, their corresponding species and an identifier that can be used to find species for which a tool was unable to select (not relevant for "all" selection) structured as follows:
```
SPECIES GENOME_FILENAME +/-
```
The second file is a `nucl_gb.accession2taxid` file stored in `root/reference_sets` which mimics the accession2taxid file in NCBI's taxonomy, but only retains accessions relevant to our usecase.

In the second step we determine the species that will be used in our experiments since using all bacterial species was computationally infeasible. For this we select an order, family and genus by considering species that have at least 100 genomes. Afterwards, remaining species were sorted according to their average intra-species distance, from highest to lowest and a taxon (order, family, genus) was chosen if it was the first to have 5 species. The Python script to do this can be found in `scripts/bacteria/select_taxa.py` and uses the `TaxTree.py` helper script as well as the NCBI taxonomy files in `root/taxdmp`. The output of this script further consists of all species that were part of the chosen taxa (per taxon) which is stored in `root/selected_${taxon}_genomes.tsv`.

In the final step we create files similar to the `all.tsv` for every selection, across all experiments (as well as for all available species). We do this by calling the `generate_selection_files.py` Python script in `scripts/bacteria`:
```bash
python -u root/scripts/bacteria/generate_selection_files.py --genomes root/genomes --selection root/selections/${method}/${threshold} --filename ${method}_${threshold} --output root/reference_sets
```
Similar to before, this will create tab-delimited files called `${method}_${threshold}.tsv` in `root/reference_sets` as well as the subfolders thereof for the experiments in order to only include the selections for considered species. \
**NOTE**: This should be run for every selection method and threshold in order to proceed to building the profiling indexes!

## Index building
With all the reference sets,w e can now build the profiling indices. In our work we used:
- Kraken2 (v2.1.3) + Bracken (v1.0.0)
- BWA (v0.7.18) + DUDes (v0.10.0)
- Centrifuge (v1.0.4.2)

### Kraken2 + Bracken
