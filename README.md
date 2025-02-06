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
In our experiments, we ran every tool on a "per-taxon basis", meaning that we made a selection for every taxon (e.g. species or SARS-CoV-2 lineage) independent of the selections for other species. When using the scripts provided in this github we assume the following general folder structure (species would be replaced by lineage in case of SARS-CoV-2):

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
|   │   ├── species_1
|   │   └── species_2
│   ├── method_2
|   │   ├── threshold_1
|   │   |   ├── species_1
|   │   |   └── species_2
|   │   ├── threshold_2
|   │   |   ├── species_1
|   │   |   └── species_2
└── scripts
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
We first generate MASH sketches for all genomes, which are then used to estimate (dis)similarity between them. From hereon we assume that the user is in the `root` directory as shown in the folder structure in the Organization section. Note that we run these commands FOR EVERY SPECIES/LINEAGE, but only show the command for a single species/lineage here. If a species/lineage only has a single available assembly/genome, then the selection is omitted, and the single available genome is always used.
```bash
# In the bacterial setting
mash sketch -s 1000 -S 123456 -k 21 -o genomes/species_X/mash_sketches.msh genomes/species_X/*.fa
mash triangle genomes/species_X/mash_sketches.msh > genomes/species_X/mash_distances.dist

# In the SARS-CoV-2 setting
mash sketch -s 5000 -S 123456 -k 31 -o genomes/lineage_X/mash_sketches.msh genomes/lineage_X/*.fa
mash triangle genomes/lineage_X/mash_sketches.msh > genomes/lineage_X/mash_distances.dist
```
This results in a MASH sketch for every genome (located in the `mash_sketches.msh` files), as well as a lower triangular distance matrix between all genomes of a species (`mash_distances.dist`).

#### Distance matrix conversion
As GGRaSP requires the input to consist of a full matrix, rather than the lower triangular matrix outputted by `mash triangle`, we use the `convert_matrix.py` Python script for converting the distances matrices.
```bash
python scripts/convert_matrix.py --matrix genomes/species_X/mash_distances.dist --output genomes/species_X
```
This produces a `converted_matrix.mat` file with the full distance matrix, that can be removed after running GGRaSP.

### Centroid selection
To obtain the centroid genome for every species, we used the `run_centroid.py` script provided in the scripts folder, storing the resulting selection in the `selections` folder:
```bash
python scripts/run_centroid.py --matrix genomes/species_X/mash_distances.dist --output selections/centroid/species_X
```
After running, the selection is found in `selections/centroid/species_X/centroid`.

### Hierarchical clustering selection (dRep-based approach)
#### Bacteria
To obtain the hierarchical clustering selection for the bacterial genomes we used the `run_hierarchical.py` script with thresholds set to 0.01, 0.03, 0.05 which correspond to similarities of 99%, 97% and 95% respectively:
```bash
python scripts/run_hierarchical.py --matrix genomes/species_X/mash_distances.dist --threshold 0.01 --output selections/hierarchical/0.99/species_X
python scripts/run_hierarchical.py --matrix genomes/species_X/mash_distances.dist --threshold 0.03 --output selections/hierarchical/0.97/species_X
python scripts/run_hierarchical.py --matrix genomes/species_X/mash_distances.dist --threshold 0.05 --output selections/hierarchical/0.95/species_X
```
After running, this produces 2 selections for every threshold: one with single-linkage clustering, and one with complete-linkage clustering which can be found in `selections/hierarchical/THRESHOLD/species_X/single-linkage_(1-THRESHOLD)` and `selections/hierarchical/THRESHOLD/species_X/complete-linkage_(1-THRESHOLD)` respectively.

#### SARS-CoV-2
Due to the high similarity between SARS-CoV-2 genomes, instead of using fixed thresholds we instead calculate percentile-based distance thresholds for 1, 5, 10, 25, 50, 90 and 99 percentile distances using the `run_hierarchical_sc2.py` script:
```bash
python scripts/run_hierarchical_sc2.py --matrix genomes/lineage_X/mash_distances.dist --output selections/hierarchical/lineage_X
```
The resulting selections can then be found in `selections/hierarchical/lineage_X/single-linkage_PERCENTILE` and `selections/hierarchical/lineage_X/complete-linkage_PERCENTILE` for every percentile.

### GGRaSP selection
As stated before, GGRaSP requires the full distance matrices rather than (lower) triangular matrices, so make sure to run the `convert_matrix.py` script before running the following:
```bash
Rscript scripts/run_ggrasp.R genomes/species_X/converted_matrix.mat selections/ggrasp/species_X
```
The resulting selection can then be found in `selections/ggrasp/species_X/ggrasp`. **IMPORTANT**: in our experiments, GGRaSP failed to select genomes several times, and for various different reasons. In the scenario where GGRaSP (or any other tool for that matter) failed to select, we randomly picked a single representative genome for every species/lineage where it failed.

### MeShClust, Gclust and VSEARCH selection
MeShClust, Gclust and VSEARCH require a single multi-fasta input file containing all genomes from which to select. As several of the genomes in our bacterial experiments had more than one entry in the corresponding fasta file (e.g. multiple chromosomes), these were first concatenated to a single entry in the order in which they appear in the fasta file. Then all (concatenated) genomes were added to a single fasta file using the `concatenate_genomes.py` script:
```bash
python scripts/concatenate_genomes.py --genomes genomes/species_X --output genomes/species_X
```
Running this results in a fasta file called `all_genomes.fasta`. Note that the script assumes that all non-concatenated genomes are stored in the same input directory ending in the `*.fa` extension.

#### MeShClust selection
We run MeShClust directly on the generated `all_genomes.fasta` files to obtain selections. For the bacterial experiments, we ran MeShClust with thresholds of 95%, 97% and 99%, whereas for the SARS-CoV-2 experiments we ran MeShClust with thresholds of 95% and 99% (since 99% is the maximum allowed threshold):
```bash
meshclust -d genomes/species_X/all_genomes.fasta -o selections/meshclust/THRESHOLD/species_X/meshclust_THRESHOLD -t 0.99 #for a threshold of 99% substitute THRESHOLD with 0.99 etc.
```

#### Gclust selection
Similar to MeShClust, Gclust requires the `all_genomes.fasta` file to perform selection, however, Gclust also requires the genomes in the file to be sorted by length (longest to shortest). They have provided a perl script for this which can be run as follows:
```bash
perl sortgenome.pl --genomes-file genomes/species_X/all_genomes.fasta --sortedgenomes-file genomes/species_X/all_genomes_sorted.fasta
```
Running as such will generate the sorted `all_genomes_sorted.fasta` file and stores it in the same directory as the individual genome files for the species/lineage of interest. Afterwards, we run Gclust with similarity threshold of 95%, 97% and 99% for the bacterial selections, and 95%, 99% and 99.9% for the SARS-CoV-2 selections. Additionally, in the SARS-CoV-2 selections we omitted the `-both` flag.
```bash
gclust -minlen 41 -both -chunk 400 -ext 1 -sparse 4 -nuc -loadall -rebuild -memiden THRESHOLD genomes/species_X/all_genomes_sorted.fasta > selections/gclust/THRESHOLD/species_X/gclust_0.THRESHOLD #for a threshold of 99%, substitute THRESHOLD with 99 etc.
```

#### VSEARCH selection
We tried running VSEARCH for the bacterial setting but this turned out to take too long. Therefore, we only ran VSEARCH for the SARS-CoV-2 setting:
```bash
vsearch --cluster_fast genomes/species_X/all_genomes.fasta --centroids selections/vsearch/THRESHOLD/species_X/vsearch_THRESHOLD --id THRESHOLD --iddef 0 --qmask none
```
In contrast to all other tools so far, VSEARCH outputs a fasta file containing the full sequences of all the selected genomes, rather than the ids of the selected sequences.

