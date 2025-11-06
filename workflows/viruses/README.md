# Workflow for SARS-CoV-2
This document details all of the steps required to reproduce the viral results in the manuscript.

## Data
The SARS-CoV-2 genomes and metadata were downloaded on June 12, 2022 from [GISAID](https://gisaid.org) using the accession ids or EPI-SET identifier that can be found on [Zenodo](https://zenodo.org/records/17524569). All other genomes/sequences that were used can also be found on Zenodo.

## Organization
In our experiments, we ran every tool on a "per-lineage" basis, leading to a selection for every SARS-CoV-2 lineage that is independent of the selections for other lineages. Throughout this workflow we assume the following general folder structure:
```
root
├── sequences.fasta <- downloaded from GISAID
├── metadata.tsv    <- downloaded from GISAID
├── genomes
│   ├── global
│   |   ├── lineage_1
|   |   │   ├── accession_1.fa
|   |   │   └── accession_2.fa
│   |   └── lineage_2
|   |       └── accession_3.fa
│   ├── country
│   └── state
├── selections
│   ├── global
│   |   ├── method_1
│   |   ├── method_2_threshold_1
│   |   ├── method_2_threshold_2
│   ├── country
|   └── state
├── reference_sets
│   ├── global
│   |   ├── method_1
│   |   ├── method_2-threshold_1
|   |   └── method_2-threshold_2
|   ├── country
|   └── state
├── indexes
│   ├── global
│   |   ├── method_1
│   |   ├── method_2-threshold_1
|   |   └── method_2-threshold_2
|   ├── country
|   └── state
├── samples
│   ├── sample_1
│   |   ├── wwsim_B.1.1.7_ab1_1.fq
|   |   ├── wwsim_B.1.1.7_ab1_2.fq
│   |   ├── wwsim_B.1.1.7_ab10_1.fq
|   |   └── wwsim_B.1.1.7_ab10_2.fq
├── estimations
│   ├── global
│   |   ├── sample_1
|   |   |   ├── method_1_ab1_predictions
|   |   |   ├── method_2-threshold_1_ab1_predictions
|   |   |   ├── method_2-threshold_2_ab1_predictions
|   |   |   ├── method_1_ab10_predictions
|   |   |   ├── method_2-threshold_1_ab10_predictions
|   |   |   └── method_2-threshold_2_ab10_predictions
|   ├── country
|   └── state
└── scripts
```

## Data collection and pre-processing
For our viral experiments, we downsampled to at most 1,000 genomes per lineage, filtering out genomes shorter than 25,000bp or with an N-content above 0.1%. We did this using the scripts provided in the [VLQ GitHub](https://github.com/baymlab/wastewater_analysis/tree/main). Throughout, we assume that the scripts in this (as well as those for other methods) are available!

The downsampling and filtering step can be applied by running the `preprocessing.sh` bash script that uses the VLQ scripts. In addition to the filtering criteria already mentioned, this script only retains genomes from January 1st, 2021 up until March 31st, 2021, and it additionally creates selections retaining only genomes from the state Connecticut, and the country USA. All resulting genomes can be found in `root/genomes/${location_type}` where they are stored in subfolders containing all the genomes for a particular lineage in individual fasta files.

## Genome selection
The next step in our experiments is to select reference genomes to be used by VLQ's profiling pipeline. In the viral experiments we used the following approaches:
- MASH (v2.3) + picking a centroid per lineage
- MASH (v2.3) + picking based on hierarchical clustering (single & complete linkage)
- MASH (v2.3) + GGRaSP (v1.3)
- MeShClust (Identity v1.2, MeShClust v2.0)
- Gclust (v1.0)
- VSEARCH (v2.28.1)
- VLQ

### MASH
We first generated MASH sketches for all genomes, which are used to estimate (dis)similarities. For all lineages with at least 2 genomes and for all experiments (global, country, state), we now run the following commands:
```bash
mash sketch -p 16 -s 5000 -S 123456 -k 31 -o root/genomes/${location_type}/${lineage}/mash_sketches.msh root/genomes/${location_type}/${lineage}/*.fa
mash triangle -p 16 root/genomes/${location_type}/${lineage}/mash_sketches.msh > root/genomes/${location_type}/${lineage}/mash_distances.dist
```
This creates both a MASH sketch for every genome (stored in the `mash_sketches.msh` files) using a sketch size of 5,000 and a k-mer size of 31, as well as a lower triangular distance matrix capturing the dissimilarity between all genomes within a lineage (stored in the `mash_distances.dist` files).

#### Distance matrix conversion
Since GGRaSP requires a full distance matrix instead of a lower triangular matrix, we used the `scripts/convert_matrix.py` Python script for converting the distance matrices as follows:
```bash
python -u scripts/convert_matrix.py --matrix root/genomes/${location_type}/${lineage}/mash_distances.dist --output root/genomes/${location_type}/${lineage}
```
This produces a `converted_matrix.mat` file with the full distance matrix that can be removed after running GGRaSP

### Centroid selection
To obtain the centroid genome for every lineage, we used the `scripts/run_centroid.py` Python script as follows:
```bash
python -u scripts/run_centroid.py --matrix root/genomes/${location_type}/${lineage}/mash_distances.dist --output root/selections/${location_type}/centroid/${lineage}
```
Afterwards, the selection can be found in `root/selections/${location_type}/centroid/${lineage}` for every location and lineage.

### Hierarchical clustering selection (dRep-derived approach)
Due to the high similarity between SARS-CoV-2 genomes of the same lineage, we used dynamic cut-off thresholds based on the 1st, 5th, 10th, 25th, 50th, 90th and 99th distance percentiles. We run the `scripts/virus/run_hierarchical.py` Python script as follows:
```bash
python -u scripts/virus/run_hierarchical.py --matrix root/genomes/${location_type}/${lineage}/mash_distances.dist --output root/selections/${location_type}/hierarchical/${lineage}
```
This produces selections that can be found in `root/selections/${location_type}/hierarchical/${lineage}/{single-linkage, complete-linkage}_${percentile}` for every value of percentile (1, 5, 10, 25, 50, 90, 99).

### GGRaSP selection
After running the `scripts/convert_matrix.py` script, GGRaSP can be run by calling:
```bash
Rscript scripts/run_ggrasp.R root/genomes/${location_type}/${lineage}/converted_matrix.mat root/selections/${location_type}/ggrasp/${lineage}
```
This uses the `scripts/run_ggrasp.R` R script which runs GGRaSP using the automatic threshold selection strategy, and the resulting selections can be found in `root/selections/${location_type}/ggrasp/${lineage}` for every location type and lineage. \
**NOTE**: In our experiments, GGRaSP failed to select reference genomes for several lineages, often for unclear reasons. 

### MeShClust, Gclust and VSEARCH selection
MeShClust, Gclust and VSEARCH all require input to be in the form of a single multi-fasta file containing all genomes to select from. For this we used the `scripts/concatenate_genomes.py` Python script as follows:
```bash
python -u scripts/concatenate_genomes.py --genomes root/genomes/${location_type}/${lineage} --output root/genomes/${location_type}/${lineage}
```
This creates a single multi-fasta file called `all_genomes.fasta` for every genome in a lineage. \
**NOTE**: this script assumes that all genomes of a lineage are stored in the same input directy, and that they end with the `.fa` extension.

#### MeShClust selection
As the maximum allowed similarity threshold in MeShClust is limited to 99%, we ran MeShClust with thresholds of 95% and 99% as follows:
```bash
meshclust -d root/genomes/${location_type}/${lineage}/all_genomes.fasta -o root/selections/${location_type}/meshclust/0.95/${lineage}/meshclust_0.95 -t 0.95 -c 64
meshclust -d root/genomes/${location_type}/${lineage}/all_genomes.fasta -o root/selections/${location_type}/meshclust/0.99/${lineage}/meshclust_0.99 -t 0.99 -c 64
```
This creates a file called `meshclust_${threshold}` in the `root/selections/${location_type}/meshclust/${threshold}/${lineage}` folder containing the clusters as well as assigned cluster representatives.

#### Gclust selection
To run Gclust, all genomes in the `all_genomes.fasta` file have to be sorted according to length (longest to shortest). For this, we used the perl script provided with Gclust:
```bash
perl sortgenome.pl --genomes-file root/genomes/${location_type}/${lineage}/all_genomes.fasta --sortedgenomes-file root/genomes/${location_type}/${lineage}/all_genomes_sorted.fasta
```
Afterwards, we ran Gclust with similarity thresholds of 95%, 99% and 99.9% with parameters as described in the [paper](https://academic.oup.com/gpb/article/17/5/496/7229746):
```bash
gclust -threads 64 -minlen 41 -chunk 400 -ext 1 -sparse 4 -nuc -loadall -rebuild -memiden 95 root/genomes/${location_type}/${lineage}/all_genomes_sorted.fasta > root/selections/${location_type}/gclust/0.95/${lineage}/gclust_0.95
gclust -threads 64 -minlen 41 -chunk 400 -ext 1 -sparse 4 -nuc -loadall -rebuild -memiden 99 root/genomes/${location_type}/${lineage}/all_genomes_sorted.fasta > root/selections/${location_type}/gclust/0.99/${lineage}/gclust_0.99
gclust -threads 64 -minlen 41 -chunk 400 -ext 1 -sparse 4 -nuc -loadall -rebuild -memiden 99.9 root/genomes/${location_type}/${lineage}/all_genomes_sorted.fasta > root/selections/${location_type}/gclust/0.999/${lineage}/gclust_0.999
```
**NOTE**: For the SARS-CoV-2 experiments we omitted the `-both` flag.

#### VSEARCH selection
We ran VSEARCH with the `--cluster_fast` option using the CD-HIT definition for sequence similarity, again using similarity thresholds of 95%, 99% and 99.9%:
```bash
vsearch --cluster_fast root/genomes/${location_type}/${lineage}/all_genomes.fasta --centroids root/selections/${location_type}/vsearch/0.95/${lineage}/vsearch_0.95 --id 0.95 --iddef 0 --qmask none --threads 64
vsearch --cluster_fast root/genomes/${location_type}/${lineage}/all_genomes.fasta --centroids root/selections/${location_type}/vsearch/0.99/${lineage}/vsearch_0.99 --id 0.99 --iddef 0 --qmask none --threads 64
vsearch --cluster_fast root/genomes/${location_type}/${lineage}/all_genomes.fasta --centroids root/selections/${location_type}/vsearch/0.999/${lineage}/vsearch_0.999 --id 0.999 --iddef 0 --qmask none --threads 64
```
In contrast to other selection tools, this outputs a fasta file containing the full selected genomes, rather than identifiers of sequences which we will account for later.

### VLQ selection
For selecting SARS-CoV-2 genomes we also use the selection methodology provided in the VLQ pipeline. However, instead of the default implementation, we use an adaptation of the `call_variants.sh` bash script which operates on individual single lineages. The adapted script can be found in `scripts/virus/run_vlq.sh` and was run for every lineage in every experiment as follows:
```bash
bash scripts/virus/run_vlq.sh ${lineage} ${ref_genome} root/genomes/${location_type} root/selections/${location_type}/vlq
```
This script uses minimap2 to call variants against a reference genome, which in our case is the MN908947.3 reference. After running for all lineages, the `root/selections/${location_type}/vlq` folder will contain a `${lineage}_merged.frq` file, a `${lineage}_merged.pi` and `${lineage}_merged.vcf.gz` file for every lineage.

To perform VLQ's selection steps, make sure to run the previous command for ALL lineages before proceeding, making sure that the `root/selections/${location_type}/vlq` folder contains all genomes (from every lineage). The selection can then by obtained by running the `select_samples.py` Python script from the VLQ pipeline as follows:
```bash
python -u select_samples.py -m root/metadata.tsv -f root/sequences.fasta -o root/selections/${location_type}/vlq --vcf root/selections/${location_type}/vlq/*_merged.vcf.gz --freq root/selections/${location_type}/vlq/*_merged.fq --max_per_lineage 1000
```
This results in a `sequences.fasta` and a `metadata.tsv` file in the `root/selections/${location_type}/vlq` folder containing all of the genome sequences and metaddata for the selected sequences (for all lineages) respectively.

## Post-processing
After all the tools are run, the selected reference genomes should all be located in the `root/selections` folder for all experiments and lineages. From here we perform two steps: first is generating an index file that stores all genomes (as their filename) and their corresponding lineage. This can be done by running the `generate_all_selection.py` Python script available in the `scripts` folder:
```bash
python -u scripts/generate_all_selection.py --genomes root/genomes/${location_type} --output root/reference_sets/${location_type}
```
This creates a tab-delimited file called `all.tsv` in the `root/reference_sets/${location_type}` folder containing all of the filenames, their corresponding lineage and an identifier that can be used to find lineages for which a tool was unable to select (not relevant for "all" selection) structured as follows:
```
LINEAGE GENOME_FILENAME    +/-
```

Afterwards, we run the `generate_selection_files.py` Python script in `scripts/virus` to create similar files for the other selections:
```bash
python -u scripts/virus/generate_selection_files.py --genomes root/genomes/${location_type} --selection root/selections/${location_type}/${method}/${threshold} --filename ${method}_${threshold} --output root/reference_sets/${location_type}
```
Running this will create a tab-delimited file called `${method}_${threshold}.tsv`, similar to `all.tsv`. Additionally, if a method was unable to produce a selection for a lineage, this script will randomly select a genome from that lineage in order to represent it in the reference set.\
**NOTE**: This should be run for every experiment, selection method (including `all`) and threshold in order to proceed to building the kallisto indexes!

## Index building
With all the reference sets, we can now build kallisto indexes. In our experiments we followed the VLQ pipeline with kallisto v0.44.0 using the `compile_kallisto.sh` script in `scripts/virus`:
```bash
bash scripts/virus/compile_kallisto.sh ${method}_${threshold} root/reference_sets/${location_type} root/indexes/${location_type} root/genomes/${location_type}
```
This script first collects all reference genomes in a selection and temporarily stores them in a single multi-fasta file called `all_genomes.fasta`. Afterwards, it builds the corresponding kallisto index in `root/indexes/${location_type}/${method}_${threshold}.idx`, and removes the temporary file. \
**NOTE**: This, again, should be run for all experiments and reference sets!

## Profiling
The first step for obtaining the profiling results is to simulate reads. For our experiments, we simulated reads from SARS-CoV-2 genomes sampled on 30th of April, 2021 from the state Connecticut as well as a corresponding B.1.1.7 genome (*EPI_ISL_2159878*). We did this by simulating independent samples with B.1.1.7 at (relative) abundances of 1%, 10%, 20%, ..., 100%. For the simulations we ran the following script that uses ART v2016.06.05:
```bash
mkdir -p root/samples
metadata="root/metadata.tsv"
sequences="root/sequences.fasta"
state="Connecticut"

mkdir -p root/samples
for seed in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do \
    python -u scripts/virus/create_benchmarks_with_seed.py --voc_perc 1,10,20,30,40,50,60,70,80,90,100 -m $metadata -fr $sequences -fv root/B.1.1.7_sequence.fasta -o root/samples/$seed --total_cov 100 -s "North America / USA / ${state}" -d 2021-04-30 --seed $seed

    for file in root/samples/$seed/*; do
        if ! [[ "$file"  =~ ^(root/samples/$seed/wwsim_|root/samples/$seed/metadata.tsv|root/samples/$seed/sequences.fasta).*$ ]]; 
        then 
            echo "$file is removed"
            rm $file
        fi
    done
done
```
This script uses the `create_benchmarks_with_seed.py` Python script in `scripts/virus` to create 20 samples for every abundance level of B.1.1.7, and stores the corresponding paired-end reads (~9,000 reads, 150bp readlength, 250bp fragment length, 10bp stdev) in `root/samples`, where every sample is its own folder containing the samples for every abundance of B.1.1.7.

With the simulated metagenomic samples, we call kallisto and the `output_abundances.py` Python script from VLQ as follows:
```bash
kallisto quant -b 0 -i root/indexes/${location_type}/${method}_${threshold}.idx -o root/estimations/${location_type}/${sample}/${method}_${threshold} root/samples/${sample}/wwsim_B.1.1.7_ab${abundance}_1.fastq root/samples/${sample}/wwsim_B.1.1.7_ab${abundance}_2.fastq
python -u output_abundances.py -m 0.1 -o root/estimations/${location_type}/${sample}/${method}_${threshold}_predictions.tsv --metadata root/metadata.tsv root/estimations/${location_type}/${sample}/${method}_${threshold}/abundance.tsv
```
This will store the direct output of kallisto in `root/estimations/${location_type}/${sample}/${method}_${threshold}/abundance.tsv` and the filtered (at 0.1% abundance) lineage abundance estimations in `root/estimations/${location_type}/${sample}/${method}_${threshold}_predictions.tsv`.

## Analysis
To analyse the results obtained we will assume that the workflow was ran as described above. Additionally for most of the resource usage monitoring, we assume that commands were run with `/usr/bin/time/` and that the output is stored somewhere.

### Reference set comparisons
In `scripts/virus/analysis_reference_sets.ipynb` we provide a Jupyter notebook that details the steps we took to obtain the results and figures (using matplotlib and seaborn) regarding the reference set comparisons.