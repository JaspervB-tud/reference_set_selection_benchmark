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
├── simulation_genomes
│   ├── genus_experiments
|   │   ├── simulation_species_1
│   |   |   ├── simulation_accession_1.fna
│   |   |   └── simulation_accession_2.fna
|   │   └── simulation_species_2
│   ├── family_experiments
│   └── order_experiments
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
│   |   |   ├── sample_1.fq <- forward reads
|   |   |   └── sample_2.fq <- reverse reads
│   |   ├── sample_2
│   |   |   ├── sample_1.fq
|   |   |   └── sample_2.fq
│   ├── family_experiments
|   └── order_experiments
├── estimations
│   ├── genus_experiments
│   |   ├── profiler_1
│   |   |	├── sample_1
│   |   |   |	├── method_1_predictions
|   |   |   |	├── method_2-threshold_1_predictions
|   |   |   |	└── method_2-threshold_2_predictions
│   |   ├── profiler_2
│   |   └── profiler_3
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

## Post-processing
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
**NOTE**: This should be run for every selection method (including `all`) and threshold in order to proceed to building the profiling indexes!

## Index building
With all the reference sets, we can now build the profiling indices. In our work we used:
- Kraken2 (v2.1.3) + Bracken (v1.0.0)
- BWA (v0.7.18) + DUDes (v0.10.0)
- Centrifuge (v1.0.4.2)

### Kraken2 + Bracken
For Kraken2+Bracken we have provided a simple bash script called `compile_kraken2-bracken.sh` that gathers all target genomes, and builds a combined index. The script can be called for every experiment as follows:
```bash
bash scripts/bacteria/compile_kraken2-bracken.sh ${method}_${threshold} root/reference_sets/${experiment} root/indexes/${experiment}/bracken root/genomes root/taxonomy root/reference_sets/nucl_gb.accession2taxid
```
This creates a combined Kraken2 and Bracken index for a given method-threshold combination and experiment, which is located in `/root/indexes/${experiment}/bracken/${method}_${threshold}`.

### Centrifuge
For Centrifuge we have a similar bash script called `compile_centrifuge.sh` which we call like:
```bash
bash scripts/bacteria/compile_centrifuge.sh ${method}_${threshold} root/reference_sets/${experiment} root/indexes/${experiment}/centrifuge root/genomes root/taxonomy root/reference_sets/nucl_gb.accession2taxid
```

### BWA + DUDes
Finally, for BWA with DUDes we use the script called `compile_bwa-dudes.sh` which can be called in a similar fashion:
```bash
bash scripts/bacteria/compile_bwa-dudes.sh ${method}_${threshold} root/reference_sets/${experiment} root/indexes/${experiment}/dudes root/genomes root/taxonomy root/reference_sets/nucl_gb.accession2taxid
```

## Profiling
The first step for obtaining the profiling results is to simulate reads. For our bacterial experiments, we simulated paired-end reads using ART v2016.06.05 from the top 10 (sorted according to completeness) genome assemblies that were "Scaffold" or "Chromosome" status (see Zenodo). Additionally, we store the files called `${experiment}_species.txt` and `${taxid}_accessions.txt` (see `files` folder) that contain the species and accessions used for simulating reads. Using these files, we simulated independent samples with approximately equal abundance (in terms of the number of reads) for every species using the following Python script using ART v2016.06.05:
```python
import subprocess
import os
import random
import shutil
from pathlib import Path
import math
import multiprocessing
from Bio import SeqIO

def generate_commands():
	"""
	This function will generate bash commands to perform the ART simulations using the parameters defined as constants below.
	"""
	NUM_SAMPLES 	= 10
	NUM_READS 		= 5_000_000 #approximate number of reads per sample
	FRAGMENT_SIZE 	= 270
	FRAGMENT_STD	= 20
	READ_LENGTH		= 150
	TECHONOLOGY		= "HS25"
	OUTPUT_PREFIX 	= f"root/samples"

	command_sets = [] #a single command set per sample -> this way an entire sample can be generated and processed in a single go (using 1 core)
	seed = 1

	for experiment_type in ["order_experiments", "family_experiments", "genus_experiments"]:
		print(f"Generating reads for {experiment_type}")
		# Retrieve list of species
		species = []
        with open(f"files/{experiment_type}/{experiment_type}_species.txt", "r") as f_in:
			next(f_in) #skip header
			for line in f_in:
				line = line.strip().split("\t")
				species.append(line[0])
		species.sort()
		# Retrieve list of accessions per species
		accessions_per_species = {}
		length_per_accession = {}
		for cur_species in species:
			accessions_per_species[cur_species] = []
            with open(f"files/{experiment_type}/{cur_species}_accessions.txt", "r") as f_in:
				for line in f_in:
					line = line.strip()
					accessions_per_species[cur_species].append(line)
					length_per_accession[line] = 0
                    with open(f"root/simulation_genomes/${experiment_type}/{cur_species}/{line}", "r") as fasta_file:
						for record in SeqIO.parse(fasta_file, "fasta"):
							length_per_accession[line] += len(record.seq)
			accessions_per_species[cur_species].sort()
		# Loop over samples
		for sample in range(1, NUM_SAMPLES+1):
			commands = []
			# Generate equal abundance samples
			num_reads_per_species = {cur_species: math.ceil(NUM_READS / 5) for cur_species in species}
			for cur_species in species:
				for cur_accession in accessions_per_species[cur_species]:
					num_reads_per_accession = math.ceil(num_reads_per_species[cur_species] * 0.1) #10 genomes for every species -> 1/10 abundance for all of them
                    output_path = f"{OUTPUT_PREFIX}/{experiment_type}/{sample}/{cur_species}" #this is output for reads per species
					output_path_object = Path(output_path)
					output_path_object.mkdir(parents=True, exist_ok=True)
					# We exclude the alignment files simply due to the space occupation of these files
					commands.append([
						"art_illumina",
						"-p",
                        "-i", f"root/simulation_genomes/{experiment_type}/{cur_species}/{cur_accession}",
						"-l", f"{READ_LENGTH}",
						"-f", f"{num_reads_per_accession * FRAGMENT_SIZE / length_per_accession[cur_accession]:.5f}", #this is to guarantee uniform coverage
						"-ss", "HS25",
						"-m", f"{FRAGMENT_SIZE}",
						"-s", f"{FRAGMENT_STD}",
						"-o", f"{output_path}/{cur_accession}_seed-{seed}_",
						"-rs", f"{seed}",
						"-na"
					])
					seed += 1
			command_sets.append(commands)
    return command_sets

def run_subprocess(command_set):
	"""
	The idea is that every command_set corresponds to a sample, and thus this function will:
		1) generate the required simulated reads
		2) combine the reads into a single file (since reads will be generated separately for every genome)
		3) remove the initial .fq files and only retain the combined files
	"""
	fwd_paths = []
	rev_paths = []
	reads_per_species = {}
	# Generate read sets
	for command in command_set:
		fwd_paths.append(f"{command[-4]}1.fq")
		rev_paths.append(f"{command[-4]}2.fq")
		subprocess.run(command)
		cur_species = f"{command[-4].split('/')[-2]}"
		cur_readcount = int(subprocess.run(["grep", "-c", "@", fwd_paths[-1]], capture_output=True, text=True).stdout.strip())
		if cur_species not in reads_per_species:
			reads_per_species[cur_species] = 0
		reads_per_species[cur_species] += cur_readcount
	# Combine forward reads and remove old files
	with open(f"{'/'.join(command[-4].split('/')[:-2])}/sample_1.fq", "w") as out_handle:
		for fq_file in fwd_paths:
			with open(fq_file, "r") as in_handle:
				SeqIO.write(SeqIO.parse(in_handle, "fastq"), out_handle, "fastq")
	# Combine reverse reads and remove old files
	with open(f"{'/'.join(command[-4].split('/')[:-2])}/sample_2.fq", "w") as out_handle:
		for fq_file in rev_paths:
			with open(fq_file, "r") as in_handle:
				SeqIO.write(SeqIO.parse(in_handle, "fastq"), out_handle, "fastq")
	# Remove folders with isolated reads per species
	prefix = "/".join(fwd_paths[0].split("/")[:-2])
	for species in reads_per_species:
		subprocess.run(["rm", "-r", f"{prefix}/{species}"])

if __name__ == "__main__":
    # Currently, this uses 16 cores to process all the ART commands
    NUM_CORES = 16
    commands = generate_commands()
    with multiprocessing.Pool(NUM_CORES) as pool:
        pool.map(run_subprocess, commands)
```
This script creates 10 mixed metagenomic samples per experiment with approximately 5,000,000 paired-end reads of 150bp (270bp fragments, 10bp stdev) and equal read-based abundance for all species, which are stored in `root/samples/${experiment}/sample_${sample}/sample_{1,2}.fq`. With the simulated samples, we can call the profilers to estimate the relative abundances of species using the constructed reference sets.

### Kraken2 + Bracken
Running Kraken2 and Bracken is done in two subsequent steps: first running Kraken2 and then running Bracken to get species-level abundance estimates. This was done for every sample and method-threshold combination as follows:
```bash
kraken2 --db root/indexes/${experiment}/bracken/${method}_${threshold} --threads 16 --report root/estimations/${experiment}/bracken/sample_${sample}/${method}_${threshold}.kreport --paired root/samples/${experiment}/sample_${sample}/sample_1.fq root/samples/${experiment}/sample_${sample}/sample_2.fq > root/estimations/${experiment}/bracken/sample_${sample}/${method}_${threshold}.kraken #generates Kraken2 output
bracken -d root/indexes/${experiment}/bracken/${method}_${threshold} -i root/estimations/${experiment}/bracken/sample_${sample}/${method}_${threshold}.kreport -o root/estimations/${experiment}/bracken/sample_${sample}/${method}_${threshold}.bracken -r 150 -l S #generates Bracken output
```
This produces several output files, of which the `${method}_${threshold}.bracken` file is the most relevant as it will contain the species-level abundances (although the `${method}_${threshold}.kreport` files can be used to infer the number of unclassified reads).

### Centrifuge
Centrifuge can be run per sample, reference set and experiment as follows:
```bash
centrifuge -p 16 -x root/indexes/${experiment}/centrifuge/${method}_${threshold}/index/${method}_${threshold} -1 root/samples/${experiment}/sample_${sample}/sample_1.fq -2 root/samples/${experiment}/sample_${sample}/sample_2.fq -S root/estimations/${experiment}/centrifuge/sample_${sample}/${method}_${threshold}.sam --report-file root/estimations/${experiment}/centrifuge/sample_${sample}/${method}_${threshold}.report
```
This produces an alignment-like file (`${method}_${threshold}.sam`) as well as a report file (`${method}_${threshold}.report`) which will contain the abundance estimates.

### BWA + DUDes
DUDes also requires multiple steps. First we use BWA-mem to align reads to the reference index:
```bash
bwa mem -t 16 -v 1 root/indexes/${experiment}/dudes/${method}_${threshold}/bwa_index/bwa root/samples/${experiment}/sample_${sample}/sample_1.fq root/samples/${experiment}/sample_${sample}/sample_2.fq > root/estimations/${experiment}/dudes/sample_${sample}/${method}_${threshold}.sam
```
After running BWA-mem we filter out unaligned reads and finally run DUDes on the filtered alignment file:
```bash
# Filter unaligned and calculate stats
samtools view -F 4 -h root/estimations/${experiment}/dudes/sample_${sample}/${method}_${threshold}.sam > root/estimations/${experiment}/dudes/sample_${sample}/${method}_${threshold}_filtered.sam
# Run DUDes
dudes -s root/estimations/${experiment}/dudes/sample_${sample}/${method}_${threshold}_filtered.sam -d root/indexes/${experiment}/dudes/${method}_${threshold}/dudes_index/dudes.npz -o root/estimations/${experiment}/dudes/sample_${sample}/${method}_${threshold}_dudes -l species
```
This results in a `${method}_${threshold}_dudes` file that contains the final abundance estimates.

## Analysis
To analyse the results obtained we will assume that the workflow was ran as described above. Additionally, for resource usage monitoring, we assume that commands were run with `/usr/bin/time` and that the output is stored somewhere.

### Reference set comparisons
In `scripts/bacteria/analysis_reference_sets.ipynb` we provide a Jupyter notebook that details the steps we took to obtain the results and figures (using matplotlib and seaborn) regarding the reference set comparisons. 

### Accuracy comparisons
In `scripts/bacteria/analysis_accuracy.ipynb` we provide a Jupyter notebook that details the steps we took to obtain the results and figures (using matplotlib and seaborn) regarding the accuracy calculations.

### Runtime comparisons
In `scripts/bacteria/analsys_runtime.ipynb` we provide a Jupyter notebook that details the steps we took to obtain the results and figures (using matplotlib and seaborn) regarding the runtime calculations. \
**NOTE**: for the bacterial results these strongly rely on the `/usr/bin/time` command!