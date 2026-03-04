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
│   ├── strain_experiments
|   │   ├── simulation_strain_1
│   |   |   ├── simulation_accession_1.fna
│   |   |   └── simulation_accession_2.fna
|   │   └── simulation_strain_2
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
│   ├── 500000
│	│   ├── genus_experiments
│	│   |   ├── sample_1
│	│   |   |   ├── sample_1.fq <- forward reads
│	|   |   |   └── sample_2.fq <- reverse reads
│	│   |   ├── sample_2
│	│   |   |   ├── sample_1.fq
│	|   |   |   └── sample_2.fq
│	│   ├── family_experiments
│	|   └── order_experiments
|   └── 5000000
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
For our bacterial experiments, we used all available "Complete Genome" status from NCBI Assembly (now NCBI Genomes) downloaded on October 9th, 2024. In this, we only included species and strains that had at least one such assembly for our analysis. Throughout this, we assume that the individual genomes are located in the `root/genomes/${species,strain}` folders for all corresponding species and *Escherichia coli* strains. \
**NOTE**: although we assume that all individual genomes are stored in `root/genomes/${taxid}`, both for strain taxids, as well as species taxids.

## Genome selection
The first step in our experiments is to perform reference genome selection. For the bacterial experiments we used the following approaches:
- MASH (v2.3) + picking the medoid per species/strain
- MASH (v2.3) + picking based on hierarchical clustering (single & complete linkage)
- MASH (v2.3) + GGRaSP (v1.3)
- MeShClust (Identity v1.2, MeShClust v2.0)
- Gclust (v1.0)
Unfortunately, due to scaling problems, we were unable to run VSEARCH for these experiments.

### MASH
We first generated MASH sketches for all genomes, which are used to estimate (dis)similarities. For all species with at least 2 genomes, and for all species-level experiments (genus-based, family-based, order-based), we ran the following commands:
```bash
mash sketch -p ${threads} -s 1000 -S 123456 -k 21 -o root/genomes/${species}/mash_sketches.msh root/genomes/${species}/*.fa
mash triangle -p ${threads} root/genomes/${species}/mash_sketches.msh > root/genomes/${species}/mash_distances.dist
```
For the strain-level experiments, we ran the same commands, except we replace `${species}` with `{strain}`, and we use a sketch size of 5,000 and k-mer size of 31. Doing this creates both a MASH sketch for every genome (stored in the `mash_sketches.msh` files), as well as a lower triangular distance matrix capturing the (dis)similarity between all genomes within a species/strain (stored in the `mash_distances.dist` files).

#### Distance matrix conversion
Since GGRaSP requires a full distance matrix instead of a lower triangular matrix, we used the `scripts/convert_matrix.py` Python script for converting the distance matrices as follows:
```bash
python -u scripts/convert_matrix.py --matrix root/genomes/${species,strain}/mash_distances.dist --output root/genomes/${species,strain}
```
This produces a `converted_matrix.mat` file with the full distance matrix that can be removed after running GGRaSP.

### Medoid selection
To obtain the medoid genome for every species, we used the `scripts/run_medoid.py` Python script as follows:
```bash
python -u scripts/run_medoid.py --matrix root/genomes/${species,strain}/mash_distances.dist --output root/selections/medoid/${species,strain}
```
Afterwards, the selection can be found in `root/selections/medoid/${species,strain}` for every bacterial species and *E. coli* strain.

### Hierarchical clustering selection (dRep-derived approach)
Unlike the viral experiments, we use fixed similarity thresholds in the bacterial experiments. In the species-level experiments, we used fixed thresholds of 95%, 97% and 99% similarity, and for the strain-level experiments we used fixed threholds of 95%, 99% and 99.9% similarity. To obtain the hierarchical clustering-based selections, we run the `scripts/bacteria/run_hierarchical.py` Python script as follows:
```bash
python -u scripts/bacteria/run_hierarchical.py --matrix root/genomes/${species,strain}/mash_distances.dist --threshold ${threshold} --output root/selections/hierarchical/${species,strain} # for a threshold of 99%, use 0.99
```
After running, selections that can be found in `root/selections/hierarchical/${species,strain}/{single-linkage, complete-linkage}_${threshold}`.

### GGRaSP selection
After running the `scripts/convert_matrix.py` script, GGRaSP can be run by calling:
```bash
Rscript scripts/run_ggrasp.R root/genomes/${species,strain}/converted_matrix.mat root/selections/ggrasp/${species,strain}
```
This uses the `scripts/run_ggrasp.R` R script which runs GGRaSP using the automatic threshold selection strategy, and the resulting selections can be found in `root/selections/ggrasp/${species,strain}` for every species/strain. \
**NOTE**: In our experiments, GGRaSP failed to select reference genomes for several species and strains, often for unclear reasons. When this happened, we represented the missing species/strains with a single, randomly chosen genome from the corresponding species/strain.

### MeShClust, Gclust and VSEARCH selection
MeShClust, Gclust and VSEARCH all require input to be in the form of a single multi-fasta file containing all genomes to select from. For this we used the `scripts/concatenate_genomes.py` Python script as follows:
```bash
python -u scripts/concatenate_genomes.py --genomes root/genomes/${species,strain} --output root/genomes/${species,strain}
```
This concatenates genomes consisting of multiple compartments (e.g. chromosomes) and creates a single multi-fasta file called `all_genomes.fasta` for every genome in a species/strain. \
**NOTE**: this script assumes that all genomes of a species/strain are stored in the same input directy, and that they end with the `.fa` extension.

#### MeShClust selection
We ran MeShClust with thresholds of 95%, 97% and 99% for species-level experiments, and 95% and 99% for strain-level experiments (99.9% was excluded as MeShClust does not allow thresholds above 99%) as follows:
```bash
meshclust -d root/genomes/${species,strain}/all_genomes.fasta -o root/selections/meshclust/${threshold}/${species,strain}/meshclust_${threshold} -c ${threads} -t ${threshold} # for a similarity threshold of 99%, use 0.99
```
This creates a file called `meshclust_${threshold}` in the `root/selections/meshclust/${threshold}/${species,strain}` folder containing the clusters as well as assigned cluster representatives.

#### Gclust selection
To run Gclust, all genomes in the `all_genomes.fasta` file have to be sorted according to length (longest to shortest). For this, we used the perl script provided with Gclust:
```bash
perl sortgenome.pl --genomes-file root/genomes/${species,strain}/all_genomes.fasta --sortedgenomes-file root/genomes/${species,strain}/all_genomes_sorted.fasta
```
Afterwards, we ran Gclust with similarity thresholds of 95%, 97% and 99% for the species-level experiments, and 95%, 99% and 99.9% for the strain-level experiments, using the parameters described in the [Gclust paper](https://academic.oup.com/gpb/article/17/5/496/7229746):
```bash
gclust -threads ${threads} -minlen 41 -both -threads ${threads} -chunk 400 -ext 1 -sparse 4 -nuc -loadall -rebuild -memiden ${threshold} root/genomes/${species,strain}/all_genomes_sorted.fasta > root/selections/gclust/0.${threshold}/${species,strain}/gclust_0.${threshold}.clusters # for a similarity threshold of 99%, use 99
```
This creates a file called `gclust_0.${threshold}.clusters` in `root/selections/gclust/0.${threshold}/${species,strain}`, similar to the output created by MeShClust.

## Post-processing
After running the tools, the selected reference genomes are all located in the `root/selection` folder for all species/strain taxids. The next step is to create selection manifests for every experiment, which document which genomes have been selected by every method-threshold combination. This can be done by running the `prepare_files.py` Python script in the `root/scripts/bacteria` folder as follows:
```bash
python scripts/bacteria/prepare_files.py --genomes_folder root/genomes --selections_folder root/selections --taxdmp_folder root/taxdmp --taxdmp_output_folder root/taxdmp_restricted
```
This creates multiple files. First is a tab-separated file per method-threshold combination, with each line represeting a selected genome and structured as:
```
TAXID GENOME_FILENAME +/-
```
with the first entry the taxid, second entry a the filename containing the genome filename, and the third entry identifying if a selection was succesful, or if it was randomly selected afterwards.

In addition to the selection files, this script also creates surrogate taxdmp files analogue to NCBI's taxdmp files, but only containing the species and strains considered here. Additionally, it has a nodes.dmp file with a custom taxonomy for the *E. coli* strains, such that every strain is a direct descendant of the corresponding species taxonomic node, which is necessary for the profiling tools.

## Selection of taxa for experiments
To select the taxa we considered here, we ran the `scripts/bacteria/select_taxa.py` script which selects an order, family and genus by consdiering species with at least 100 genomes, and focusing on intra-species diversity (see manuscript). The selected taxa and accessions from which we simulated reads for every experiment can be found in the files folder of this repository. The exception to this is the strain experiments, for which we used strains from an existing mock community sample (see manuscript).

## Index building
With all the reference sets, we can now build the profiling indices. In our work we used:
- Kraken2 (v2.1.3) + Bracken (v1.0.0)
- BWA (v0.7.18) + DUDes (v0.10.0)
- Centrifuge (v1.0.4.2)

### Kraken2 + Bracken
For Kraken2+Bracken we have provided a simple bash script called `compile_kraken2-bracken.sh` that gathers all target genomes, and builds a combined index. The script can be called for every experiment as follows:
```bash
bash scripts/bacteria/compile_kraken2-bracken.sh ${method}_${threshold} root/reference_sets/${experiment} root/indexes/${experiment}/bracken root/genomes root/taxdmp root/taxdmp_restricted/full.accession2taxid
```
This creates a combined Kraken2 and Bracken index for a given method-threshold combination and experiment, which is located in `/root/indexes/${experiment}/bracken/${method}_${threshold}`.
**NOTE**: for the strain-level experiments, we pass `root/taxdmp_restricted` instead of `root/taxdmp` to copy the custom taxonomy.

### Centrifuge
For Centrifuge we have a similar bash script called `compile_centrifuge.sh` which we call like:
```bash
bash scripts/bacteria/compile_centrifuge.sh ${method}_${threshold} root/reference_sets/${experiment} root/indexes/${experiment}/centrifuge root/genomes root/taxdmp root/taxdmp_restricted/partial.accession2taxid
```
**NOTE**: for the strain-level experiments, we pass `root/taxdmp_restricted` instead of `root/taxdmp` to copy the custom taxonomy.

### BWA + DUDes
Finally, for BWA with DUDes we use the script called `compile_bwa-dudes.sh` which can be called in a similar fashion:
```bash
bash scripts/bacteria/compile_bwa-dudes.sh ${method}_${threshold} root/reference_sets/${experiment} root/indexes/${experiment}/dudes root/genomes root/taxdmp root/taxdmp_restricted/full.accession2taxid
```
**NOTE**: for the strain-level experiments, we pass `root/taxdmp_restricted` instead of `root/taxdmp` to copy the custom taxonomy.

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
import argparse

def generate_commands(readcount):
	"""
	This function will generate bash commands to perform the ART simulations using the parameters defined as constants below.
	"""
	NUM_SAMPLES 	= 10
	NUM_READS 		= readcount #approximate number of reads per sample
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
                    output_path = f"{OUTPUT_PREFIX}/{experiment_type}/{readcount}/{sample}/{cur_species}" #this is output for reads per species
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
	parser = argparse.ArgumentParser()
	parser.add_argument("--readcount", type=int, required=True)
	args = parser.parse_args()

	num_reads = args.readcount

    # Currently, this uses 16 cores to process all the ART commands
    NUM_CORES = 16
    commands = generate_commands(num_reads)
    with multiprocessing.Pool(NUM_CORES) as pool:
        pool.map(run_subprocess, commands)
```
We ran this script with both `--readcount 500000` and `--readcount 5000000` to create 10 mixed metagenomic samples per experiment with approximately 500,000 paired-end and 5,000,000 reads of 150bp (270bp fragments, 10bp stdev) and equal read-based abundance for all species, respectively. The resulting reads are stored in `root/samples/{500000,5000000}/${experiment}/sample_${sample}/sample_{1,2}.fq`. With the simulated samples, we can call the profilers to estimate the relative abundances of species using the constructed reference sets.

### Kraken2 + Bracken
Running Kraken2 and Bracken is done in two subsequent steps: first running Kraken2 and then running Bracken to get species-level abundance estimates. This was done for every sample and method-threshold combination as follows:
```bash
kraken2 --db root/indexes/${experiment}/bracken/${method}_${threshold} --threads ${threads} --report root/estimations/{500000,5000000}/${experiment}/bracken/sample_${sample}/${method}_${threshold}.kreport --paired root/samples/{500000,5000000}/${experiment}/sample_${sample}/sample_1.fq root/samples/{500000,5000000}/${experiment}/sample_${sample}/sample_2.fq > root/estimations/{500000,5000000}/${experiment}/bracken/sample_${sample}/${method}_${threshold}.kraken #generates Kraken2 output
bracken -d root/indexes/${experiment}/bracken/${method}_${threshold} -i root/estimations/{500000,5000000}/${experiment}/bracken/sample_${sample}/${method}_${threshold}.kreport -o root/estimations/{500000,5000000}/${experiment}/bracken/sample_${sample}/${method}_${threshold}.bracken -r 150 -l S #generates Bracken output
```
This produces several output files, of which the `${method}_${threshold}.bracken` file is the most relevant as it will contain the species-level abundances (although the `${method}_${threshold}.kreport` files can be used to infer the number of unclassified reads). **NOTE**: for the strain-level experiments we instead run Bracken with `-l S1`.

### Centrifuge
Centrifuge can be run per sample, reference set and experiment as follows:
```bash
centrifuge -p ${threads} -x root/indexes/${experiment}/centrifuge/${method}_${threshold}/index/${method}_${threshold} -1 root/samples/{500000,5000000}/${experiment}/sample_${sample}/sample_1.fq -2 root/samples/{500000,5000000}/${experiment}/sample_${sample}/sample_2.fq -S root/estimations/{500000,5000000}/${experiment}/centrifuge/sample_${sample}/${method}_${threshold}.sam --report-file root/estimations/{500000,5000000}/${experiment}/centrifuge/sample_${sample}/${method}_${threshold}.report
```
This produces an alignment-like file (`${method}_${threshold}.sam`) as well as a report file (`${method}_${threshold}.report`) which will contain the abundance estimates.

### BWA + DUDes
DUDes also requires multiple steps. First we use BWA-mem to align reads to the reference index:
```bash
bwa mem -t ${threads} -v 1 root/indexes/${experiment}/dudes/${method}_${threshold}/bwa_index/bwa root/samples/{500000,5000000}/${experiment}/sample_${sample}/sample_1.fq root/samples/{500000,5000000}/${experiment}/sample_${sample}/sample_2.fq > root/estimations/{500000,5000000}/${experiment}/dudes/sample_${sample}/${method}_${threshold}.sam
```
After running BWA-mem we filter out unaligned reads and finally run DUDes on the filtered alignment file:
```bash
# Filter unaligned and calculate stats
samtools view -F 4 -h root/estimations/{500000,5000000}/${experiment}/dudes/sample_${sample}/${method}_${threshold}.sam > root/estimations/{500000,5000000}/${experiment}/dudes/sample_${sample}/${method}_${threshold}_filtered.sam
# Run DUDes
dudes -s root/estimations/{500000,5000000}/${experiment}/dudes/sample_${sample}/${method}_${threshold}_filtered.sam -d root/indexes/${experiment}/dudes/${method}_${threshold}/dudes_index/dudes.npz -o root/estimations/{500000,5000000}/${experiment}/dudes/sample_${sample}/${method}_${threshold}_dudes -l species
```
This results in a `${method}_${threshold}_dudes` file that contains the final abundance estimates. **NOTE**: for the strain-level experiments, we instead run DUDes with `-l strain`.

## Analysis
To analyse the results obtained we will assume that the workflow was ran as described above. Additionally, for resource usage monitoring, we assume that commands were run with `/usr/bin/time` and that the output is stored somewhere. **NOTE**: the Jupyter notebooks provided are for species-level experiments. The strain-level analysis was done in the same way, except using different thresholds, and profiling at the strain-level.

### Reference set comparisons
In `scripts/bacteria/analysis_reference_sets.ipynb` we provide a Jupyter notebook that details the steps we took to obtain the results and figures (using matplotlib and seaborn) regarding the reference set comparisons. 

### Accuracy comparisons
In `scripts/bacteria/analysis_accuracy.ipynb` we provide a Jupyter notebook that details the steps we took to obtain the results and figures (using matplotlib and seaborn) regarding the accuracy calculations.

### Runtime comparisons
In `scripts/bacteria/analsys_runtime.ipynb` we provide a Jupyter notebook that details the steps we took to obtain the results and figures (using matplotlib and seaborn) regarding the runtime calculations. **NOTE**: for the runtime analysis, we only consider the 5M readpair samples.

# Mock community samples for strain-level experiments
Here we will briefly describe the workflow for the mock community analysis in the strain-level experiments (both processing + read simulations).

**Dependencies**:
- bowtie2 v2.5.4
- fastp v1.1.0
- fastqc v0.12.1
- picard v2.20.4
- samtools v1.23
- sra-tools v3.2.1

## Pre-processing
We downloaded the SRR13355226 mock community sample, which consists of 99% human host DNA, and 1% *E. coli* DNA, as well as a human reference genome (GRCh38) by running the following steps:
```bash
ACC=SRR13355226

mkdir -p root/mock_community/fastq ## this is where reads and intermediate files will be stored
cd root/mock_community

prefetch $ACC ##requires sra_tools

## Convert to FastQ
fasterq-dump $ACC/$ACC.sra --split-files --threads ${threads} --outdir fastq

## Download human host from USCS
mkdir -p root/mock_community/human_ref
cd root/mock_community/human_ref
wget -c https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/latest/hg38.fa.gz
gunzip -c hg38.fa.gz > hg38.fa
```

Next, we run FastQC to assess the quality of reads:
```bash
mkdir -p root/mock_community/fastq/qc/raw

fastqc -t ${threads} -o root/mock_community/fastq/qc/raw root/mock_community/fastq/*.fastq
```

The result of running FastQC was the following:
```
==== SRR13355226_1_fastqc/summary.txt ====
PASS    Basic Statistics        SRR13355226_1.fastq
PASS    Per base sequence quality       SRR13355226_1.fastq
WARN    Per tile sequence quality       SRR13355226_1.fastq
PASS    Per sequence quality scores     SRR13355226_1.fastq
FAIL    Per base sequence content       SRR13355226_1.fastq
FAIL    Per sequence GC content SRR13355226_1.fastq
PASS    Per base N content      SRR13355226_1.fastq
PASS    Sequence Length Distribution    SRR13355226_1.fastq
PASS    Sequence Duplication Levels     SRR13355226_1.fastq
PASS    Overrepresented sequences       SRR13355226_1.fastq
FAIL    Adapter Content SRR13355226_1.fastq
==== SRR13355226_2_fastqc/summary.txt ====
PASS    Basic Statistics        SRR13355226_2.fastq
FAIL    Per base sequence quality       SRR13355226_2.fastq
WARN    Per tile sequence quality       SRR13355226_2.fastq
PASS    Per sequence quality scores     SRR13355226_2.fastq
FAIL    Per base sequence content       SRR13355226_2.fastq
FAIL    Per sequence GC content SRR13355226_2.fastq
PASS    Per base N content      SRR13355226_2.fastq
PASS    Sequence Length Distribution    SRR13355226_2.fastq
PASS    Sequence Duplication Levels     SRR13355226_2.fastq
PASS    Overrepresented sequences       SRR13355226_2.fastq
FAIL    Adapter Content SRR13355226_2.fastq
```

This shows that adapters need to be filtered, and that tails need to be trimmed. Additionally, we perform minimal read length filtering by running FastP:
```bash
mkdir -p root/mock_community/fastq/qc/fastp

BASEDIR=root/mock_community/fastq
fastp \
	-i ${BASEDIR}/SRR13355226_1.fastq \
    -I ${BASEDIR}/SRR13355226_2.fastq \
    -o ${BASEDIR}/qc/fastp/SRR13355226_1.fastq \
    -O ${BASEDIR}/qc/fastp/SRR13355226_2.fastq \
    -h ${BASEDIR}/qc/fastp/fastp_report.html \
    -j ${BASEDIR}/qc/fastp/fastp_report.json \
    -w ${SLURM_CPUS_PER_TASK} \
    --detect_adapter_for_pe \
    --cut_tail \
    --cut_tail_mean_quality 20 \
    --length_required 50 
```

The next step is to map reads to the human genome in order to 1: de-host the sample and 2: estimate read fragment statistics for simulations:
```bash
BASEDIR=root/mock_community
ACC=SRR13355226

## Build index
bowtie2-build ${BASEDIR}$/human_ref/hg38.fa ${BASEDIR}/human_ref/bowtie2/hg38

mkdir -p ${BASEDIR}/dehosting/dehosted
mkdir -p ${BASEDIR}/dehosting/host

bowtie2 \
	--very-sensitive \
	-x ${BASEDIR}/human_ref/bowtie2/hg38 \
	-1 ${BASEDIR}/fastq/qc/${ACC}_1.fastq \
	-2 ${BASEDIR}/fastq/qc/${ACC}_2.fastq \
	-p ${threads} \
	--un-conc ${BASEDIR}/dehosting/dehosted/dehosted.fastq \
	-S - \
	| samtools view -@ ${threads} -b - \
	| samtools sort -@ ${threads} -o ${BASEDIR}/dehosting/host/host.sorted.bam -

samtools index -@ ${threads} ${BASEDIR}/dehosting/host/host.sorted.bam

samtools fastq ${BASEDIR}/dehosting/host/host.sorted.bam \
	-f 12 -F 256 \
	-1 ${BASEDIR}/dehosting/dehosted/dehosted_unmapped_1.fastq \
	-2 ${BASEDIR}/dehosting/dehosted/dehosted_unmapped_2.fastq \
	-0 /dev/null -s /dev/null -n
```
 
The final step in pre-processing, is to estimate the stats of the aligned reads using Picard:
```bash
BASEDIR=root/mock_community

# First retrieve proper pairs with mapping quality >= 20
samtools view -b -f 2 -q 20 -F 2304 \
	${BASEDIR}/dehosting/host/host.sorted.bam \
	| samtools sort -@ ${threads} -o ${BASEDIR}/dehosting/host/host.proper.q20.bam -

# Index
samtools index -@ ${threads} ${BASEDIR}/dehosting/host/host.proper.q20.bam

# Picard
picard CollectInsertSizeMetrics \
	I=${BASEDIR}/dehosting/host/host.proper.q20.bam \
    O=${BASEDIR}/dehosting/host/insert_metrics.txt \
    H=${BASEDIR}/dehosting/host/insert_hist.pdf \
    M=0.5
```

From this we find that the average fragment size is 265.34 with a standard deviation of 104.33, which we will use as input for the mock sample simulations that we will describe below.

Below is the Python script that we used to generate the mock community samples:
```python
import subprocess
import os
import random
import shutil
from pathlib import Path
import math
import argparse
import multiprocessing
from Bio import SeqIO

def main():
	"""
	This script will generate simulated reads for the strain-level experiments.
	We generate 3 sets of samples:
		3. 10 samples with properties similar to the mock community sample
	"""
	parser = argparse.ArgumentParser()
	# Overall parameters
	parser.add_argument("--genomes_folder", type=str, required=True)
	# Mock sample parameters
	parser.add_argument("--num_mock_samples", type=int, default=10)
	parser.add_argument("--num_mock_reads", type=int, default=5_000_000)
	parser.add_argument("--mock_fragment_size", type=int, default=270)
	parser.add_argument("--mock_fragment_std", type=int, default=20)
	parser.add_argument("--mock_read_length", type=int, default=151)
	parser.add_argument("--mock_technology", type=str, default="HS25")
	parser.add_argument("--mock_output_prefix", type=str, required=True)
	# Coverage sample parameters
	parser.add_argument("--cov_output_prefix", type=str, required=True)
	# Temporaryy directory for ART output (will be removed at the end of the script)
	parser.add_argument("--tmp_folder", type=str, required=True)
	args = parser.parse_args()

	# Constants (strains, abundances for mock sample, etc.)
	STRAINS = [316401, 364106, 386585, 331111] #H10407, UTI89, Sakai, E24377A
	ABUNDANCES = {
		316401: "0.800",	#H10407
		364106: "0.150", 	#UTI89
		386585: "0.049", 	#Sakai
		331111: "0.001", 	#E24377A
	}

	# Find genome lenghts (which is used to determine fold coverage)
	length_per_strain = {} #only one genome per strain, so we can just sum the lengths of all contigs for that genome
	for strain in STRAINS:
		length_per_strain[strain] = 0
		with open(f"{args.genomes_folder}/{strain}.fna", "r") as f_in:
			for record in SeqIO.parse(f_in, "fasta"):
				length_per_strain[strain] += len(record.seq)

	############### Generate mock community samples ###############
	print("Generating mock community samples...")
	seed = args.num_mock_reads #use a different seed from the baseline experiments
	# Here we do not use math.ceil since we want to closely mimic the properties of the mock community sample
	num_reads_per_strain = {
		strain: args.num_mock_reads * float(ABUNDANCES[strain]) for strain in STRAINS
	}
	for sample in range(1, args.num_mock_samples+1):
		output_path = f"{args.tmp_folder}/samples_mock/{sample}"
		output_path_object = Path(output_path)
		output_path_object.mkdir(parents=True, exist_ok=True)
		# Generate commands for current sample
		fwd_paths = []
		rev_paths = []
		for s_idx, strain in enumerate(STRAINS):
			"""
			Below we define the fold coverage as 2 x (read length) x (number of reads) / (genome length) since this is more consistent with
			the definition internally used by ART. For all other experiments, where the fold coverage is less relevant (here we
			want to closely mimic the properties of the mock community sample), we will stick to the fragment size x number of reads / genome length definition, 
			which we use for all other experimental settings as well.
			"""
			command = [
				"art_illumina",
				"-p", #paired-end
				"-i", f"{args.genomes_folder}/{strain}.fna",
				"-l", str(args.mock_read_length),
				"-f", f"{num_reads_per_strain[strain] * args.mock_read_length * 2 / length_per_strain[strain]:.5f}", #calculate fold coverage based on number of reads, read length and genome length
				"-ss", args.mock_technology,
				"-m", str(args.mock_fragment_size),
				"-s", str(args.mock_fragment_std),
				"-o", f"{output_path}/{strain}_",
				"-rs", str(seed),
				"-na"
			]
			fwd_paths.append(f"{output_path}/{strain}_1.fq")
			rev_paths.append(f"{output_path}/{strain}_2.fq")
			subprocess.run(command, check=True) #runs command and raises an error if it fails
			seed += 1
		
		# Combine forward reads
		output_path = f"{args.mock_output_prefix}/samples_mock/{sample}"
		os.makedirs(output_path, exist_ok=True)
		with open(f"{output_path}/sample_1.fq", "w") as f_out:
			subprocess.run(["cat", *fwd_paths], stdout=f_out, check=True)
		# Combine reverse reads
		with open(f"{output_path}/sample_2.fq", "w") as f_out:
			subprocess.run(["cat", *rev_paths], stdout=f_out, check=True)
		# Remove temporary folder with individual .fq files
		shutil.rmtree(f"{args.tmp_folder}/samples_mock/{sample}")

		print(f"Finished generating mock community sample {sample}/{args.num_mock_samples}")

if __name__ == "__main__":
	main()
```

The script was called as:
```bash
python -u generate_mock_samples.py \
    --genomes_folder ${PREFIX}/genomes \
    --num_mock_samples 10 \
    --num_mock_reads 726777 \
    --mock_fragment_size 265 \
    --mock_fragment_std 104 \
    --mock_read_length 150 \
    --mock_technology "HS25" \
    --mock_output_prefix ${PREFIX} \
```
where we assume that ${PREFIX} points to a folder that contains a subfolder called `genomes` with 1 genome for every strain in the mock community (see manuscript).

The analysis carried out for the mock community sample (both real and simulated) was identical to the other simulated samples, using the analysis scripts adapted for strain-level experiments.