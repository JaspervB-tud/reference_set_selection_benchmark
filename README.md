# Benchmarking the impact of reference genome selection on taxonomic profiling accuracy
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

The manuscript itself can be found on [bioRxiv](https://doi.org/10.1101/2025.02.07.637076).

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
├── reference_sets
├── indexes
│   ├── profiler_1
│   |   ├── method_1
│   |   ├── method_2-threshold_1
|   |   └── method_2-threshold_2
|   ├── profiler_2
│   |   ├── method_1
│   |   ├── method_2-threshold_1
|   |   └── method_2-threshold_2
├── samples
│   ├── sample_1
│   |   ├── sample_1_1.fq
|   |   ├── sample_1_2.fq
|   |   ├── sample_1_groundtruth_folds.tsv
|   |   └── sample_1_groundtruth_reads.tsv
├── estimations
│   ├── sample_1
│   |   ├── profiler_1
|   |   └── profiler_2
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

### VLQ selection
For selecting SARS-CoV-2 genomes, we also used the selection methodology provided in the VLQ pipeline, however, we use a slightly adapted version of the `call_variants.sh` script which now operates on a single lineage. The script should be called as follows for every lineage:
```bash
bash scripts/run_vlq.sh LINEAGE REF_GENOME genomes/ full_path_to_selections_folder/vlq
```
Here, LINEAGE is the name of the lineage (e.g. X in most of our example code so far), REF_GENOME is the reference genome against which variants are called (MN908947.3 for our experiments). Running the script will generate a `LINEAGE_merged.frq` file, `LINEAGE_merged.pi` file and a `LINEAGE_merged.vcf.gz` file for the lineage on which the script is called. 

In order to perform VLQ's selection steps, this should be run on ALL lineages before proceeding, making sure that the `selections/vlq` folder contains the required files for all lineages. Afterwards we call the `select_samples.py` script from VLQ which requires a single fasta file containing all genomes (from every lineage), and a corresponding metadata file. Assuming that these are stored in the `genomes` folder, we run:
```bash
python select_samples.py -m genomes/metadata.tsv -f genomes/all_genomes.fasta -o selections/vlq --vcf selections/vlq/*_merged.vcf.gz --freq selections/vlq/*_merged.fq --max_per_lineage 1000
```
This generates both a `sequences.fasta` and `metadata.tsv` file in the `selections/vlq` folder, containing the genome sequences and the metadata respectively, for every lineage.

### Post processing
After running the tools, all selected genomes should be located in the `selections` folder. At this point we used the `generate_all_selection.py` script to write a summary file with an "all" selection that includes every genome. Furthermore, this will also generate an accession2taxid file in the same format as the accession2taxid files on NCBI Taxonomy. This will be passed to the taxonomic profilers in order to match the ids of individual sequences to the corresponding species taxids. In principle this can be mapped to more specific taxids (e.g. sub-species or strain level if available), but in our experiments we performed experiments at the species level. This can be ran as follows:
```bash
python scripts/generate_all_selection.py --genomes genomes --output reference_sets --a2t
```
This generates two files in the `reference_sets` folder: a `nucl_gb.accession2taxid` file which contains sequence id to taxonomy id mappings, and a file called `all.tsv` which contains the filenames of all genome files, for every species. Every line in this tab delimited file contains a genome, and is structured as:
```
SPECIES_TAXID   FILENAME    +/-
```
The final column contains a "+" if the corresponding genome was selected, and a "-" if it was randomly selected in case the selection algorithm failed to perform a selection, or if there was only a single available genome.

Alternatively, one can run:
```bash
python scripts/generate_all_selection.py --genomes genomes --output reference_sets
```
for the SARS-CoV-2 setting, since in that case the accession2taxid mapping file is not necessary due to the availability of metadata.

After running the `generate_all_selection.py` script, we now process every selection to produce a similar output which is later used to fetch all the genomes when building the profiling index and when comparing reference sets. For this we run the following:
```bash
python scripts/generate_selection_files.py --genomes genomes --selection selections/METHOD/THRESHOLD --filename METHOD_THRESHOLD --output reference_sets
```
Running this will produce an output file called `METHOD_THRESHOLD.tsv` (with METHOD and THRESHOLD substituted for the chosen method and threshold combination), formatted as the `all.tsv` file, in the `reference_sets` folder. This script should be run for all selection methods and thresholds to produce all the corresponding selection files.

### Comparing reference sets
To compare the generated sequence sets for every method (both for bacteria and SARS-CoV-2) we used the generated `METHOD_THRESHOLD.tsv` files and calculated for every pair of methods (including the "all" selection) the containment index (in both directions) only considering non-singleton species for which both methods were able to produce a selection (i.e. both methods will have a "+" in the final column of the `.tsv` file).

## Index building
Having prepared all the selections and related output files, these can now be used to build the indexes for taxonomic profilers. In our work we used:
- Kraken2 (v2.1.3) + Bracken (v1.0.0)
- BWA (v0.7.18) + DUDes (v0.10.0)
- Centrifuge (v1.0.4.2)
in the bacterial setting, and
- VLQ (kallisto v0.44.0)
for the SARS-CoV-2 setting. Below we will describe how we built the indexes for every profiler.

### Bacteria
#### Kraken2 + Bracken
For Kraken2 in combination with Bracken we have provided a simple bash script called `compile_kraken2-bracken.sh` which gathers all target genomes, and builds a combined Kraken2 and Bracken index with all of the target genomes. The script can be called as follows:
```bash
bash scripts/compile_kraken2-bracken.sh METHOD_THRESHOLD reference_sets indexes/kraken2-bracken genomes taxonomy
```
This assumes that the taxonomy folder is an effective copy of the NCBI taxdmp folder (version used in manuscript was downloaded on September 24, 2024) such that it contains at least the `names.dmp` and `nodes.dmp` files, and that the created `nucl_gb.accession2taxid` file is in the `reference_sets` folder. After running, the resulting index can be found in the `indexes/kraken2-bracken/METHOD_THRESHOLD` where METHOD and THRESHOLD are placeholders for the method and threshold respectively.

#### Centrifuge
For Centrifuge we can use the provided bash script called `compile_centrifuge.sh` which effectively mimics the behavior of the `compile_kraken2-bracken.sh` script:
```bash
bash scripts/compile_centrifuge.sh METHOD_THRESHOLD reference_sets indexes/centrifuge genomes taxonomy
```

#### BWA + DUDes
Finally, for BWA and DUDes we use the script called `compile_bwa-dudes.sh`, which does the same:
```bash
bash scripts/compile_bwa-dudes.sh METHOD_THRESHOLD reference_sets indexes/bwa-dudes genomes taxonomy
```

### SARS-CoV-2
#### Kallisto
For Kallisto, which forms the backbone of the VLQ pipeline, we have provided a script called `compile_kallisto.sh` which can be ran as:
```bash
bash scripts/compile_kallisto.sh METHOD_THRESHOLD reference_sets indexes/kallisto genomes
```

## Profiling
In our manuscript we used all profilers to estimate a taxonomic profile of simulated paired-end metagenomic samples. In both the Bacteria and SARS-CoV-2 cases these reads were 150bp long, but fragments sizes differ (see manuscript). The simulated reads themselves are available on [Zenodo](https://doi.org/10.5281/zenodo.14727633).

### Bacteria
#### Kraken2 + Bracken
Running Kraken2 and Bracken is done in two steps: first running Kraken2, and running Bracken at a desired taxonomic level (species here) afterwards. Assuming the set-up provided here this can be done for every METHOD_THRESHOLD as follows:
```bash
kraken2 --db indexes/kraken2-bracken/METHOD_THRESHOLD --report estimations/sample_1/kraken2-bracken/METHOD_THRESHOLD.kreport --paired samples/sample_1/sample_1_1.fq samples/sample_1/sample_1_2.fq > estimations/sample_1/kraken2-bracken/METHOD_THRESHOLD.kraken
bracken -d indexes/kraken2-bracken/METHOD_THRESHOLD -i estimations/sample_1/kraken2-bracken/METHOD_THRESHOLD.kreport -o estimations/sample_1/kraken2-bracken/METHOD_THRESHOLD.bracken -r 150 -l S
```

#### Centrifuge
Centrifuge can be run as follows:
```bash
centrifuge -x indexes/centrifuge/METHOD_THRESHOLD -1 samples/sample_1/sample_1_1.fq -2 samples/sample_1/sample_1_2.fq -S estimations/sample_1/centrifuge/METHOD_THRESHOLD.sam --report-file estimations/sample_1/centrifuge/METHOD_THRESHOLD.report
```

#### BWA + DUDes
DUDes relies on an alignment file which is used to perform taxonomic classification. We have used BWA MEM for this:
```bash
bwa mem -v 1 indexes/bwa-dudes/METHOD_THRESHOLD/bwa_index/bwa samples/sample_1/sample_1_1.fq -2 samples/sample_1/sample_1_2.fq > estimations/sample_1/bwa-dudes/METHOD_THRESHOLD.sam
```

After running BWA MEM, we filter out unaligned reads, and we store the number of reads that were aligned for later processing. Finally, we run DUDes on the filtered alignment file
```bash
# Filter unaligned and calculate stats
samtools view -F 4 -h estimations/sample_1/bwa-dudes/METHOD_THRESHOLD.sam > estimations/sample_1/bwa-dudes/METHOD_THRESHOLD_filtered.sam
samtools view -c -f 4 estimations/sample_1/bwa-dudes/METHOD_THRESHOLD.sam estimations/sample_1/bwa-dudes/METHOD_THRESHOLD_num-unaligned.txt
samtools view -c -F 4 estimations/sample_1/bwa-dudes/METHOD_THRESHOLD.sam estimations/sample_1/bwa-dudes/METHOD_THRESHOLD_num-aligned.txt
# Run DUDes
dudes -s estimations/sample_1/bwa-dudes/METHOD_THRESHOLD_filtered.sam -d indexes/bwa-dudes/METHOD_THRESHOLD/dudes_index/dudes.npz -o estimations/sample_1/bwa-dudes/METHOD_THRESHOLD_dudes -l species
```

### SARS-CoV-2
In order to obtain SARS-CoV-2 lineage abundance estimates we first run kallisto, followed by the `output_abundances.py` script from VLQ:
```bash
kallisto quant -b 0 -i indexes/kallisto/METHOD_THRESHOLD/index.idx -o estimations/sample_1/kallisto/METHOD_THRESHOLD samples/sample_1/sample_1_1.fq samples/sample_1/sample_1_2.fq
python output_abundances.py -m 0.1 -o estimations/sample_1/kallisto/METHOD_THRESHOLD/predictions.tsv --metadata genomes/metadata.tsv estimations/sample_1/kallisto/METHOD_THRESHOLD/abundance.tsv
```

## Analyzing results
Accuracy metrics were analyzed using the `analyze_bacteria.ipynb` and `analyze_sc2.ipynb` Jupyter notebooks available in the `scripts` folder. In the bacteria case, `analyze_bacteria.ipynb` directly calculates the number of unclassified reads since this is directly available in the output of Bracken. For BWA in combination with DUDes, we instead calculated them from the produced `_num-unaligned` files and for Centrifuge they were obtained from the generated `.sam` files using:
```bash
cut -f2 "estimations/sample_1/centrifuge/METHOD_THRESHOLD.sam" | awk '$1 == "unclassified" {count++} END {printf "%d, ", count}'
```
In the SARS-CoV-2 case, the percentage of unaligned reads can directly be obtained from the `run_info.json` file that is produced when running kallisto.
