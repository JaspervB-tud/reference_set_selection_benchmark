The goal of selection is to find a subset of genomes that collectively represent all of the available genomes for a species. We have done this with the following approaches:

- Centroid selection based on MASH distances
- Hierarchical clustering based on MASH distances (as in dRep)
- GGRaSP based on MASH distances
- Gclust
- MeShClust

After selecting reference genomes, we assess the performance of taxonomic profilers based on the selected reference genomes. In our work we have done this for:

- Kraken2 + Bracken
- Centrifuge
- BWA + DUDes

In what follows we will provide the steps taken to generate the results in our manuscript.

## Genome selection
### MASH-based selection
We first generate MASH sketches, which are used to estimate the (dis)similarity between genomes. For the example provided here, all genomes are decompressed and stored in the "genomes_per_species" folder, where every subfolder is labeled according to the (fake) taxid of the species and contains the genomes of the species.

Assuming that the user is in this folder, we would now run the following two lines for every species (we only show the lines for running it for species "0"):

```bash
mash sketch -s 1000 -S 123456 -k 21 -o genomes_per_species/0/mash_sketches.msh genomes_per_species/0/*.fa
mash triangle genomes_per_species/0/mash_sketches.msh > genomes_per_species/0/mash_distances.dist
```

After doing so, we obtain a MASH sketch for every genome (located in the `mash_sketches.msh` files), as well as a lower triangular distance matrix between all genomes of a species.

#### Hierarchical clustering selection
To obtain the hierarchical clustering selection, we can now run `run_hierarchical.py`, storing the results in the "selections" folder. In the manuscript, we ran the script with both single-linkage and complete-linkage linkages, and thresholds equal to 0.01, 0.03 and 0.05 which correspond to similarities of 99%, 97% and 95% respectively. The following snippet shows how this was done (for species 0) with all three thresholds:

```bash
python ../../scripts-bacteria/run_hierarchical.py --matrix genomes_per_species/0/mash_distances.dist --threshold 0.01 --output selections/0
python ../../scripts-bacteria/run_hierarchical.py --matrix genomes_per_species/0/mash_distances.dist --threshold 0.03 --output selections/0
python ../../scripts-bacteria/run_hierarchical.py --matrix genomes_per_species/0/mash_distances.dist --threshold 0.05 --output selections/0
```
After running, the selections can be found in `./selections/0/TYPE-linkage_THRESHOLD`.

#### Centroid selection
The centroid selection can be obtained by running `run_centroid.py`, and storing the results in the "selections" folder again.
```bash
python ../../scripts-bacteria/run_centroid.py --matrix genomes_per_species/0/mash_distances.dist --output selections/0
```
After running, the selection can be found in `./selections/0/centroid`.

#### GGRaSP selection
GGRaSP expects the output to be formatted as a full matrix, rather than a (lower) triangular matrix. To this end we can use the `convert_matrix.py` script to convert the distance matrix obtained from MASH. Afterwards, we use the `run_ggrasp.R` script to run GGRaSP. NOTE: GGRaSP is prone to fail for various different reasons. In the scenario where GGRaSP failed (or any other tool for that matter), we randomly picked a single representative for every taxon where it failed in order to include the corresponding taxon in the reference database.
```bash
python ../../scripts-bacteria/convert_matrix.py --matrix genomes_per_species/0/mash_distances.dist --output genomes_per_species/0
Rscript ../../scripts-bacteria/run_ggrasp.R 0
```
After running, the selection can be found in `./selections/0/ggrasp`.

### Gclust and MeShClust
Both Gclust and MeShClust assume that the input consists of a single multi-fasta file where every entry corresponds to a genome. Since several genomes in our experiments had more than one entry in their fasta file (extra-chromosomal DNA for example), we therefore first call the script `concatenate_genomes.py` in order to create the required input file. The entries of a fasta file were concatenated in the order in which they appeared.
```bash
python ../../scripts-bacteria/concatenate_genomes.py --genomes genomes_per_species/0 --output genomes_per_species/0
```

#### Gclust
In addition to requiring that all genomes are provided in a single multi-fasta file, Gclust also requires users to run a sorting script that sorts the genomes from longest to shortest:
```bash
perl sortgenome.pl --genomes-file genomes_per_species/0/all_genomes.fasta --sortedgenomes-file genomes_per_species/0/all_genomes_sorted.fasta
```
After sorting we run Gclust for similarities 95%, 97% and 99% respectively as follows:
```bash
gclust -minlen 41 -both -chunk 400 -ext 1 -sparse 4 -nuc -loadall -rebuild -memiden 95 genomes_per_species/0/all_genomes_sorted.fasta > selections/0/gclust_0.95
gclust -minlen 41 -both -chunk 400 -ext 1 -sparse 4 -nuc -loadall -rebuild -memiden 97 genomes_per_species/0/all_genomes_sorted.fasta > selections/0/gclust_0.97
gclust -minlen 41 -both -chunk 400 -ext 1 -sparse 4 -nuc -loadall -rebuild -memiden 99 genomes_per_species/0/all_genomes_sorted.fasta > selections/0/gclust_0.99
```
This uses the same parameters as used for the results in Table 3 of the original Gclust paper, and after running the selections can be found in `./selections/0/gclust_THRESHOLD`.

#### MeShClust
MeShClust can be ran directly on the multi-fasta file that includes all of the genomes. We ran MeShClust using default parameters and similarity thresholds of 95%, 97% and 99% respectively as follows:
```bash
meshclust -d genomes_per_species/0/all_genomes.fasta -o selections/0/meshclust_0.95 -t 0.95
meshclust -d genomes_per_species/0/all_genomes.fasta -o selections/0/meshclust_0.97 -t 0.97
meshclust -d genomes_per_species/0/all_genomes.fasta -o selections/0/meshclust_0.99 -t 0.99
```
After running, the selections can be found in `./selections/0/meshclust_THRESHOLD`.