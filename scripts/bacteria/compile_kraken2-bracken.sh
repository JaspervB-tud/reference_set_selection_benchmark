METHOD=$1
SELECTION_FOLDER=$2
INDEX_FOLDER=$3
GENOMES_FOLDER=$4
TAXONOMY_FOLDER=$5
A2T_LOC=$6

# Create output folders
mkdir -p "${INDEX_FOLDER}/${METHOD}"
mkdir -p "${INDEX_FOLDER}/${METHOD}/taxonomy"
mkdir -p "${INDEX_FOLDER}/${METHOD}/library"
mkdir -p "${INDEX_FOLDER}/${METHOD}/genomes"

# Copy taxonomy files
cp ${TAXONOMY_FOLDER}/* "${INDEX_FOLDER}/${METHOD}/taxonomy"
cp ${A2T_LOC} "${INDEX_FOLDER}/${METHOD}/taxonomy"

# Add genomes to the database and remove individual genomes afterwards, saving only the aggregated genome file
while IFS=$'\t' read -r SPECIES SEQUENCE col3; do
    cp "${GENOMES_FOLDER}/${SPECIES}/${SEQUENCE}" "${INDEX_FOLDER}/${METHOD}/genomes/"
    cat "${INDEX_FOLDER}/${METHOD}/genomes/${SEQUENCE}" >> "${INDEX_FOLDER}/${METHOD}/genomes/all_genomes.fasta"
    rm "${INDEX_FOLDER}/${METHOD}/genomes/${SEQUENCE}"
done < "${SELECTION_FOLDER}/${METHOD}.tsv"

# Add to Kraken2 library (default parameters, number of threads can be specified by user using --threads flag)
kraken2-build --add-to-library "${INDEX_FOLDER}/${METHOD}/genomes/all_genomes.fasta" --db "${INDEX_FOLDER}/${METHOD}"
# Remove fasta, no longer needed here
rm -r ${INDEX_FOLDER}/${METHOD}/genomes
# Build index and cleanup afterwards
kraken2-build --build --db "${INDEX_FOLDER}/${METHOD}"
#/tudelft.net/staff-umbrella/refsetopt/Bracken/bracken-build -d "${INDEX_FOLDER}/${METHOD}" -l 150 #we work with reads of 150bp
kraken2-build --clean --db "${INDEX_FOLDER}/${METHOD}"