METHOD=$1
SELECTION_FOLDER=$2
INDEX_FOLDER=$3
GENOMES_FOLDER=$4
TAXONOMY_FOLDER=$5

# Create output folders
mkdir -p "${INDEX_FOLDER}/${METHOD}"
mkdir -p "${INDEX_FOLDER}/${METHOD}/bwa_index"
mkdir -p "${INDEX_FOLDER}/${METHOD}/dudes_index"
mkdir -p "${INDEX_FOLDER}/${METHOD}/genomes"

# Add genomes to the database and remove individual genomes afterwards, saving only the aggregated genome file
while IFS=$'\t' read -r SPECIES SEQUENCE col3; do
    cp "${GENOMES_FOLDER}/${SPECIES}/${SEQUENCE}" "${INDEX_FOLDER}/${METHOD}/genomes/"
    cat "${INDEX_FOLDER}/${METHOD}/genomes/${SEQUENCE}" >> "${INDEX_FOLDER}/${METHOD}/genomes/all_genomes.fasta"
    rm "${INDEX_FOLDER}/${METHOD}/genomes/${SEQUENCE}"
done < "${SELECTION_FOLDER}/${METHOD}.tsv"

# Build BWA index
bwa index -b 50000000 "${INDEX_FOLDER}/${METHOD}/genomes/all_genomes.fasta" \
    -p "${INDEX_FOLDER}/${METHOD}/bwa_index/bwa"

# Build DUDes index
dudesdb -m 'av' -f "${INDEX_FOLDER}/${METHOD}/genomes/all_genomes.fasta" \
    -n "${TAXONOMY_FOLDER}/nodes.dmp" \
    -a "${TAXONOMY_FOLDER}/names.dmp" \
    -g "${SELECTION_FOLDER}/nucl_gb.accession2taxid" \
    -o "${INDEX_FOLDER}/${METHOD}/dudes_index/dudes"

# Remove genomes after building
rm -r "${INDEX_FOLDER}/${METHOD}/genomes"