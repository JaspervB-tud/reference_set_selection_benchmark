METHOD=$1
SELECTION_FOLDER=$2
INDEX_FOLDER=$3
GENOMES_FOLDER=$4
TAXONOMY_FOLDER=$5

# Create output folders
mkdir -p "${INDEX_FOLDER}/${METHOD}"
mkdir -p "${INDEX_FOLDER}/${METHOD}/index"
mkdir -p "${INDEX_FOLDER}/${METHOD}/genomes"
mkdir -p "${INDEX_FOLDER}/${METHOD}/taxonomy"

# Centrifuge requires a particular format for the taxonomy mapping (seqid \t taxid)
awk 'NR > 1 {print $2 "\t" $3}' "${SELECTION_FOLDER}/nucl_gb.accession2taxid" > "${INDEX_FOLDER}/${METHOD}/taxonomy/nucl_gb.accession2taxid"

# Add genomes to the database and remove individual genomes afterwards, saving only the aggregated genome file
while IFS=$'\t' read -r SPECIES SEQUENCE col3; do
    cp "${GENOMES_FOLDER}/${SPECIES}/${SEQUENCE}" "${INDEX_FOLDER}/${METHOD}/genomes/"
    cat "${INDEX_FOLDER}/${METHOD}/genomes/${SEQUENCE}" >> "${INDEX_FOLDER}/${METHOD}/genomes/all_genomes.fasta"
    rm "${INDEX_FOLDER}/${METHOD}/genomes/${SEQUENCE}"
done < "${SELECTION_FOLDER}/${METHOD}.tsv"

# Build index
centrifuge-build --conversion-table ${INDEX_FOLDER}/${METHOD}/taxonomy/nucl_gb.accession2taxid \
    --taxonomy-tree ${TAXONOMY_FOLDER}/nodes.dmp \
    --name-table ${TAXONOMY_FOLDER}/names.dmp \
    "${INDEX_FOLDER}/${METHOD}/genomes/all_genomes.fasta" "${INDEX_FOLDER}/${METHOD}/index/"

# Remove genomes after building
rm -r "${INDEX_FOLDER}/${METHOD}/genomes"