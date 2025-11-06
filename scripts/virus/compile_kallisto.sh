METHOD=$1
SELECTION_FOLDER=$2
INDEX_FOLDER=$3
GENOMES_FOLDER=$4

# Create output folders
mkdir -p "${INDEX_FOLDER}/${METHOD}/genomes"

line_num=0
# Add genomes to the database and remove individual genomes afterwards, saving only the aggregated genome file
while IFS=$'\t' read -r LINEAGE SEQUENCE col3; do
    # Skip header line
    ((line_num++))
    if [ $line_num -eq 1 ]; then
        continue
    fi
    cp "${GENOMES_FOLDER}/${LINEAGE}/${SEQUENCE}" "${INDEX_FOLDER}/${METHOD}/genomes/"
    cat "${INDEX_FOLDER}/${METHOD}/genomes/${SEQUENCE}" >> "${INDEX_FOLDER}/${METHOD}/genomes/all_genomes.fasta"
    rm "${INDEX_FOLDER}/${METHOD}/genomes/${SEQUENCE}"
done < "${SELECTION_FOLDER}/${METHOD}.tsv"

# Build kallisto index
kallisto index -i "${INDEX_FOLDER}/${METHOD}.idx" "${INDEX_FOLDER}/${METHOD}/genomes/all_genomes.fasta"

# Remove genomes
rm -r "${INDEX_FOLDER}/${METHOD}"