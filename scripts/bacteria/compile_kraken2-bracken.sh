METHOD=$1
REFSET_FOLDER=$2
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
cp ${A2T_LOC} "${INDEX_FOLDER}/${METHOD}/taxonomy/nucl_gb.accession2taxid"

line_num=0
# Add genomes to the database and remove individual genomes afterwards, saving only the aggregated genome file
while IFS=$'\t' read -r SPECIES SEQUENCE col3; do
    # Skip header line
    ((line_num++))
    if [ $line_num -eq 1 ]; then
        continue
    fi
    cp "${GENOMES_FOLDER}/${SPECIES}/${SEQUENCE}" "${INDEX_FOLDER}/${METHOD}/genomes/"
    if [[ "${SEQUENCE}" == *.gz ]]; then # unzip if necessary
        gunzip "${INDEX_FOLDER}/${METHOD}/genomes/${SEQUENCE}"
        SEQUENCE="${SEQUENCE%.gz}"
    fi
    cat "${INDEX_FOLDER}/${METHOD}/genomes/${SEQUENCE}" >> "${INDEX_FOLDER}/${METHOD}/genomes/all_genomes.fasta"
    rm "${INDEX_FOLDER}/${METHOD}/genomes/${SEQUENCE}"
done < "${REFSET_FOLDER}/${METHOD}.tsv"

# Add to Kraken2 library (default parameters, number of threads can be specified by user using --threads flag)
kraken2-build --add-to-library "${INDEX_FOLDER}/${METHOD}/genomes/all_genomes.fasta" --db "${INDEX_FOLDER}/${METHOD}" --threads 32
# Remove fasta, no longer needed here
rm -r ${INDEX_FOLDER}/${METHOD}/genomes
# Build index
kraken2-build --build --db "${INDEX_FOLDER}/${METHOD}"
bracken-build -d "${INDEX_FOLDER}/${METHOD}" -l 150 -t 32
# Clean up intermediate files
kraken2-build --clean --db "${INDEX_FOLDER}/${METHOD}"