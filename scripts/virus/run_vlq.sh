LINEAGE=$1
REF_GENOME=$2
BASE_DIR=$3
OUTPUT_DIR=$4

# Make directory for storing reference genomes
mkdir -p ${OUTPUT_DIR}

# Move to directory
cd "${BASE_DIR}/${LINEAGE}"
for FASTA in *.fa; do
    minimap2 -c -x asm20 --end-bonus 100 -t 20 --cs $REF_GENOME $FASTA 2>${FASTA%.fa}.paftools.log | sort -k6,6 -k8,8n > ${FASTA%.fa}.paf && paftools.js call -s ${FASTA%.fa} -L 100 -f $REF_GENOME ${FASTA%.fa}.paf > ${FASTA%.fa}.vcf 2>>${FASTA%.fa}.paftools.log;
    bgzip -f ${FASTA%.fa}.vcf;
    bcftools index -f ${FASTA%.fa}.vcf.gz;
done
cd "${OUTPUT_DIR}"
# Count samples
SAMPLE_COUNT=$(ls ${BASE_DIR}/${LINEAGE}/*.vcf.gz | wc -l);
echo "There are ${SAMPLE_COUNT} genomes to select from"
if [[ ${SAMPLE_COUNT} -eq 1 ]]; then
    cp ${BASE_DIR}/${LINEAGE}/*.vcf.gz ${OUTPUT_DIR}/${LINEAGE}_merged.vcf.gz;
else
    bcftools merge -o "${OUTPUT_DIR}/${LINEAGE}_merged.vcf.gz" -O z -0 ${BASE_DIR}/${LINEAGE}/*.vcf.gz;
fi
vcftools --gzvcf "${OUTPUT_DIR}/${LINEAGE}_merged.vcf.gz" --out "${OUTPUT_DIR}/${LINEAGE}_merged" --site-pi;
vcftools --gzvcf "${OUTPUT_DIR}/${LINEAGE}_merged.vcf.gz" --out "${OUTPUT_DIR}/${LINEAGE}_merged" --freq;
echo "Done processing lineage ${LINEAGE}, removing .paf, .paftools.log, .vcf.gz and .vcf.gz.csi files"
find "${BASE_DIR}/${LINEAGE}" -type f \( -name "*.paf" -o -name "*.paftools.log" -o -name "*.vcf.gz" -o -name "*.vcf.gz.csi" \) -delete
