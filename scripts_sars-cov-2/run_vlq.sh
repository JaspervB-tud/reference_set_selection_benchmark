REF_GENOME=$1
FASTA=$2
OUTPUT_DIR=$3
LINEAGE=$4
BASE_DIR=$5
THREADS=$6

minimap2 -c -x asm20 --end-bonus 100 -t $THREADS --cs $REF_GENOME $FASTA 2>${FASTA%.fa}.paftools.log | sort -k6,6 -k8,8n > ${FASTA%.fa}.paf && paftools.js call -s ${FASTA%.fa} -L 100 -f $REF_GENOME ${FASTA%.fa}.paf > ${FASTA%.fa}.vcf 2>>${FASTA%.fa}.paftools.log;
bgzip -f ${FASTA%.fa}.vcf;
bcftools index -f ${FASTA%.fa}.vcf.gz;
bcftools merge -o "${OUTPUT_DIR}/${LINEAGE}_merged.vcf.gz" -O z -0 ${BASE_DIR}/${LINEAGE}/*.vcf.gz;
vcftools --gzvcf "${OUTPUT_DIR}/${LINEAGE}_merged.vcf.gz" --out "${OUTPUT_DIR}/${LINEAGE}_merged" --site-pi;
vcftools --gzvcf "${OUTPUT_DIR}/${LINEAGE}_merged.vcf.gz" --out "${OUTPUT_DIR}/${LINEAGE}_merged" --freq;
find "${BASE_DIR}/${LINEAGE}" -type f \( -name "*.paf" -o -name "*.paftools.log" -o -name "*.vcf.gz" -o -name "*.vcf.gz.csi" \) -delete
python3 select_samples.py -m metadata.tsv -f sequences.fasta -o ${OUTPUT_DIR} --vcf ${OUTPUT_DIR}/*_merged.vcf.gz --freq ${OUTPUT_DIR}/*_merged.frq
