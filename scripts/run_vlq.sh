#!/bin/sh
#SBATCH --qos=medium         # Request Quality of Service. Default is 'short' (maximum run time: 4 hours)
#SBATCH --time=8:00:00      # Request run time (wall-clock). Default is 1 minute
#SBATCH --ntasks=1          # Request number of parallel tasks per job. Default is 1
#SBATCH --cpus-per-task=20   # Request number of CPUs (threads) per task. Default is 1 (note: CPUs are always allocated to jobs per 2).
#SBATCH --mem=32GB          # Request memory (MB) per node. Default is 1024MB (1GB). For multiple tasks, specify --mem-per-cpu instead
#SBATCH --mail-type=ALL     # Set mail type to 'END' to receive a mail when the job finishes. 
#SBATCH --output=sbatch_output/prepare_genomes_america_global%j.out # Set name of output log. %j is the Slurm jobId
#SBATCH --error=sbatch_output/prepare_genomes_america_global__%j.err # Set name of error log. %j is the Slurm jobId

##/usr/bin/scontrol show job -d "$SLURM_JOB_ID"  # check sbatch directives are working
LINEAGE=$1
REF_GENOME=$2
BASE_DIR=$3
OUTPUT_DIR=$4
PATH_TO_VLQ=$5

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
# Run selection script
python3 "${PATH_TO_VLQ}/pipeline/select_samples.py" \
    -m "${BASE_DIR}/${LINEAGE}/metadata.tsv" \
    -f "${BASE_DIR}/${LINEAGE}/all_genomes.fasta" \
    -o ${OUTPUT_DIR} \
    --vcf ${OUTPUT_DIR}/*_merged.vcf.gz \
    --freq ${OUTPUT_DIR}/*_merged.frq \
    --max_per_lineage 1000
