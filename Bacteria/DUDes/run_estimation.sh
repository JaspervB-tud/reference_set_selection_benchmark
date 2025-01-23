INDEX_LOC=$1
OUTPUT_LOC=$2
THREADS=$3
READS_PREFIX=$4

bwa mem -t $THREADS -v 1 $INDEX_LOC "${READS_PREFIX}_1.fq" "${READS_PREFIX}_2.fq" > $OUTPUT_LOC/output.sam
# Filter unaligned reads
samtools view -F 4 -h $OUTPUT_LOC/output.sam > $OUTPUT_LOC/output_filtered.sam
# Calculate unaligned reads
samtools view -c -f 4 $OUTPUT_LOC/output.sam > $OUTPUT_LOC/num_unaligned.txt
dudes -s $OUTPUT_LOC/output_filtered.sam -d $INDEX_LOC -o $OUTPUT_LOC/output.report -l species -t $THREADS