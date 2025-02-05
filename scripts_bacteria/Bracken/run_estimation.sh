INDEX_LOC=$1
OUTPUT_LOC=$2
THREADS=$3
READS_PREFIX=$4

kraken2 -db $INDEX_LOC -threads $THREADS -report $OUTPUT_LOC/output.kreport --paired "${READS_PREFIX}_1.fq" "${READS_PREFIX}_2.fq" > $OUTPUT_LOC/output.kraken
bracken -d $INDEX_LOC -i $OUTPUT_LOC/output.kreport -o $OUTPUT_LOC/bracken_output.bracken -r 150 -l S