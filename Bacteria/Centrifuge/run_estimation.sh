INDEX_LOC=$1
OUTPUT_LOC=$2
THREADS=$3
READS_PREFIX=$4

centrifuge -p $THREADS -x $INDEX_LOC -1 "${READS_PREFIX}_1.fq" -2 "${READS_PREFIX}_2.fq" -S $OUTPUT_LOC/output.sam --report-file $OUTPUT_LOC/output.report
# Determine number of unmapped reads
cut -f2 $OUTPUT_LOC/output.sam | awk '$1 == "unclassified" {count++} END {printf"%d, ", count}; done | sed 's/,$/\n/' > $OUTPUT_LOC/unclassified