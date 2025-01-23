
INDEX_LOC=$1
OUTPUT_LOC=$2
READS_PREFIX=$3
METADATA_LOC=$4

# Run kallisto
kallisto quant -t 16 -b 0 -i $INDEX_LOC -o $OUTPUT_LOC "${READS_PREFIX}_1.fq" "${READS_PREFIX}_2.fq"
# Run VLQ
output_abundances.py -m 0.1 -o $OUTPUT_LOC/predictions.tsv -metadata $METADATA_LOC $OUTPUT_LOC/abundances.tsv