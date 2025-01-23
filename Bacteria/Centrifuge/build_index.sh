INPUT_GENOMES=$1
OUTPUT_LOC=$2
TXDMP=$3
THREADS=$4

centrifuge-build -p $THREADS --conversion-table $TXDMP/nucl.accession2taxid --taxonomy-tree $TXDMP/nodes.dmp  --name-table $TXDMP/names.dmp $INPUT_GENOMES $OUTPUT_LOC