INPUT_GENOMES=$1
OUTPUT_LOC=$2
TXDMP=$3

bwa index -b 50000000 $INPUT_GENOMES -p $OUTPUT_LOC
dudesdb -m 'av' -f $INPUT_GENOMES -n $TXDMP/nodes.dmp -a $TXDMP/names.dmp -g $TXDMP/nucl.accession2taxid -o $OUTPUT_LOC

