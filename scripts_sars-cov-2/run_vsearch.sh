INPUT_GENOMES=$1
OUTPUT_FILE=$2
THREADS=$3
THRESHOLD=$4

# iddef 0 corresponds to CD-HIT similarity scores
vsearch --cluster_fast $INPUT_GENOMES --centroids $OUTPUT_FILE --id ${THRESHOLD} --iddef 0 --qmask none --threads $THREADS