INPUT_GENOMES=$1
OUTPUT_FILE=$2
THREADS=$3

mash sketch -s 1000 -S 123456 -k 21 -p 8 -o $OUTPUT_FILE $INPUT_GENOMES