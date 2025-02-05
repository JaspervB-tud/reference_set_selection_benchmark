INPUT_GENOMES=$1
OUTPUT_FILE=$2
THREADS=$3
THRESHOLD=$4

# Parameters are largely taken from original publication (table 3): https://academic.oup.com/gpb/article/17/5/496/7229746
gclust -minlen 41 -threads $THREADS -chunk 400 -ext 1 -sparse 4 -nuc -loadall -rebuild -memiden $THRESHOLD $INPUT_GENOMES > $OUTPUT_FILE