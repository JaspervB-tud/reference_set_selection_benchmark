INPUT_GENOMES=$1
OUTPUT_LOC=$2
THREADS=$3

makedir -p $OUTPUT_LOC
kraken2-build --add-to-library $INPUT_GENOMES --db $OUTPUT_LOC --threads $THREADS
kraken2-build --build --db $OUTPUT_LOC --threads $THREADS
bracken-build -d $OUTPUT_LOC -t $THREADS -l 150
kraken2-build --clean --db $OUTPUT_LOC --threads $THREADS