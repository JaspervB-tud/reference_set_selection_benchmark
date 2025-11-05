## This script assumes that the GISAID data (sequences.fasta, metadata.tsv) is already downloaded in the corresponding root folder.
## Furthermore, we assume that the VLQ pipeline was cloned in root/scripts/virus/wastewater_analysis

NCONTENT=0.001
NGENOMES=1000
STARTDATE="2021-01-01"
ENDDATE="2021-03-31"
SEQUENCES="root/sequences.fasta"
METADATA="root/metadata.tsv"

for reference_set in "Connecticut,Connecticut,state"  "USA,USA,country" "North America,North_America,global"; do 
    IFS=',' read -r location folder_name location_type <<< "$reference_set"
    echo "Now processing ${location}: ${location_type}"
    mkdir -p root/genomes/${location_type}

    if [[ "${location_type}" == "global" ]]; then
        python3 -u root/scripts/virus/wastewater_analysis/pipeline/preprocess_references.py -m $METADATA -f $SEQUENCES -o "root/genomes/${location_type}" --startdate $STARTDATE --enddate $ENDDATE --max_N_content $NCONTENT -k $NGENOMES
    elif [[ "${location_type}" == "country" ]]; then
        python3 -u root/scripts/virus/wastewater_analysis/pipeline/preprocess_references.py -m $METADATA -f $SEQUENCES -o "root/genomes/${location_type}" --startdate $STARTDATE --enddate $ENDDATE --country $location --max_N_content $NCONTENT -k $NGENOMES
    elif [[ "${location_type}" == "state" ]]; then
        python3 -u root/scripts/virus/wastewater_analysis/pipeline/preprocess_references.py -m $METADATA -f $SEQUENCES -o "root/genomes/${location_type}" --startdate $STARTDATE --enddate $ENDDATE --state $location --max_N_content $NCONTENT -k $NGENOMES
    fi
done