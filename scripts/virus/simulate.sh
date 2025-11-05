mkdir -p root/samples
metadata="root/metadata.tsv"
sequences="root/sequences.fasta"
state="Connecticut"

mkdir -p root/samples
for seed in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do \
    python -u root/scripts/virus/create_benchmarks_with_seed.py --voc_perc 1,10,20,30,40,50,60,70,80,90,100 -m $metadata -fr $sequences -fv root/B.1.1.7_sequence.fasta -o root/samples/$seed --total_cov 100 -s "North America / USA / ${state}" -d 2021-04-30 --seed $seed

    for file in root/samples/$seed/*; do
        if ! [[ "$file"  =~ ^(root/samples/$seed/wwsim_|root/samples/$seed/metadata.tsv|root/samples/$seed/sequences.fasta).*$ ]]; 
        then 
            echo "$file is removed"
            rm $file
        fi
    done
done