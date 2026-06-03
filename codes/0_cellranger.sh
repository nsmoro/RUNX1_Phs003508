#!/bin/bash

module load cellranger/7.1.0

samplelist="./"

for f in "$samplelist"/*; 
do
file=$(basename "$f")
#echo "$file"
base_name="${file%%_*}"
y=$base_name
#echo $y
ids="${y}_counts"
#echo $ids
cellranger count --id=$ids --fastqs=./ --sample=$y --transcriptome=/home/refdata-gex-GRCh38-2024-A/
done

