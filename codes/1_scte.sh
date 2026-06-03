#!/bin/bash

module load samtools/1.15.1

#SAMTOOLS

base_dir="/home/files"

for sub_dir in "$base_dir"/*_counts; do
  mkdir -p "$sub_dir/outs/filtered"

  for bam_file in "$sub_dir/outs"/*.bam; do
    filename=$(basename "$bam_file" .bam)
    output_file="$sub_dir/outs/filtered/${filename}_filtered.bam"

    samtools view "$bam_file" -h | \
    awk '/^@/ || /CB:/' | \
    samtools view -h -b > "$output_file"
  done
done

#scTE - using filtered bam file

base_dir="/home/files"

cd "$base_dir"

for sub_dir in "$base_dir"/*_counts; do
    sub_dir_name=$(basename "$sub_dir")
    prefix=${sub_dir_name%%_*}
    for outs_dir in "$sub_dir"/outs; do
        for filt_dir in "$outs_dir"/filtered; do
            for bam_file in "$filt_dir"/*.bam; do
                output_dir="$filt_dir/${prefix}_out"
                mkdir -p "$output_dir"
                scTE -i "$bam_file" -o "$output_dir" -x /home/mm10.exclusive.idx --hdf5 True -CB CB -UMI UB
            done
        done
    done
done

#scTE - using unfiltered bam file

base_dir="/home/files"

for sub_dir in "$base_dir"/*_counts; do
    sub_dir_name=$(basename "$sub_dir")
    prefix=${sub_dir_name%%_*}
    outs_dir="$sub_dir/outs"
    bam_file="$outs_dir/possorted_genome_bam.bam"
    barcode_file="$outs_dir/raw_feature_bc_matrix/barcodes.tsv.gz"
    output_dir="$outs_dir/${prefix}_out"

    mkdir -p "$output_dir"
    gunzip -c "$barcode_file" > "$output_dir/barcodes.tsv"

    bsub -J scTE_${prefix} -n 8 -R "rusage[mem=8000]" -oo "$output_dir/${prefix}.log" \
    "cd $output_dir && scTE -i $bam_file -o $prefix -x /home/mm10.exclusive.idx --hdf5 True --expect-cells 99999"
done

