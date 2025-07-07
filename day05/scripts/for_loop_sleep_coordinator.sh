#!/bin/bash

indir=/scratch/Shares/public/sread2025/data_files/day5/fastq/for_loops_fastq

for filepath in "$indir"/*.fastq
do 
echo "$filepath"
sbatch --export=filepath="$filepath" for_loop_sleep_parallel.sbatch
done