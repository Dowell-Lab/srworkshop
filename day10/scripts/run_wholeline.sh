#!/bin/bash

indir=/scratch/Users/ispe1418/disease_datasets/ptsd/scripts/
wholefileoflines=${indir}ptsd_download_links
outdir=/scratch/Users/ispe1418/disease_datasets/ptsd/fastq/

mkdir -p "$outdir"

echo "Getting the SRR Accession Numbers"

nlines=$(wc -l < "$wholefileoflines")
echo "Total lines: $nlines"

# Submit a single array job, capped at N concurrent tasks
MAX_CONCURRENT=5   # adjust as needed

sbatch --array=1-"$nlines"%${MAX_CONCURRENT} \
       --export=outdir="$outdir",infile="$wholefileoflines" \
       acommandsbatch.sbatch



echo "DONE!"
