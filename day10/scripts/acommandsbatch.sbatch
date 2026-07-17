#!/bin/bash

#SBATCH --job-name=get_data           # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=50mb
#SBATCH --output=/scratch/Users/ispe1418/disease_datasets/ptsd/eando/%x_%A_%a.out
#SBATCH --error=/scratch/Users/ispe1418/disease_datasets/ptsd/eando/%x_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ispe1418@colorado.edu

# This script is meant to be run as a job array.
# Each array task handles one line (one wget) from $infile.

echo "Run on: $(hostname)"
echo "Run from: $(pwd)"
echo "Script:  $0"
date

echo "Array job ID:  $SLURM_ARRAY_JOB_ID"
echo "Task ID:       $SLURM_ARRAY_TASK_ID"
echo "Requested CPUs: $SLURM_CPUS_ON_NODE"

# Variables coming from --export in sbatch: outdir, infile
echo "outdir:  $outdir"
echo "infile:  $infile"

mkdir -p "$outdir"

whichline="$SLURM_ARRAY_TASK_ID"

# Get the wget line for this task
line_content=$(sed -n "${whichline}p" "$infile")
if [ -z "$line_content" ]; then
    echo "No content for line $whichline in $infile, exiting."
    exit 0
fi

# Extract URL and filename
url=${line_content##* }
fname=${url##*/}
fullpath="${outdir}/${fname}"
#fullpath="${outdir%/}/${fname}" this is if you don't like the double slash in your path

echo "URL:     $url"
echo "Target:  $fullpath"

# If file already exists, skip
if [ -f "$fullpath" ]; then
    echo "File already exists, skipping download: $fullpath"
    exit 0
fi

# Add a rate limit (adjust as needed)
#this just relaces wget in the line with wget --limit-rate=200k because the host server was locking us out for downloading to much
cmd="${line_content/wget /wget --limit-rate=200k }"

echo "Running: $cmd"
cd "$outdir"
eval "$cmd" #this runs to wget command in my orginal line_content

echo "DONE"
