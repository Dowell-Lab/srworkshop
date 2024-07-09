#!/bin/bash
#SBATCH --output=/scratch/Users/<username>/workshop-day3/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<username>/workshop-day3/eofiles/%x_%j.err
#SBATCH --mail-user=<email address>
#SBATCH --mail-type=END
#SBATCH --time=00:10:00
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=256Mb

# Load our module in the proper version
module load fastqc/0.11.5

# Change to our data directory
cd /scratch/Users/<username>/workshop-day3/data

# Run fastqc on the file we downloaded
fastqc SRR062641.filt.fastq.gz

# Alternately, instead of the last two steps, we can run fastqc on the absolute path from any directory
#fastqc /scratch/Users/<username>/workshop-day3/data/SRR062641.filt.fastq.gz

# Optional command to direct output to results folder
#fastqc -o /scratch/Users/<username>/workshop-day3/results/ /scratch/Users/<username>/workshop-day3/data/SRR062641.filt.fastq.gz
