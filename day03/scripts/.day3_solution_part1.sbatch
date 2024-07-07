#!/bin/bash
#SBATCH --output=/scratch/Users/<username>/workshop-day3/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<username>/workshop-day3/eofiles/%x_%j.err
#SBATCH --mail-user=<email address>
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=256Mb

# Make the directory if it doesn't exist, -p won't complain if it does
mkdir -p /scratch/Users/<username>/workshop-day3/data
# Change into the directory
cd /scratch/Users/<username>/workshop-day3/data
# Download our data
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/sequence_read/SRR062641.filt.fastq.gz
