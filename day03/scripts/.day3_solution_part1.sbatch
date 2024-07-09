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

# Original rsync command to copy the file that doesn't exist
# Commented out so that it won't actually run, but if this file did exist, the command would work
#rsync /scratch/Shares/public/sread2024/data_files/day3/SRR062641.filt.fastq.gz /scratch/Users/<username>/workshop-day3/data/

# Download our data
cd /scratch/Users/<username>/workshop-day3/data/
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/sequence_read/SRR062641.filt.fastq.gz
