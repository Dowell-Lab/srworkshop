#!/bin/bash
#SBATCH --output=/scratch/Users/<username>/homework-day3/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<username>/homework-day3/eofiles/%x_%j.err
#SBATCH --mail-user=<email address>
#SBATCH --mail-type=END
#SBATCH --time=00:20:00
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1Gb

# Load our module in the proper version
module load fastqc/0.11.5

date

# Change to our data directory
cd /scratch/Users/<username>/homework-day3/
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR792/SRR792473/SRR792473_1.fastq.gz

date

# Run fastqc on the file we downloaded
fastqc SRR792473_1.fastq.gz

date

# Remove the fastq file
rm SRR792473_1.fastq.gz

# Back up data to home
rsync -r ../homework-day3 ~/