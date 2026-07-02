#!/bin/bash
#SBATCH --job-name=<JOB_NAME>                         # Job Name
#SBATCH --mail-type=ALL                               # Mail (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOU@EMAIL.COM>                   # Your email address
#SBATCH --nodes=1                                     # Number of nodes requested, should always be 1
#SBATCH --ntasks=1                                    # Number of CPUs (processor cores/tasks)
#SBATCH --mem=256mb                                   # Memory limit, formatted 256mb, 8gb, etc.
#SBATCH --time=00:10:00                               # Time limit hrs:min:sec
#SBATCH --partition=short                             # Partition/queue requested on server
#SBATCH --output=/scratch/Users/<YOURUSERNAME>/workshop-day3/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<YOURUSERNAME>/workshop-day3/eofiles/%x_%j.err

# Original rsync command to copy the file that doesn't exist
# Commented out so that it won't actually run, but if this file did exist, the command would work
#rsync /scratch/Shares/public/sread2025/data_files/day3/SRR062641.filt.fastq.gz /scratch/Users/<username>/workshop-day3/data/

# Download our data
cd /scratch/Users/<username>/workshop-day3/data/
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/sequence_read/SRR062641.filt.fastq.gz
