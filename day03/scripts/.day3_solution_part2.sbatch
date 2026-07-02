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
