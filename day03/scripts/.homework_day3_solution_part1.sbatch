#!/bin/bash
#SBATCH --job-name=<JOB_NAME>                         # Job Name
#SBATCH --mail-type=ALL                               # Mail (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOU@EMAIL.COM>                   # Your email address
#SBATCH --nodes=1                                     # Number of nodes requested, should always be 1
#SBATCH --ntasks=1                                    # Number of CPUs (processor cores/tasks)
#SBATCH --mem=256mb                                   # Memory limit, formatted 256mb, 8gb, etc.
#SBATCH --time=00:10:00                               # Time limit hrs:min:sec
#SBATCH --partition=short                             # Partition/queue requested on server
#SBATCH --output=/scratch/Users/<YOURUSERNAME>/homework-day3/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<YOURUSERNAME>/homework-day3/eofiles/%x_%j.err

date
sleep 120
date
