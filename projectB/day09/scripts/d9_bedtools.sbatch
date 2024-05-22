#!/bin/bash
#SBATCH --job-name=d9_bedtools                              # Job name
#SBATCH --mail-type=NONE                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<email>                        # Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:10:00                                 # Time limit hrs:min:sec
#SBATCH --mem=2gb                                       # Memory limit
#SBATCH --output=/scratch/Users/<username>/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<username>/eofiles/%x_%j.err

################## JOB INFO ########################################

printf "\nBAM Directory: $BAMS"
printf "\nOutput Directory: $OUT"
printf "\nOutput Rootname: $FILENAME"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

################### LOAD MODULES ###################################

module load bedtools/2.25.0

################## SET VARIABLES ###################################

chrom_sizes='/scratch/Shares/public/genomes/hisatfiles/hg38/hg38.chrom.sizes'
chr21_genes=''
BACH1_peaks=''
results='/scratch/Users/<username>/day9/bedtools'

mkdir -p "$results"

################## SCRIPT ##########################################

# Write your bedtools/file manipulation commands here


