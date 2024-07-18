#!/bin/bash
#SBATCH --job-name=d9_bedtools                          # Job name
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<EMAIL>        # Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:10:00                                 # Time limit hrs:min:sec
#SBATCH --mem=2gb                                       # Memory limit
#SBATCH -p short
#SBATCH --output=/scratch/Users/<username>/day9/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<username>/day9/eofiles/%x_%j.err


################## USER INPUT ######################################

chrom_sizes='/scratch/Shares/public/genomes/hisatfiles/hg38/hg38.chrom.sizes'
ANNOT='</path/to/annotation/file.bed>'
CHIP_DIR='/scratch/Shares/public/sread2024/cookingShow/day8b/whole_genome/macs/'
HCT116_DMSO_PEAKS=${CHIP_DIR}/HCT116_DMSO_peaks_clean.narrowPeak
HCT116_NUTLIN_PEAKS=${CHIP_DIR}/HCT116_NUTLIN_peaks_clean.narrowPeak
results='/scratch/Users/<username>/day9/bedtools_results/'
TSS_ANNOT=${results}/hg38_refseq_tssregion.bed

mkdir -p ${results}
################## JOB INFO ########################################

printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

################### LOAD MODULES ###################################

module load bedtools/2.28.0


################## SCRIPT ##########################################

# Write your bedtools/file manipulation commands here


# --- Which p53 peaks are associated with genes (not intergenic)


# --- Which p53 peeks are entirely within genes 


# --- Which p53 peaks are intergenic 


# --- Which genes overlap with p53 peaks? 


# --- Which genes have p53 ChIP sites in their promoters



###---- Now do for Nutlin treated HCT116 Samples ----###




###---- Jacquard Calculations ----###
