#!/bin/bash
#SBATCH --job-name=d9_bedtools                          # Job name
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=meaghan.courvan@colorado.edu        # Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:10:00                                 # Time limit hrs:min:sec
#SBATCH --mem=2gb                                       # Memory limit
#SBATCH --output=/scratch/Shares/dowell/meco9877/andrysik2017/%x_%j.out
#SBATCH --error=/scratch/Shares/dowell/meco9877/andrysik2017/%x_%j.err


################## USER INPUT ######################################

chrom_sizes='/scratch/Shares/public/genomes/hisatfiles/hg38/hg38.chrom.sizes'
ANNOT='/scratch/Shares/dowell/meco9877/andrysik2017/rnaseq/bedtools_demo/hg38_refseq_cleaned.sorted.bed'
HCT116_DE='/scratch/Shares/dowell/meco9877/andrysik2017/rnaseq/output/counts/de_output/andrysik2017_hct116_nutlin_dmso.bed'
CHIP_DIR='/scratch/Shares/dowell/rutendo/projects/sread/sread2024/day8/cooking_show/macs/w_ctrl'
HCT116_DMSO_PEAKS=${CHIP_DIR}'/HCT116_DMSO_peaks_clean.narrowPeak'
HCT116_NUTLIN_PEAKS=${CHIP_DIR}'/HCT116_NUTLIN_peaks_clean.narrowPeak'
results='/your/results/directory'
TSS_ANNOT=${results}/

mkdir -p "$results"
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
