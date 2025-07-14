#!/bin/bash
#SBATCH --job-name=gene
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Mary.A.Allen@colorado.edu                     # Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:05:00                                 # Time limit hrs:min:sec
#SBATCH --partition=short                               # Partition/queue requested on server    
#SBATCH --mem=3mb                                       # Memory limit
#SBATCH --output=/scratch/Users/maallen3/day6/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/maallen3/day6/eofiles/%x_%j.err

################################################################                                                                                                                      
##Load required modules 
################################################################
module load R/4.3.1

################################################################ 
##Set paths
################################################################
rscript_path=/path/to/step1_genes_R_submit_aws.R

################################################################ 
##Job information
################################################################ 
printf "\nR code being run: ${rscript_path}"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou have requested $SLURM_CPUS_ON_NODE core(s).\n"

################################################################                  
## Run Script 
################################################################
#using the Rscript command to run the script
Rscript ${rscript_path} $agenename

echo Job finished at `date +"%T %a %d %b %Y"`
