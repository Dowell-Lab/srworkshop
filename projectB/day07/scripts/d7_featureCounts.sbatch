#!/bin/bash

#SBATCH --job-name=<NAME OF JOB>				 # Job name
#SBATCH --mail-type=ALL				 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOUR E-MAIL ADDRESS>		 # Where to send mail
#SBATCH --nodes=1					 # Number of cores job will run on
#SBATCH --ntasks=4					 # Number of CPU (processers, tasks)
#SBATCH --time=00:10:00					 # Time limit hrs:min:sec
#SBATCH --partition short				 # Job queue
#SBATCH --mem=4mb					 # Memory limit
#SBATCH --output=/YOUR/EOFILES/PATH/%x_%j.out
#SBATCH --error=/YOUR/EOFILES/PATH/%x_%j.err

################### SET VARIABLES ######################################

FEATURECOUNTS=/PATH/TO/YOUR/d7_featureCounts.R

########################################################################
################## PRINT JOB INFO ######################################

printf "Sample ID: $ROOTNAME"
printf "\nDirectory: $PROJECT"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

#######################################################################
################## LOAD MODULES #######################################
module load R/4.3.1

Rscript $FEATURECOUNTS
