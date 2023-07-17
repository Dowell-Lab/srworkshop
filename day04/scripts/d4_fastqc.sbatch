#!/bin/bash
#SBATCH --job-name=<JOB_NAME>                           # Job name
#SBATCH --mail-type=FAIL                                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOUR_EMAIL>              		# Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:15:00                                 # Time limit hrs:min:sec
#SBATCH --partition=short                               # Partition/queue requested on server    
#SBATCH --mem=2gb                                       # Memory limit
#SBATCH --output=/scratch/Users/<USERNAME>/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<USERNAME>/eofiles/%x_%j.err

################ SET REQUIRED VARIABLES ##################################################
## the fastq files will be used as input to fastqc. 
## output will be a fastqc file used to assess quality

INDIR=/scratch/Workshop/sread/data_files/day4/fastq
OUTDIR=/scratch/Users/<USERNAME>/fastqc

FILENAME=chr21Ethan_repA.RNA.end1

##########################################################################################
############### LOAD REQUIRED MODULES ####################################################

module load fastqc/0.11.5

##########################################################################################
################## PRINT JOB INFO ########################################################

printf "\nfastq Directory: $INDIR"
printf "\nOutput Directory: $OUTDIR"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

##########################################################################################
################## RUN JOB ###############################################################
mkdir -p ${OUTDIR}

fastqc \
 ${INDIR}/${FILENAME}.fastq \
 -o ${OUTDIR}
                                                                                
echo Job finished at `date +"%T %a %d %b %Y"`
