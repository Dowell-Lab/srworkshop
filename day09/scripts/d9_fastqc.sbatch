/scratch/Users/jessicatwes/#!/bin/bash
#SBATCH --job-name=fastqc                               # Job name
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<email>                      	# Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:05:00                                 # Time limit hrs:min:sec
#SBATCH --mem=50mb                                     # Memory limit
#SBATCH --output=/scratch/Users/<<YourUsername>/eofiles/d9_fastqc.out
#SBATCH --error=/scratch/Users/<YourUsername>/eofiles/d9_fastqc.err

################ SET REQUIRED VARIABLES ##################################################
## the fastq files will be used as input to fastqc. 
## output will be a fastqc file used to assess quality

INDIR=/scratch/Users/<YourUsername>/day9/fastq
OUTDIR=/scratch/Users/<YourUsername>/day9/qc/fastqc

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
#mkdir -p ${OUTDIR}

for FILENAME in `ls $INDIR`; do

    fastqc \
	${INDIR}/${FILENAME} \
	-o ${OUTDIR}
        
done
                                                                        
echo Job finished at `date +"%T %a %d %b %Y"`
