#!/bin/bash
#SBATCH --job-name=<JOB_NAME>		       	# Job name
#SBATCH --mail-type=ALL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOUR_EMAIL>           	# Where to send mail
#SBATCH --nodes=1                               # Number of nodes requested
#SBATCH --ntasks=8                              # Number of CPUs (processor cores/tasks)
#SBATCH --mem=8gb                               # Memory limit
#SBATCH --time=01:30:00                        	# Time limit hrs:min:sec
#SBATCH --partition=short                       # Partition/queue requested on server
#SBATCH --output=/scratch/Users/<USERNAME>/eofiles/%x.%j.out  
#SBATCH --error=/scratch/Users/<USERNAME>/eofiles/%x.%j.err

##########################################################################################
################## DISPLAY JOB CONTEXT ###################################################

echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

################ ASSIGNS PATH VARIABLES ##################################################
## the fastq files will be used as input to fastqc and trimmomatic
## trimmed reads will then be passed on to the mapping step

INDIR=/scratch/Shares/public/sread2022/data_files/day4
FASTQ=${INDIR}/fastq

OUTDIR=/scratch/Users/<USERNAME>/day4
FASTQC=${OUTDIR}/fastqc
TRIM=${OUTDIR}/trimmomatic

FILENAME=chr21Eric_repA

mkdir -p ${OUTDIR}
mkdir -p ${FASTQC}
mkdir -p ${TRIM}

##########################################################################################
############### LOAD REQUIRED MODULES ####################################################
module load fastqc/0.11.5
module load trimmomatic/0.36

##########################################################################################
################## RUN PIPELINE ##########################################################
##1: Run fastqc on the samples (here run on example file ${FILENAME}.RNA.end1.fastq)
fastqc ${FASTQ}/${FILENAME}.RNA.end1.fastq -o ${FASTQC}
fastqc ${FASTQ}/${FILENAME}.RNA.end2.fastq -o ${FASTQC}

##2: Trim FASTQ Files                                                                                                                                                                                                                                                
java -jar /opt/trimmomatic/0.36/trimmomatic-0.36.jar PE \
 -threads 8 \
 -phred33 \
 -trimlog ${TRIM}/trimlog \
 ${FASTQ}/${FILENAME}.RNA.end1.fastq ${FASTQ}/${FILENAME}.RNA.end2.fastq \
 ${TRIM}/${FILENAME}.RNA.end1.trimmed.fastq ${TRIM}/${FILENAME}.RNA.end1.unpaired.fastq \
 ${TRIM}/${FILENAME}.RNA.end2.trimmed.fastq ${TRIM}/${FILENAME}.RNA.end2.unpaired.fastq \
 ILLUMINACLIP:/opt/trimmomatic/0.36/adapters/TruSeq3-PE.fa:2:30:10 \
 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

##3: Check Post-Trimming QC stats
fastqc ${TRIM}/*.trimmed.fastq -o ${FASTQC}

echo Job finished at `date +"%T %a %d %b %Y"`
