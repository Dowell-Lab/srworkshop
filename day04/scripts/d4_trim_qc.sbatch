#!/bin/bash
#SBATCH --job-name=<JOB_NAME>		       	# Job name
#SBATCH --mail-type=END                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOUR_EMAIL>           	# Where to send mail
#SBATCH --nodes=1                           # Number of nodes requested
#SBATCH --ntasks=2                          # Number of CPUs (processor cores/tasks)
#SBATCH --mem=1gb                           # Memory limit
#SBATCH --time=00:10:00                     # Time limit hrs:min:sec
#SBATCH --partition=short                 # Partition/queue requested on server
#SBATCH --output=/scratch/Users/<USERNAME>/workshop-day4/eofiles/%x.%j.out  
#SBATCH --error=/scratch/Users/<USERNAME>/workshop-day4/eofiles/%x.%j.err

################## Display job context ################################
## This script takes fastq files as input to fastqc and trimmomatic
## Trimmed reads will then be passed on to the mapping step
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

################## Assign variables ###################################
fastq=/scratch/Shares/public/sread/data_files/day4/fastq
outdir=/scratch/Users/<USERNAME>/workshop-day4
fastqc="$outdir"/fastqc
trim="$outdir"/trimmomatic

filebase=chr21Eric_repA.RNA  # Prefix of the files to analyze
adapter_file=/data/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa # Adapter sequences

mkdir -p $outdir
mkdir -p $fastqc
mkdir -p $trim

################## Load modules/software ##############################
module load <fastqc module> # Check actual module name on command line

# Trimmomatic is a program that we run the executable manually with Java
# Set this executable path as a variable
trimmomatic_exec=/data/software/Trimmomatic-0.39/trimmomatic-0.39.jar

################## Run the pipeline ###################################
##1: Run fastqc on the samples (here run on example file "$filebase".end1.fastq)
fastqc -o $fastqc "$fastq"/"$filebase".end1.fastq
fastqc -o $fastqc "$fastq"/"$filebase".end2.fastq

##2: Trim FASTQ Files
java -jar $trimmomatic_exec PE \
 -threads $SLURM_NTASKS \
 -phred33 \
 -trimlog "$trim"/trimlog \
 "$fastq"/"$filebase".end1.fastq "$fastq"/"$filebase".end2.fastq \
 "$trim"/"$filebase".end1.trimmed.fastq "$trim"/"$filebase".end1.unpaired.fastq \
 "$trim"/"$filebase".end2.trimmed.fastq "$trim"/"$filebase".end2.unpaired.fastq \
 ILLUMINACLIP:"$adapter_file":2:30:10 \
 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

##3: Check Post-Trimming QC stats
fastqc -o "$fastqc" "$trim"/*.trimmed.fastq

echo Job finished at `date +"%T %a %d %b %Y"`
