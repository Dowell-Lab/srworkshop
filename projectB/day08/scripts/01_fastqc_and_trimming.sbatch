#!/bin/bash                              
#SBATCH --job-name=<JOB ID>                                  # Job name
#SBATCH --mail-type=ALL                                      # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<E-MAIL>                                 # Where to send mail
#SBATCH --nodes=1                                            # Numbers of nodes
#SBATCH --ntasks=1                                           # Number of tasks
#SBATCH --cpus-per-task=2                                    # Number CPU
#SBATCH --time=00:30:00 # Time limit                         # Time limit hrs:min:sec
#SBATCH --partition=short                                  # Partition/queue requested on server
#SBATCH --mem=500mb                                          # Memory limit
#SBATCH --output=/scratch/Users/<username>/day8/eofiles/%x_%j.out # Where to save out files             
#SBATCH --error=/scratch/Users/<username>/day8/eofiles/%x_%j.err  # Where to save error files

echo "BEGIN: Preprocessing................."

################################################################                                                                                                                      
##Load required modules 
################################################################

module load fastqc/0.11.5
module load trimmomatic/0.36
module load java/1.8.0_101

################################################################
##Initialize folders 
################################################################
wd=/scratch/Users/<YourUsername>/day8
qc=${wd}/qc
trim=${wd}/trim
fastq=/scratch/Shares/public/sread2025/data_files/day8b/fastq

#make folders that do not exist
mkdir -p ${qc} ${trim} 

################################################################
##Job information
################################################################
printf "\nFASTQ Directory: ${fastq}"
printf "\nQC Directory: ${qc}"
printf "\nTrimmed FASTQ Directory: ${trim}"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou have requested $SLURM_CPUS_ON_NODE core(s).\n"

################################################################
##Run script                                                     
################################################################
#create qc subfolders to store intermediate files
mkdir -p ${qc}/trim_log
mkdir -p ${qc}/fastqc

#here were are using loops 
#CHALLENGE: what do these commands do (ls, tr, cut)? 
for sample in `ls ${fastq} | tr '.' '\t' | cut -f 1`; 
do
    #1: trim fastq file using trimmomatic 
    java -jar /opt/trimmomatic/0.36/trimmomatic-0.36.jar SE \
	-threads 2 \
	-phred33 \
	-trimlog \
	${qc}/trim_log/${sample}.trimlog \
	${fastq}/${sample}.fastq.gz \
	${trim}/${sample}.trimmed.fastq.gz \
	ILLUMINACLIP:/opt/trimmomatic/0.36/adapters/TruSeq3-SE.fa:2:30:10 \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    #2: check quality of the fastq files after trimming 
    #CHALLENGE: Compare them to the raw input 
    fastqc --threads=4 ${trim}/${sample}.trimmed.fastq.gz -o ${qc}/fastqc

done

echo "DONE: Preprocessing................."
echo Job finished at `date +"%T %a %d %b %Y"`
