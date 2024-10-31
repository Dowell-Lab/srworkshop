#!/bin/bash
#SBATCH --job-name=<JOB ID>                                    # Job name
#SBATCH --mail-type=ALL                                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<E-MAIL>                                   # Where to send mail 
#SBATCH --nodes=1                                              # Numbers of nodes
#SBATCH --ntasks=1 # Num of CPU                                # Number of CPU (tasks) 
#SBATCH --time=00:20:00 # Time limit                           # Time limit hrs:min:sec
#SBATCH --partition=short                                    # Partition/queue requested on server
#SBATCH --mem=100mb                                            # Memory limit 
#SBATCH --output=/scratch/Users/<username>/day8/eofiles/%x_%j.out   # Where to save out files
#SBATCH --error=/scratch/Users/<username>/day8/eofiles/%x_%j.err    # Where to save error files

echo "BEGIN: Preprocessing................."

################################################################                                                                                                                      
##Load required modules 
################################################################

module load preseq/2.0.3

################################################################
##Initialize folders 
################################################################
wd=/scratch/Users/<YourUsername>/day8
qc=${wd}/qc
bam=${wd}/bam
trim=${wd}/trim

#make folders that do not exist
mkdir -p ${qc} ${bam}

################################################################
##Job information
################################################################
printf "\nBAM Directory: ${bam}"
printf "\nQC Directory: ${qc}"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou have requested $SLURM_CPUS_ON_NODE core(s).\n"

################################################################
##Run script                                                     
################################################################

#5:check for complexity of files                         
##this program takes as input all your mapped bam files   
mkdir -p ${qc}/preseq

for sample in `ls ${trim} | tr '.' '\t' | cut -f 1`;
do
    ##i: generate the complexity curve
    # -s is the step size (for larger samples the default is 1million) 
    preseq c_curve -B -s 10000 \
	-o ${qc}/preseq/${sample}.c_curve.txt \
	${bam}/${sample}.sorted.bam

    ##ii: predict the complexity curve of a sequencing library using an initial experiment
    preseq lc_extrap -B \
	-o ${qc}/preseq/${sample}.lc_extrap.txt \
	${bam}/${sample}.sorted.bam

done

#6: combine all the QC metrics in one summary file with multiqc
multiqcdir=${wd}/multiqc
mkdir -p ${multiqcdir}

source /scratch/Shares/public/sread2024/generalpythonvenv/bin/activate
multiqc ${qc} -o ${multiqcdir}

echo "DONE: Preprocessing................."
echo Job finished at `date +"%T %a %d %b %Y"`
