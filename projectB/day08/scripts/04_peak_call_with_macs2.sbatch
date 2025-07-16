#!/bin/bash
#SBATCH --job-name=<JOB ID>                                  # Job name
#SBATCH --mail-type=ALL                                      # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<E-MAIL>                                 # Where to send mail
#SBATCH --nodes=1                                            # Numbers of nodes
#SBATCH --ntasks=1                                           # Number of CPU (tasks)
#SBATCH --time=00:30:00                                      # Time limit hrs:min:sec
#SBATCH --partition=short                                  # Partition/queue requested on server    
#SBATCH --mem=2gb                                            # Memory limit
#SBATCH --output=/scratch/Users/<username>/day8/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<username>/day8/eofiles/%x_%j.err

echo "BEGIN: Preprocessing................."

################### LOAD MODULES ###################################
module load bedtools/2.25.0
source /scratch/Shares/public/sread2025/generalpythonvenv/bin/activate

################## SET VARIABLES ###################################
#INDIR is where my bams file are stored. OUTDIR is where I want my output from running MACS to go
WD=/scratch/Users/<YourUsername>/day8
INDIR=${WD}/bam
OUTDIR=${WD}/macs
mkdir -p ${OUTDIR} 

#Blacklist region bed file 
BLACKLIST=/scratch/Shares/public/sread2025/data_files/day8b/chip_blacklist/blacklist_hg38.bed

#Samples (HCT116 samples are below)
INPUT=SRR4090089_chr21
DMSO=SRR4090090_chr21
NUTLIN=SRR4090091_chr21

#Prefix of your filenames
CELL=HCT116
FILENAME_DMSO=${CELL}_DMSO
FILENAME_NUTLIN=${CELL}_NUTLIN

################## JOB INFO ########################################

printf "\nBAM Directory: $INDIR"
printf "\nOutput Directory: $OUTDIR"
printf "\nOutput Rootname: $CELL"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

###################################################################
##INPUT vs DMSO
###################################################################
echo macs3 HCT116 DMSO
date

#### Call peaks with controls
# If you want to get broad peaks you can use the flag --broad
macs3 callpeak \
 -c ${INDIR}/${INPUT}.sorted.bam \
 -t ${INDIR}/${DMSO}.sorted.bam \
 --outdir ${OUTDIR} \
 -n ${FILENAME_DMSO} \
 -g hs \
 -B \
 -q 0.01 \

echo macs3 done
date

#### Removing ENCODE Blacklist regions
echo removing blacklist regions
date

bedtools intersect \
 -a ${OUTDIR}/${FILENAME_DMSO}_peaks.narrowPeak \
 -b ${BLACKLIST} \
 -v \
 > ${OUTDIR}/${FILENAME_DMSO}_peaks_clean.narrowPeak

echo blacklist regions removed
date

####################################################################
##INPUT vs NUTLIN
####################################################################
echo macs3 HCT116 NUTLIN
date

#### Call peaks with controls
# If you want to get broad peaks you can use the flag --broad
macs3 callpeak \
 -c ${INDIR}/${INPUT}.sorted.bam \
 -t ${INDIR}/${NUTLIN}.sorted.bam \
 --outdir ${OUTDIR} \
 -n ${FILENAME_NUTLIN} \
 -g hs \
 -B \
 -q 0.01 \

echo macs3 done
date

#### Removing ENCODE Blacklist regions
echo removing blacklist regions
date

bedtools intersect \
 -a ${OUTDIR}/${FILENAME_NUTLIN}_peaks.narrowPeak \
 -b ${BLACKLIST} \
 -v \
 > ${OUTDIR}/${FILENAME_NUTLIN}_peaks_clean.narrowPeak

echo blacklist regions removed

echo "DONE: Preprocessing................."
echo Job finished at `date +"%T %a %d %b %Y"`


