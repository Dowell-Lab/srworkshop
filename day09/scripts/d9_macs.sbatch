#!/bin/bash
#SBATCH --job-name=<JOB_NAME>                           # Job name
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOUR_EMAIL>                        # Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:30:00                                 # Time limit hrs:min:sec
#SBATCH --partition=compute                             # Partition/queue requested on server    
#SBATCH --mem=2gb                                       # Memory limit
#SBATCH --output=/scratch/Users/<USERNAME>/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<USERNAME>/eofiles/%x_%j.err


################### LOAD MODULES ###################################

module load python/2.7.14
module load python/2.7.14/MACS/2.1.1
module load bedtools/2.25.0

################## SET VARIABLES ###################################

#INDIR is where my bams file are stored. OUTDIR is where I want my output from running MACS to go
INDIR='/path/to/bam/file'
OUTDIR='/path/to/your/macs/output/file'
mkdir -p ${OUTDIR}/w_ctrl ${OUTDIR}/no_ctrl

#Blacklist region bed file 
BLACKLIST='/scratch/Shares/dowell/sread/data_files/day9/blacklist_hg38/problematic_regions_hg38.bed'

#Prefix of your filename
FILENAME='BACH1'

################## JOB INFO ########################################

printf "\nBAM Directory: $BAMS"
printf "\nOutput Directory: $OUT"
printf "\nOutput Rootname: $FILENAME"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

#######################################################################

echo macs2
date
date

#### Call peaks with controls
# If you want to get broad peaks you can use the flag --broad
macs2 callpeak \
 -c ${INDIR}/${FILENAME}.input.chr21.sorted.bam \
 -t ${INDIR}/${FILENAME}.chr21.sorted.bam \
 --outdir ${OUTDIR}/w_ctrl \
 -n ${FILENAME} \
 -g hs \
 -B \
 -q 0.01 \

echo macs2 done
date
date

#### Removing ENCODE Blacklist regions
echo removing blacklist regions
date
date

bedtools intersect \
 -a ${OUTDIR}/w_ctrl/${FILENAME}_peaks.narrowPeak \
 -b ${BLACKLIST} \
 -v \
 > ${OUTDIR}/w_ctrl/${FILENAME}_peaks_clean.narrowPeak

echo blacklist regions removed
date
date


