#!/bin/bash
#SBATCH --job-name=d9_preprocess                        # Job name
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<email>                      	# Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:05:00                                 # Time limit hrs:min:sec
#SBATCH --mem=50mb                                     # Memory limit
#SBATCH --output=/scratch/Users/<username>/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<username>/eofiles/%x_%j.err

################ SET REQUIRED VARIABLES ##################################################
## the fastq files will be used as input to fastqc. 
## output will be a fastqc file used to assess quality

INDIR=/scratch/Users/<YourUsername>/day9/fastq
OUTDIR=/scratch/Users/<YourUsername>/day9/qc
BAM=/scratch/Users/<YourUsername>/day9/bam
SAM=/scratch/Users/<YourUsername>/day9/sam
INDICES=/scratch/Shares/public/genomes/hisatfiles/hg38/HISAT2/genome
OUTDIR=/scratch/Users/<YourUsername>/day9/qc

##########################################################################################
############### LOAD REQUIRED MODULES ####################################################

module load fastqc/0.11.5
module load hisat2/2.1.0
module load samtools/1.8

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

################## FASTQC
## Use fastqc to check read quality
## The fastq files will be used as input to fastqc. 
## Output will be a fastqc file used to assess quality

mkdir -p ${OUTDIR}/fastqc

for FILENAME in `ls $INDIR`; do
    
    fastqc \
	${INDIR}/${FILENAME} \
	-o ${OUTDIR}/fastqc
        
done
echo fastqc complete

################## HISAT2
## Align reads to reference genome using HISATI2
## We are mapping reads to the genome, but since this is ChIP-seq,
## we will ignore splicing aware features for the mapper

HISATQC=${OUTDIR}/hisat_mapstats
mkdir -p ${HISATQC}

for FILENAME in `ls $INDIR | grep 'BACH1' | tr '.' '\t' | cut -f 1`; do

    hisat2 -p 4 \ 
        --very-sensitive \ 
        --no-spliced-alignment \ 
        -x ${INDICES} \ 
        -U ${INDIR}/${FILENAME}.fastq \ 
        --new-summary \ 
        > ${SAM}/${FILENAME}.sam \ 
        2> ${HISATQC}/${FILENAME}.hisat2_mapstats.txt

    ###convert SAM to BAM
    samtools view -@ 4 -h -bS -o ${BAM}/${FILENAME}.bam ${SAM}/${FILENAME}.sam

    ###sort BAM file
    samtools sort -@ 4 ${BAM}/${FILENAME}.bam > ${BAM}/${FILENAME}.sorted.bam
    rm ${BAM}/${FILENAME}.bam

    ###index sorted bam file
    samtools index ${BAM}/${FILENAME}.sorted.bam  ${BAM}/${FILENAME}.sorted.bam.bai

done

echo hisat2 alignment and bam to sam complete

echo Job finished at `date +"%T %a %d %b %Y"`
