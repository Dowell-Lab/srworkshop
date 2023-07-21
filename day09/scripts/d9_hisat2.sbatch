#!/bin/bash
#SBATCH --job-name=hisat2                               # Job name
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<your_email>                 	# Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=4                                      # Number of CPU (tasks)
#SBATCH --time=00:10:00                                 # Time limit hrs:min:sec
#SBATCH --mem=1gb                                     # Memory limit
#SBATCH --output=/scratch/Users/<username>/eofiles/d9_hisat2.out
#SBATCH --error=/scratch/Users/<username>/eofiles/d9_hisat2.err

################ SET REQUIRED VARIABLES ##################################################
## the fastq files will be used as input to fastqc. 
## output will be a fastqc file used to assess quality

INDIR=/scratch/Users/<username>/day9/fastq
BAM=/scratch/Users/<username>/day9/bam
SAM=/scratch/Users/<username>/day9/sam
QC=/scratch/Users/<username>/day9/qc/hisat_mapstats
INDICES=/scratch/Shares/public/genomes/hisatfiles/hg38/HISAT2/genome

##########################################################################################
############### LOAD REQUIRED MODULES ####################################################

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
#mkdir -p ${OUTDIR}

for FILENAME in `ls $INDIR | tr '.' '\t' | cut -f 1`; do

    hisat2 -p 4 \
	--very-sensitive \
	--no-spliced-alignment \
	-x ${INDICES} \
	-U ${INDIR}/${FILENAME}.fastq \
	--new-summary \
	> ${SAM}/${FILENAME}.sam \
	2> ${QC}/${FILENAME}.hisat2_mapstats.txt
        
    ###convert SAM to BAM                                    
    samtools view -@ 4 -h -bS -o ${BAM}/${FILENAME}.bam ${SAM}/${FILENAME}.sam

    ###sort BAM file            
    samtools sort -@ 4 ${BAM}/${FILENAME}.bam > ${BAM}/${FILENAME}.sorted.bam
    rm ${BAM}/${FILENAME}.bam
    
    ###index sorted bam file                                
    samtools index ${BAM}/${FILENAME}.sorted.bam  ${BAM}/${FILENAME}.sorted.bam.bai

done
                                                                        
echo Job finished at `date +"%T %a %d %b %Y"`
