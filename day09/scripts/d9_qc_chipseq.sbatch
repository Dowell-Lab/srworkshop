#!/bin/bash
#SBATCH --job-name=d9_qc                                # Job name
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<your_email>                 	# Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:15:00                                 # Time limit hrs:min:sec
#SBATCH --mem=50mb                                     # Memory limit
#SBATCH --output=/scratch/Users/<username>/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<username>/eofiles/%x_%j.err

################ SET REQUIRED VARIABLES ##################################################
## the fastq files will be used as input to fastqc. 
## output will be a fastqc file used to assess quality

INDIR=/scratch/Users/<username>/day9/bam
OUTDIR=/scratch/Users/<username>/day9/qc

##########################################################################################
############### LOAD REQUIRED MODULES ####################################################

module load preseq/2.0.3
module load samtools/1.8

##########################################################################################
################## PRINT JOB INFO ########################################################

printf "\nBAM Directory: $INDIR"
printf "\nOutput Directory: $OUTDIR"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

##########################################################################################
################## RUN JOB ###############################################################
## Preseq to check for complexity of files
## this program takes as input all your mapped bam files

mkdir -p ${OUTDIR}

################### Preseq
preseqdir=${OUTDIR}/preseq
mkdir -p ${preseqdir}

for FILENAME in `ls $INDIR | grep 'BACH1' | grep -v 'bai'`; do

    # generate the complexity curve
    # -s is the step size (for larger samples the default is 1million)
    preseq c_curve -B -s 10000 \
	-o ${preseqdir}/${FILENAME}.c_curve.txt \
	${INDIR}/${FILENAME}

    # predict the complexity curve of a sequencing library using an initial experiment 
    preseq lc_extrap -B \
	-o ${preseqdir}/${FILENAME}.lc_extrap.txt \
	${INDIR}/${FILENAME}
        
done
echo preseq complete

################## RSeQC
## Read distribution is one of the module in RSeQC that provides information of read distribution
## Not running because we the refseq file is not on the server and too large for github repo

#rseqcdir=${OUTDIR}/rseqc
#mkdir -p ${rseqcdir}
#module load python/2.7.14/rseqc/2.6.4

#for FILENAME in `ls ${BAMDIR} | grep -v 'bai'`;do

#    read_distribution.py --input-file=${BAMDIR}/${FILENAME} --refgene=<hg38_refseq.gtf> > ${rseqcdir}/${FILENAME}.read_dist.txt

#done
#echo rseqc report completed

################## MultiQC                                              
## Multiqc to summarize all the QC analysis (including fastQC and HISTAT2)
multiqcdir=${OUTDIR}/multiqc
mkdir -p ${multiqcdir}

## Activate virtual environment
source /data/multiqc-venv/bin/activate

multiqc ${OUTDIR} -o ${multiqcdir}

echo multiqc complete
echo Job finished at `date +"%T %a %d %b %Y"`
