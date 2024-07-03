#!/bin/bash
#SBATCH --job-name=<JOB ID>                                   # Job name
#SBATCH --mail-type=ALL                                       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<E-MAIL>                                  # Where to send mail
#SBATCH --nodes=1                                             # Numbers of nodes
#SBATCH --ntasks=1                                            # Number of CPU (tasks)
#SBATCH --time=00:30:00                                       # Time limit hrs:min:sec
#SBATCH --partition=compute                                   # Partition/queue requested on server    
#SBATCH --mem=2gb                                             # Memory limit
#SBATCH --output=/scratch/Users/<username>/day8/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<username>/day8/eofiles/%x_%j.err

echo "BEGIN: Preprocessing................."
################### LOAD MODULES ###################################

module load bedtools/2.25.0
module load meme/5.1.1

################## SET VARIABLES ###################################

#INDIR is where your macs output peak files are stored. OUTDIR is where I want my output from running meme and tomtom
WD=/scratch/Users/<YourUsername>/day8
INDIR=${WD}/macs
OUTDIR=${WD}/meme
mkdir -p ${OUTDIR}

HG38_FASTA=/scratch/Shares/public/genomes/chr21/chr21.fa

FILENAME=HCT116_NUTLIN

####################################################################
#### Job information
####################################################################
printf "\nMACS Peak Output Directory: ${INDIR}"
printf "\nOutput Directory: ${OUTDIR}"
printf "\nOutput File Name: $FILENAME"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

####################################################################
#### Get fasta of peak files
####################################################################
echo convert peaks call to fasta format
date 


bedtools getfasta \
 -fi ${HG38_FASTA} \
 -bed ${INDIR}/${FILENAME}_peaks_clean.narrowPeak \
 -fo ${OUTDIR}/${FILENAME}_peaks_clean.fasta

echo fasta file done 
date

####################################################################
#### Use MEME to find motifs in fasta sequences
####################################################################

#### Run MEME command line
## URL https://meme-suite.org/meme/tools/meme
#/opt/meme/5.1.1/bin/meme \
# ${OUTDIR}/${FILENAME}_peaks_clean.fasta \
# -dna \
# -oc ${OUTDIR}/${FILENAME} \
# -nostatus \
# -time 18000 \
# -mod zoops \
# -nmotifs 5 \
# -minw 6 \
# -maxw 50 \
# -objfun classic \
# -revcomp \
# -markov_order 0

####################################################################
#### Run TOMTOM to identify the motifs 
####################################################################

### Run TOMTOM command line
## URL https://meme-suite.org/meme/tools/tomtom

#mkdir -p ${WD}/tomtom

#/opt/meme/5.1.1/bin/tomtom \
# ${OUTDIR}/${FILENAME}/meme.html \
# -no-ssc \
# -oc ${WD}/tomtom/${FILENAME} \
# -verbosity 1 \
# -min-overlap 5 \
# -dist pearson \
# -evalue \
# -thresh 10.0 \
# /scratch/Shares/dowell/motifs/HOCOMOCODatabaseFIMO/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme

#echo MEME and TOMTOM done 
#date

echo "DONE: Preprocessing................."
echo Job finished at `date +"%T %a %d %b %Y"`
