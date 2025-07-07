#!/bin/bash

INDIR=/scratch/Shares/public/sread2025/data_files/day5/fastq/for_loops_fastq/ 
OUTDIR=/scratch/Users/<username>/workshop-day5/results

# makes a new directory if it does not already exist
mkdir -p ${OUTDIR}

# loops through each file that end in .end1.fastq in the in directory
for pathandfilename in ${INDIR}/*.end1.fastq
do 
# CODE TO OPERATE ON PATH_AND_FILENAME
# first get JUST the filename prefix (not including path or .end1.fastq)
FILENAME=$(basename $pathandfilename .end1.fastq)
echo $FILENAME

# now run the sbatch with the script d5-fastq-to-bam.sbatch.
# sbatch --export submits the script to the compute cluster with the variable values specified
# clarify the variable values for the indir, rootname, and outdir
# Which variable should change with each iteration the he loop? --> FILENAME?

sbatch --export=indir=$INDIR,rootname=$FILENAME,outdir=$OUTDIR d5-fastq-to-tdf.sbatch
done


echo "DONE!"