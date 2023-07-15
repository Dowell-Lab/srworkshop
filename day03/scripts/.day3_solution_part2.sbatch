#!/bin/bash
#SBATCH --output=/scratch/Users/zmaas/workshop-day3/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zmaas/workshop-day3/e_and_o/%x_%j.err
#SBATCH --mail-user=zachary.maas@colorado.edu
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=256mb

module load fastqc
cd /scratch/Users/zmaas/workshop-day3/data
fastqc SRR062641.filt.fastq.gz
