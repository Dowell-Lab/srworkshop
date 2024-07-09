#!/bin/bash
#SBATCH --output=/scratch/Users/<username>/homework-day3/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<username>/homework-day3/eofiles/%x_%j.err
#SBATCH --mail-user=<email address>
#SBATCH --mail-type=END
#SBATCH --time=00:10:00
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=256Mb

date
sleep 120
date