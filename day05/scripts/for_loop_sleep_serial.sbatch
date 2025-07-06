#!/bin/bash
#SBATCH --job-name=loopexample                      # Job name
#SBATCH --mail-type=ALL                             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=YOUREMAIL@colorado.edu          # Where to send mail
#SBATCH --nodes=1                                   # Numbers of nodes
#SBATCH --ntasks=1                                  # Number of CPU (tasks)
#SBATCH --time=00:10:00                             # Time limit hrs:min:sec
#SBATCH --partition=short                           # Partition/queue requested on server    
#SBATCH --mem=1gb                                   # Memory limit
#SBATCH --output=/path/to/your/e_o/%x_%j.out
#SBATCH --error=/path/to/your/e_o/%x_%j.err

################ DEFINE INPUT DIRECTORY #######################################

indir=/scratch/Shares/public/sread2025/data_files/day5/fastq/for_loops_fastq

################ RUN FOR LOOP #################################################

for filepath in "$indir"/*.fastq
do
echo "$filepath"
sleep 15
done
