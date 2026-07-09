#!/bin/bash
#SBATCH --job-name=cellrangerCount_male19                           # Job name
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=18:00:00                                 # Time limit hrs:min:sec
#SBATCH --partition=short                               # Partition/queue requested on server    
#SBATCH --mem=10gb                                       # Memory limit
#SBATCH --output=/Users/<github_username>/workshop-day7/cellranger_count/eando/cellrangercount.%j.out
#SBATCH --error=/Users/<github_username>/workshop-day7/cellranger_count/eando/cellrangercount.%j.err

##########################################################################################
################## MODULE LOAD ###################################################

module load cellranger/7.2.0

##########################################################################################
################## DISPLAY JOB CONTEXT ###################################################

echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

##########################################################################################
################## RUN JOB ###############################################################


#makes an eando directory unless one already exists
mkdir -p ~/workshop-day7/cellranger_count/eando

cellranger count --id=T21BM_male19 \
 --fastqs=/scratch/Shares/public/sread/data_files/day7a/fastq/sampled_fastq \
 --transcriptome=/scratch/Shares/public/sread/cookingShow/day7a/genomes/refdata-gex-GRCh38-2020-A \
 --sample=DSOX19_1 \
 --include-introns=true \
 --localcores=1    ###must match ntasks above

                                                                                
echo Job finished at `date +"%T %a %d %b %Y"`
