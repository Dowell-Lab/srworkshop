#!/bin/bash                                                                                                                                                      

#SBATCH --job-name=get_data # Job name
#SBATCH --nodes=1                                                                                                                                                
#SBATCH --ntasks=1 # Num of CPU                                                                                                  
#SBATCH --time=24:00:00 # Time limit                                                                                                                
#SBATCH --partition short                                                                                                                                        
#SBATCH --mem=10gb # Memory                                                                                                                                 
#SBATCH --output=/scratch/Shares/dowell/rutendo/projects/sread/sread2024/day8/e_and_o/%x_%j.out                                                                                                           
#SBATCH --error=/scratch/Shares/dowell/rutendo/projects/sread/sread2024/day8/e_and_o/%x_%j.err                                                                  
#SBATCH --mail-type=FAIL                      
#SBATCH --mail-user=rutendo.sigauke@cuanschutz.edu 


################################################################################
####################### initialize directories  ################################
wd=/scratch/Shares/dowell/rutendo/projects/sread/sread2024/day8
scripts=${wd}/scripts
sra=${wd}/sra

mkdir -p ${sra}

##############################
##summary of run
##############################
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

srrs=${scripts}/andrysik2017_chip.srr

echo "GET SRRs"

for srr in `cat ${srrs}`; do

    wget ${srr}
    mv ${srr} ${sra}
done

echo "DONE!"
