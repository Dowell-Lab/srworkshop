#!/bin/bash                                                                                                                                                      

#SBATCH --job-name=get_data # Job name
#SBATCH --nodes=1                                                                                                                                                
#SBATCH --ntasks=1 # Num of CPU                                                                                                  
#SBATCH --time=24:00:00 # Time limit                                                                                                                
#SBATCH --partition short                                                                                                                                        
#SBATCH --mem=50mb # Memory                                                                                                                                 
#SBATCH --output=/path/to/e_and_o/%x_%j.out                                                                                                           
#SBATCH --error=/path/to/e_and_o/%x_%j.err                                                                  
#SBATCH --mail-type=FAIL                      
#SBATCH --mail-user=username@email.edu 


################################################################################
####################### initialize directories  ################################
SRRPATH=$srrpath

##############################
##summary of run
##############################
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

# use wget to download the files
${wholeline}

echo "DONE!"
