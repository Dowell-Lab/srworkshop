#!/bin/bash                                                                                                                                      
#SBATCH --job-name=submit_Rscript # Job name
#SBATCH --nodes=1                                                                                                                                #SBATCH --ntasks=1 # Number of CPU                                                                                                  
#SBATCH --time=00:02:00 # Time limit                                                                                                             #SBATCH --partition compute                                                                                                                      #SBATCH --mem=3mb # Memory                                                                                                                       #SBATCH --output=/path/to/e_and_o/%x_%j.out                                                                                                      #SBATCH --error=/path/to/e_and_o/%x_%j.err                                                                  
#SBATCH --mail-type=FAIL                      
#SBATCH --mail-user=you_name@email.address 


################################################################################                  
####################### initialize directories  ################################                                                                                                                                     

Rscript /path/to/Learning_R_submit_aws.R
