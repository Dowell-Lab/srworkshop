#!/bin/bash
#SBATCH --job-name=<job_name>
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<your_email>                        # Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:05:00                                 # Time limit hrs:min:sec
#SBATCH --partition=compute                             # Partition/queue requested on server    
#SBATCH --mem=3mb                                       # Memory limit
#SBATCH --output=/scratch/Users/<your_username>/day6/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<your_username>/day6/eofiles/%x_%j.err

################################################################################                  
####################### initialize directories  ################################                                                                                                                                     

Rscript /path/to/Learning_R_submit_aws.R
