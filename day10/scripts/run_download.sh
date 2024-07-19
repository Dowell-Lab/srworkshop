#!/bin/bash

################################################################################                                                                                                                        
####################### initialize directories  ################################
################################################################################                                                                                                                        
 
wd=/scratch/Users/your_user_name/day10
scripts=${wd}/scripts

srrs=${scripts}/andrysik2017_chip.srr

echo "Getting the SRR Accession Numbers"

##loop through the list of paths for SRRs                                                                                                                                                                 
for srr in `cat ${srrs}`; 
do

    # run the sbatch with the script
    # changing the paths based on the loop
    sbatch --export=srrpath=$srr 00_download_Andrysik2017ChIPseq.sbatch
    sleep 0.5s

done

echo "DONE!"
