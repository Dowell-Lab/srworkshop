#!/bin/bash

################################################################################                                                                                                                        
####################### initialize directories  ################################
################################################################################                                                                                                                        
 
indir=/Users/<username>/workshopday10/scripts/
wholefileoflines=ena-file-download-read_run-PRJNA785930-fastq_ftp-20240719-1321.sh
outdir=/Users/<username>/day10/fastqs/

mkdir -p outdir

echo "Getting the SRR Accession Numbers"

##loop through the list of paths for SRRs                                                                                                                                                                 
for wholeline in `cat ${wholefileoflines}`; 
do

    # run the sbatch with the script
    # changing the paths based on the loop
    sbatch --export=wholeline=$wholeline acommandsbatch.sbatch
    sleep 0.5s

done

echo "DONE!"
