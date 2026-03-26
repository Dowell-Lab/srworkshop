#!/bin/bash

################################################################################                                                                                                                        
####################### initialize directories  ################################
################################################################################                                                                                                                        
 
indir=/scratch/Users/ispe1418/disease_datasets/ptsd/scripts/
wholefileoflines=${indir}ptsd_download_links
outdir=/scratch/Users/ispe1418/disease_datasets/ptsd/fastq/
script_dir=/scratch/Users/ispe1418/disease_datasets/ptsd/scripts/

mkdir -p $outdir

echo "Getting the SRR Accession Numbers"


nlines=`cat $wholefileoflines | wc -l`

echo $nlines

for i in $(seq 1 $nlines); 
do echo $i; 
   sbatch --export=outdir=$outdir,infile=$wholefileoflines,whichline=$i,script_dir=$script_dir acommandsbatch.sbatch
   sleep 5  
done


echo "DONE!"
