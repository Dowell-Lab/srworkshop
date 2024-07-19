#!/bin/bash

################################################################################                                                                                                                        
####################### initialize directories  ################################
################################################################################                                                                                                                        
 
indir=/Users/<username>/workshopday10/scripts/
wholefileoflines=${indir}<shellscriptyoudownloaded>
outdir=/scratch/Users/<username>/day10/fastqs/

mkdir -p $outdir

echo "Getting the SRR Accession Numbers"


nlines=`cat $wholefileoflines | wc -l`

echo $nlines

for i in $(seq 1 $nlines); 
do echo $i; 
   sbatch --export=outdir=$outdir,infile=$wholefileoflines,whichline=$i acommandsbatch.sbatch
done


echo "DONE!"
