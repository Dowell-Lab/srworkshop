#!/bin/bash 
#SBATCH --job-name=featurecounts # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<your_email> # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=10gb # Memory limit
#SBATCH --time=10:00:00 # Time limit hrs:min:sec
#SBATCH --output=/path/to/eando/%x.%j.out # Standard output
#SBATCH --error=/path/to/eando/%x.%j.err # Standard error log


echo $rootname
echo $indir
echo $outdir



Rscript --vanilla RNAseq_featurecounts.R $rootname $indir $outdir

