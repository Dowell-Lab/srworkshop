#!/bin/bash                                                                                                                                                      

#SBATCH --job-name=get_data # Job name
#SBATCH --nodes=1                                                                                                                                                
#SBATCH --ntasks=1 # Num of CPU                                                                                                  
#SBATCH --time=24:00:00 # Time limit                                                                                                                
#SBATCH --partition short                                                                                                                                        
#SBATCH --mem=50mb # Memory                                                                                                                                 
#SBATCH --output=/scratch/Users/ispe1418/disease_datasets/ptsd/eando/%x_%j.out
#SBATCH --error=/scratch/Users/ispe1418/disease_datasets/ptsd/eando/%x_%j.err                                                                  
#SBATCH --mail-type=FAIL                      
#SBATCH --mail-user=ispe1418@colorado.edu 


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

echo $outdir
echo $infile
echo $whichline

cd $outdir

line_content=$(awk -v line="$whichline" 'NR == line {print $0}' "$infile")

# Extract the URL (last field of the line) and then its basename
url=${line_content##* }                    # everything after last space
fname=${url##*/}                           # SRR5412628.fastq.gz

# Example: start a new script for each wget
cat "${script_dir}header.txt"                >  "${outdir}shellscript${whichline}.sh"

# Define full path to the expected file
printf 'fullpath="%s"\n' "${outdir}${fname}" >> "${outdir}shellscript${whichline}.sh"

# Only run the download command if the file does NOT already exist
printf 'if [ ! -f "$fullpath" ]; then\n' >> "${outdir}shellscript${whichline}.sh"
printf '  %s\n' "$line_content" >> "${outdir}shellscript${whichline}.sh"
printf 'else\n' >> "${outdir}shellscript${whichline}.sh"
printf '  echo "Downloaded file already exists: $fullpath"\n' >> "${outdir}shellscript${whichline}.sh"
printf 'fi\n' >> "${outdir}shellscript${whichline}.sh"


sbatch ${outdir}shellscript${whichline}.sh

echo "DONE!"
