#!/bin/bash    
#SBATCH --job-name=check_data # Job name
#SBATCH --nodes=1                                                                                                                                                
#SBATCH --ntasks=1 # Num of CPU                                                                                                  
#SBATCH --time=24:00:00 # Time limit                                                                                                                
#SBATCH --partition short                                                                                                                                        
#SBATCH --mem=50mb # Memory                                                                                                                                 
#SBATCH --output=/scratch/Users/allenma/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/allenma/e_and_o/%x_%j.err                                                                  
#SBATCH --mail-type=FAIL                      
#SBATCH --mail-user=allenma@colorado.edu


indir=/Users/allenma/CNOT_and_brain/AGO2/fastqs/
outdir=/Users/allenma/CNOT_and_brain/AGO2/qc/
out=${outdir}afterdownload.md5sum.txt


for f in "${indir%/}"/*.fastq.gz; do
    [ -e "$f" ] || continue
    echo "Processing $f..."
    gunzip -c "$f" | md5sum | awk -v file="$(basename "$f")" '{print $1, file}' >> "$out"
done


