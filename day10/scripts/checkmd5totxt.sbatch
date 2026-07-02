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


source ~/jhub_venv/bin/activate

#filereport_read_run=/Users/allenma/CNOT_and_brain/AGO2/metadata/filereport_read_run_SRP233483.tsv
#indir=/Users/allenma/CNOT_and_brain/AGO2/fastqs/
#outdir=/Users/allenma/CNOT_and_brain/AGO2/qc/
#out=${outdir}afterdownload.md5sum.txt

filereport_read_run=/Users/allenma/CNOT_and_brain/EDC3/metadata/filereport_read_run_PRJNA436576.tsv
indir=/Users/allenma/CNOT_and_brain/EDC3/fastqs/
outdir=/Users/allenma/CNOT_and_brain/EDC3/qc/
out=${outdir}afterdownload.md5sum.txt




#remove the check file if it exists
rm -f $out


#get the md5 sum
for f in "${indir%/}"/*.fastq.gz; do
    [ -e "$f" ] || continue
    echo "Processing $f..."
    md5sum "$f" | awk -v file="$(basename "$f")" '{print $1, file}' >> "$out"
done

#check it agains orginals
python3 checkmd5sumagainstori.py $out $filereport_read_run 
