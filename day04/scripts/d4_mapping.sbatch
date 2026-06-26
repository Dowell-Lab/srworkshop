#!/bin/bash
#SBATCH --job-name=<JOB_NAME>			       # Job name
#SBATCH --mail-type=END                       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOUR_EMAIL_HERE>          # Where to send mail
#SBATCH --nodes=1                              # Number of nodes requested
#SBATCH --ntasks=<CORES>                             # Number of CPUs (processor cores/tasks)
#SBATCH --mem=<MEMORY>                              # Memory limit
#SBATCH --time=<WALLTIME>                        # Time limit hrs:min:sec
#SBATCH --partition=short                      # Partition/queue requested on server
#SBATCH --output=/scratch/Users/<USERNAME>/workshop-day4/eofiles/%x.%j.out       # Standard output:  the %x will be replaced with job_name and the %j by the job id
#SBATCH --error=/scratch/Users/<USERNAME>/workshop-day4/eofiles/%x.%j.err        # Standard error log

################## Display job context ################################
## This script maps the trimmed fastq files and outputs several file formats
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

################## Assign variables ###################################
datadir=/scratch/Users/<USERNAME>/workshop-day4
hisat2="$datadir"/hisat2
trim="$datadir"/trimmomatic

filebase=chr21Eric_repA.RNA  # Prefix of the files to analyze
hisat2_indexbase=/scratch/Shares/public/genomes/hisatfiles/hg38/HISAT2/genome # Genome index path/prefix

# Make directories
mkdir -p $hisat2

################## Load modules/software ##############################
# Load module for hisat2
module load hisat2/2.2.1

# On some clusters there might be a samtools module
# But samtools v1.23.1 is just installed system-wide here

################## Run the pipeline ###################################
##4: Map trimmed reads to reference genome                                                                                                                                                     
hisat2 -p $SLURM_NTASKS --very-fast -x $hisat2_indexbase \
-1 "$trim"/"$filebase".end1.trimmed.fastq \
-2 "$trim"/"$filebase".end2.trimmed.fastq \
> "$hisat2"/"$filebase".sam \
2> "$hisat2"/"$filebase".hisat2_mapstats.txt

##5: Convert SAM file to BAM file and index the BAM files
##The BAM file is a compressed binary version of a SAM file, so takes up less space
##The sorted bam also takes up less space
###5.1 - convert SAM to BAM 
samtools view -@ $SLURM_NTASKS -bS -o "$hisat2"/"$filebase".bam "$hisat2"/"$filebase".sam

###5.2 - sort bam file
samtools sort -@ $SLURM_NTASKS "$hisat2"/"$filebase".bam > "$hisat2"/"$filebase".sorted.bam

###5.3 - index sorted bam file
samtools index "$hisat2"/"$filebase".sorted.bam "$hisat2"/"$filebase".sorted.bam.bai

echo Job finished at `date +"%T %a %d %b %Y"`
