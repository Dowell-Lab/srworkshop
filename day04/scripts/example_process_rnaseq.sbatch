#!/bin/bash
#SBATCH --job-name=example_rnaseq_process       # Job name
#SBATCH --mail-type=ALL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOUR_EMAIL_HERE>           # Where to send mail
#SBATCH --nodes=1                               # Number of nodes requested
#SBATCH --ntasks=8                              # Number of CPUs (processor cores/tasks)
#SBATCH --mem=4gb                               # Memory limit
#SBATCH --time=01:30:00                        # Time limit hrs:min:sec
#SBATCH --partition=                       # Partition/queue requested on server
#SBATCH --output=       # Standard output:  the %x will be replaced with job_name and the %j by the job id
#SBATCH --error=        # Standard error log

### Displays the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

##############################
### Assigns path variables ###
##############################

## the fastq files will be used as input to fastqc and trimmomatic
## trimmed reads will then be passed on to the mapping step
FASTQ=/scratch/Shares/public/sread2024/homework_data_files/day4
DATADIR=                 #day4 scratch working directory
FASTQC=${DATADIR}/fastqc
HISAT2=${DATADIR}/hisat2
TRIM=${DATADIR}/trimmomatic

FILENAME=chr21Ethan_repA    # Prefix of the files to analyze
######################
### Loads modules   ##
######################
module load fastqc/0.11.5 
module load trimmomatic
module load hisat2/2.1.0
module load samtools/1.8

#######################
##Run the pipeline   ##
#######################
##1: Run fastqc on the samples (here run on example file ${FILENAME}.RNA.end1.fastq)

##2: Trim FASTQ Files

##3: Check Post-Trimming QC stats

##4: Map trimmed reads to reference genome

##5: Convert SAM file to BAM file and index the BAM files
##The BAM file is a compressed binary version of a SAM file, so takes up less space
##The sorted bam also takes up less space
###5.1 - convert SAM to BAM

###5.2 - sort bam file

###5.3 - index sorted bam file

##6: Backup and clean up
###6.1 - remove your SAM and unsorted BAM file

###6.2 - back up your scripts and results


## An alternate way to do step 6 is to backup everything that you need and then delete your
## entire scratch working directory. As long as the inputs are still present, you can
## regenerate the outputs from the rest of the pipeline

echo Job finished at `date +"%T %a %d %b %Y"`