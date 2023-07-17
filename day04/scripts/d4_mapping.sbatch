#!/bin/bash
#SBATCH --job-name=d4_mapping			       # Job name
#SBATCH --mail-type=ALL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOUR_EMAIL_HERE>          # Where to send mail
#SBATCH --nodes=1                              # Number of nodes requested
#SBATCH --ntasks=8                             # Number of CPUs (processor cores/tasks)
#SBATCH --mem=2gb                              # Memory limit
#SBATCH --time=01:30:00                        # Time limit hrs:min:sec
#SBATCH --partition=short                      # Partition/queue requested on server
#SBATCH --output=/scratch/Users/<USERNAME>/eofiles/%x.%j.out       # Standard output:  the %x will be replaced with job_name and the %j by the job id
#SBATCH --error=/scratch/Users/<USERNAME>/eofiles/%x.%j.err        # Standard error log

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
DATADIR=/scratch/Users/<USERNAME>/day4

HISAT2=${DATADIR}/hisat2
TRIM=${DATADIR}/trimmomatic

FILENAME=chr21Eric_repA

# Make directories

mkdir -p ${HISAT2}
mkdir -p ${TRIM}

######################
### Loads modules   ##
######################
module load hisat2/2.1.0
module load samtools/1.8

#######################
##Run the pipeline   ##
#######################
##4: Map trimmed reads to reference genome                                                                                                                                                     
hisat2 --very-fast -x /scratch/Shares/public/genomes/hisatfiles/hg38/HISAT2/genome \
-1 ${TRIM}/${FILENAME}.RNA.end1.trimmed.fastq \
-2 ${TRIM}/${FILENAME}.RNA.end2.trimmed.fastq \
> ${HISAT2}/${FILENAME}.RNA.sam \
2> ${HISAT2}/${FILENAME}.hisat2_maptstats.txt

##5: Convert SAM file to BAM file and index the BAM files
##The BAM file is a compressed binary version of a SAM file, so takes up less space
##The sorted bam also takes up less space
###5.1 - convert SAM to BAM 
samtools view -@ 8 -bS -o ${HISAT2}/${FILENAME}.RNA.bam ${HISAT2}/${FILENAME}.RNA.sam

###5.2 - sort bam file
samtools sort -@ 8 ${HISAT2}/${FILENAME}.RNA.bam > ${HISAT2}/${FILENAME}.RNA.sorted.bam

###5.3 - index sorted bam file
samtools index ${HISAT2}/${FILENAME}.RNA.sorted.bam  ${HISAT2}/${FILENAME}.RNA.sorted.bam.bai


echo Job finished at `date +"%T %a %d %b %Y"`
