#!/bin/bash
#SBATCH --job-name=<JOB_NAME>                        # Job name
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOUR_EMAIL>            	# Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=4                                      # Number of CPU (tasks)
#SBATCH --time=00:30:00                                 # Time limit hrs:min:sec
#SBATCH --partition=compute                            # Partition/queue requested on server    
#SBATCH --mem=10gb                                      # Memory limit
#SBATCH --output=/scratch/Users/<YOUR_USERNAME>/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<YOUR_USERNAME>/eofiles/%x_%j.err

################ SET REQUIRED VARIABLES #################################################


# First, we specify our file name, location, and outdir. ONLY change the outdir path to YOUR
# user directory. The rest of the script can be kept the same!

FILENAME=chr21Eric_repA.RNA
INDIR=/scratch/Shares/public/sread2023/cookingShow/day5/bam/
OUTDIR=/scratch/Users/<YOUR_USERNAME>/day5

# We'll also need to give the Genome size and info so that the file can be compressed into a TDF
genomeFasta=/scratch/Shares/public/genomes/hisatfiles/hg38/hg38.fa
genomeSizes=/scratch/Shares/public/genomes/hisatfiles/hg38/hg38.chrom.sizes

mkdir -p ${OUTDIR}

##########################################################################################
############### LOAD REQUIRED MODULES ####################################################
module load samtools/1.8
module load bedtools/2.25.0
module load igvtools/2.3.75

##########################################################################################
################## PRINT JOB INFO ########################################################

printf "\nfastq Directory: $INDIR"
printf "\nOutput Directory: $OUTDIR"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

##########################################################################################
################## RUN JOB ###############################################################

echo BedGraph and TDF

BEDGRAPHdir=${OUTDIR}/bedgraphForTdf
TDFdir=${OUTDIR}/tdf
mkdir -p ${BEDGRAPHdir} ${TDFdir}

##########################################################################################
################## PRINT JOB INFO ########################################################
# A bedgraph condenses a BAM file into a simple histogram of reads for each region in the genome.

# Because our file is paired end, we have to process each end separately. To make bedgraphs, we also need to
# process each strand separately. So, we have to make 4 different files:
# First, read1 bedgraph on the + strand. Read1 is reverse-stranded, so we have to flip it!
# Second, read1 bedgraph on the - strand. Read1 is reverse-stranded, so we have to flip it!
# Third, read2 bedgraph on the + strand
# Fourth, read2 bedgraph on the - strand
# Then, we'll concatenate them all together to make 1 big bedgraph file


# First, we extract Read1 from the bamfile. This is a little tricky to understand. SAMtools uses a code to determine read1 and read2. We give it that code
# using the -f parameter. 

# To learn more about these codes, visit http://davetang.org/wiki/tiki-index.php?page=SAMTools#Extracting_only_the_first_read_from_paired_end_BAM_files
# 0x0040 is hexadecimal for 64 (i.e. 16 * 4), which is binary for 1000000, corresponding to the first read in the pair.
samtools view \
        -h -b -f 0x0040 \
        ${INDIR}/${FILENAME}.sorted.bam \
        > ${BEDGRAPHdir}/${FILENAME}.pairfirst.bam

echo Finished pulling out the first read with flags.

# Now, we get Read2 from the bamfile. Same story, we need to use a SAMTools flag
# https://broadinstitute.github.io/picard/explain-flags.html
# 128 means second in pair, 128 in hexadecimal is 0x0080
samtools view \
        -h -b -f 0x0080 \
        ${INDIR}/${FILENAME}.sorted.bam \
        > ${BEDGRAPHdir}/${FILENAME}.pairsecond.bam

echo Finished pulling out the second read with flags.

# Now, we use genomeCoverageBed to start making bedGraphs. Remember, read1 is reverse-stranded!

genomeCoverageBed \
        -bg -split -strand - \
        -ibam ${BEDGRAPHdir}/${FILENAME}.pairfirst.bam \
        -g ${genomeSizes} \
        > ${BEDGRAPHdir}/${FILENAME}.pairfirst.pos.bed

echo Finished putting the read1 negative strand reads on the positive strand in the bedgraph file.

genomeCoverageBed \
        -bg  -split -strand + \
        -ibam ${BEDGRAPHdir}/${FILENAME}.pairfirst.bam \
        -g ${genomeSizes} | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }' \
        > ${BEDGRAPHdir}/${FILENAME}.pairfirst.neg.bed

echo Finished putting the read1 postive strand reads on the negative strand in the bedgraph file.

# Now, we repeat the process for Read2!

genomeCoverageBed \
        -bg -split -strand + \
        -ibam ${BEDGRAPHdir}/${FILENAME}.pairsecond.bam \
        -g ${genomeSizes} \
        > ${BEDGRAPHdir}/${FILENAME}.pairsecond.pos.bed

echo Finished putting the read2 postive strand reads on the postive strand in the bedgraph file.

genomeCoverageBed \
        -bg -split -strand - \
        -ibam ${BEDGRAPHdir}/${FILENAME}.pairsecond.bam \
        -g ${genomeSizes} | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }' \
        > ${BEDGRAPHdir}/${FILENAME}.pairsecond.neg.bed

echo Finished putting the read2 negative strand reads on the negative strand in the bedgraph file.

# Joining the read1 and read2 bedgraphs
    unionBedGraphs \
        -i ${BEDGRAPHdir}/${FILENAME}.pairfirst.pos.bed ${BEDGRAPHdir}/${FILENAME}.pairsecond.pos.bed \
        | awk -F '\t' {'print $1"\t"$2"\t"$3"\t"($4+$5)'} \
        > ${BEDGRAPHdir}/${FILENAME}.pos.bedGraph

    unionBedGraphs \
        -i ${BEDGRAPHdir}/${FILENAME}.pairfirst.neg.bed ${BEDGRAPHdir}/${FILENAME}.pairsecond.neg.bed \
        | awk -F '\t' {'print $1"\t"$2"\t"$3"\t"($4+$5)'} \
        > ${BEDGRAPHdir}/${FILENAME}.neg.bedGraph

echo Finished adding the postive and negative reads from both strands back together.

# Joining the + and - bedGraphs
cat \
        ${BEDGRAPHdir}/${FILENAME}.pos.bedGraph \
        ${BEDGRAPHdir}/${FILENAME}.neg.bedGraph \
        > ${BEDGRAPHdir}/${FILENAME}.bed

# Sorting the final Bedgraph
sortBed \
        -i ${BEDGRAPHdir}/${FILENAME}.bed \
        > ${BEDGRAPHdir}/${FILENAME}.BedGraph

echo Finished making the strand corrected sorted bedgraph.

# Compressing the bedGraph into a binary TDF file for easy viewing in IGV
igvtools toTDF \
        ${BEDGRAPHdir}/${FILENAME}.BedGraph \
        ${TDFdir}/${FILENAME}.tdf \
        ${genomeFasta}

echo Finished making TDF file.
echo Job finished at `date +"%T %a %d %b %Y"`
