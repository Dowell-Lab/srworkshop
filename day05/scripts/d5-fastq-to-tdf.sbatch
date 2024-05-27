#!/bin/bash
#SBATCH --job-name=fastqtotdf                        # Job name
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=YOUREMAIL@colorado.edu            	# Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=4                                      # Number of CPU (tasks)
#SBATCH --time=00:30:00                                 # Time limit hrs:min:sec
#SBATCH --partition=compute                            # Partition/queue requested on server    
#SBATCH --mem=10gb                                      # Memory limit
#SBATCH --output=/path/to/your/e_o/%x_%j.out
#SBATCH --error=/path/to/your/e_o/%x_%j.err

################ SET REQUIRED VARIABLES #################################################


FILENAME=$rootname
INDIR=$indir
BAM=$outdir/bams/
SAM=$outdir/sams/
QC=$outdir/qc/


INDICES=/scratch/Shares/public/genomes/hisatfiles/hg38/HISAT2/genome

mkdir -p $BAM $SAM $QC

##########################################################################################
############### LOAD REQUIRED MODULES ####################################################

module load hisat2/2.1.0
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

####################################################################################
##### Map the FASTQ files to get a SAM file and its compressed form (BAM file) #####
####################################################################################

# Mapping paired-end FASTQ files to make a SAM file.
hisat2 \
	--threads 4 \
	--very-sensitive \
	-x ${INDICES} \
	-1 ${INDIR}/${FILENAME}.end1.fastq \
	-2 ${INDIR}/${FILENAME}.end2.fastq \
	--new-summary \
	> ${SAM}/${FILENAME}.sam \
	2> ${QC}/${FILENAME}.hisat2_mapstats.txt
        
# Converting SAM to BAM.                               
samtools \
	view -@ 4 -h -bS \
	-o ${BAM}/${FILENAME}.bam \
	${SAM}/${FILENAME}.sam

# Sorting BAM file.         
samtools \
	sort -@ 4 \
	${BAM}/${FILENAME}.bam \
	> ${BAM}/${FILENAME}.sorted.bam

# Removing unsorted BAM file.
rm ${BAM}/${FILENAME}.bam
    
# Indexing sorted BAM file.                          
samtools \
	index \
	${BAM}/${FILENAME}.sorted.bam \
	${BAM}/${FILENAME}.sorted.bam.bai
                                                                        
echo BedGraph and TDF


####################################################################
##### Specify variables to make downstream code easier to read #####
####################################################################

BAMroot=$rootname
inDir=$BAM
outDir=$outdir
genomeFasta=/scratch/Shares/public/genomes/hisatfiles/hg38/hg38.fa
genomeSizes=/scratch/Shares/public/genomes/hisatfiles/hg38/hg38.chrom.sizes
BAM=${inDir}/${BAMroot}.sorted.bam

BEDGRAPHdir=${outDir}/bedgraphForTdf
TDFdir=${outDir}/tdf
STATSdir=${outDir}/stats
mkdir -p ${BEDGRAPHdir} ${TDFdir} ${STATSdir}




#################################################
##### Extract Reads 1 & 2 from the Bam file #####
#################################################

# First, we extract Read1 from the bamfile. 
# Because SAMtools uses a code to save information about reads, this can seem complicated. First, we need to get the filestats.
# Obtaining file stats from BAM input file.
samtools flagstat \
        -@ 8 ${BAM} \
        > ${STATSdir}/${BAMroot}.bam.flagstat \
        2>${STATSdir}/${BAMroot}.bam.flagstat.err

echo Finished obtaining BAM file stats.


# We then tell SAMtools that we want the first read of a pair via this code by giving it the code for a "first read in a pair": 40. We say this via the -f parameter. If you'd like to learn more
#   aboout these codes and how to use them, here's a cheat sheet: https://broadinstitute.github.io/picard/explain-flags.html



samtools view \
        -h -b -f 40 \
        ${BAM} \
        > ${BEDGRAPHdir}/${BAMroot}.pairfirst.bam

echo Finished pulling out the first read with flags 

# Now, we do the same thing but get the Read2 from the bamfile which has the code 128 (https://broadinstitute.github.io/picard/explain-flags.html).
samtools view \
        -h -b -f 128 \
        ${BAM} \
        > ${BEDGRAPHdir}/${BAMroot}.pairsecond.bam

echo Finished pulling out the second read with flags


####################################################################################
##### Make the bedGraphs with genomeCoverageBed for the First and Second Reads #####
####################################################################################

# This library is flipped which means that if the bam file says a read is mapped to the positive strand, 
#     that read is actually mapped to the negative strand. 
# Biologists: This is the case when scientists transform the RNA into cDNA. This allows a much more stable library, but means we have to flip everything for the first pair in the read.
# Non-biologists: Don't worry about the biologicial why, but this should remind you why it is imperative to UNDERSTAND the data you are analyzing, otherwise you can analyze it completely incorrectly but still get results.

# <ANNOTATE WHAT THE ${FILENAME} will provide >

############# First Reads In Pair (POS & NEG Strands) #############
genomeCoverageBed \
        -bg -split -strand - \
        -ibam ${BEDGRAPHdir}/${BAMroot}.pairfirst.bam \
        -g ${genomeSizes} \
        > ${BEDGRAPHdir}/${BAMroot}.pairfirst.pos.bed

echo Finished putting the read1 negative strand reads on the positive strand in the bedgraph file.

genomeCoverageBed \
        -bg  -split -strand + \
        -ibam ${BEDGRAPHdir}/${BAMroot}.pairfirst.bam \
        -g ${genomeSizes} | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }' \
        > ${BEDGRAPHdir}/${BAMroot}.pairfirst.neg.bed

echo Finished putting the read1 postive strand reads on the negative strand in the bedgraph file.

############# Second Reads In Pair (POS & NEG Strands) #############

genomeCoverageBed \
        -bg -split -strand + \
        -ibam ${BEDGRAPHdir}/${BAMroot}.pairsecond.bam \
        -g ${genomeSizes} \
        > ${BEDGRAPHdir}/${BAMroot}.pairsecond.pos.bed

echo Finished putting the read2 postive strand reads on the postive strand in the bedgraph file.

genomeCoverageBed \
        -bg -split -strand - \
        -ibam ${BEDGRAPHdir}/${BAMroot}.pairsecond.bam \
        -g ${genomeSizes} | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }' \
        > ${BEDGRAPHdir}/${BAMroot}.pairsecond.neg.bed


######################################################################
##### Combine the separate bedgraph files into ONE bedgraph file #####
######################################################################

##### Joining the read1 and read2 bedgraphs ######
# 1. unionBedGraphs combines the bedgraphs into chr"\t"start"\t"stop"\t"num_reads_overlapping_file1"\t"num_reads_overlapping_file2
# 2. The "|" takes this output and feeds into Awk. We tell Awk to print columns 1-3 with tabs in between. Then to add the last two columns together.
# 3. The ">" saves the file.

    unionBedGraphs \
        -i ${BEDGRAPHdir}/${FILENAME}.pairfirst.pos.bed ${BEDGRAPHdir}/${FILENAME}.pairsecond.pos.bed \
        | awk -F '\t' {'print $1"\t"$2"\t"$3"\t"($4+$5)'} \
        > ${BEDGRAPHdir}/${BAMroot}.pos.bedGraph

    unionBedGraphs \
        -i ${BEDGRAPHdir}/${FILENAME}.pairfirst.neg.bed ${BEDGRAPHdir}/${FILENAME}.pairsecond.neg.bed \
        | awk -F '\t' {'print $1"\t"$2"\t"$3"\t"($4+$5)'} \
        > ${BEDGRAPHdir}/${BAMroot}.neg.bedGraph

##### Joining the Positive and negative bedgraphs ######
cat \
        ${BEDGRAPHdir}/${BAMroot}.pos.bedGraph \
        ${BEDGRAPHdir}/${BAMroot}.neg.bedGraph \
        > ${BEDGRAPHdir}/${BAMroot}.bed

##### Sorting the final bedgraph ######
sortBed \
        -i ${BEDGRAPHdir}/${BAMroot}.bed \
        > ${BEDGRAPHdir}/${BAMroot}.BedGraph

echo Finished making the strand corrected sorted bedgraph.


##### Converting millionsmapped corrected Bedgraph into TDF. #####
igvtools toTDF \
        ${BEDGRAPHdir}/${BAMroot}.BedGraph \
        ${TDFdir}/${BAMroot}.tdf \
        ${genomeFasta}

echo Finished making TDF file.


echo Job finished at `date +"%T %a %d %b %Y"`
