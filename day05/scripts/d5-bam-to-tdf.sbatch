#!/bin/bash
#SBATCH --job-name=<JOB_NAME>                        # Job name
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<YOUR_EMAIL>            	# Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=4                                      # Number of CPU (tasks)
#SBATCH --time=00:30:00                                 # Time limit hrs:min:sec
#SBATCH --partition=short                            # Partition/queue requested on server    
#SBATCH --mem=10gb                                      # Memory limit
#SBATCH --output=/scratch/Users/<YOUR_USERNAME>/workshop-day5/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<YOUR_USERNAME>/workshop-day5/eofiles/%x_%j.err

################ SET REQUIRED VARIABLES #################################################


# First, we specify our file name, location, and outdir. ONLY change the outdir path to YOUR
# user directory. The rest of the script can be kept the same!

FILENAME=chr21Eric_repA.RNA
INDIR=/scratch/Shares/public/sread2025/cookingShow/day5/bam/
OUTDIR=/scratch/Users/<YOUR_USERNAME>/workshop-day5

# We'll also need to give the Genome size and info so that the file can be compressed into a TDF
genomeFasta=/scratch/Shares/public/genomes/hisatfiles/hg38/hg38.fa
genomeSizes=/scratch/Shares/public/genomes/hisatfiles/hg38/hg38.chrom.sizes

mkdir -p ${OUTDIR}

##########################################################################################
############### LOAD REQUIRED MODULES ####################################################
# These are the programs used to run our code
module load samtools/1.8
module load bedtools/2.25.0
module load igvtools/2.3.75

##########################################################################################
################## PRINT JOB INFO ########################################################
# This prints out the information about our job. It can be helpful to have this info later.
printf "\nfastq Directory: $INDIR"
printf "\nOutput Directory: $OUTDIR"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

##########################################################################################
################## GET DIRECTORIES READY #################################################
# We are setting directories to hold our Bedgraphs and TDFs once made

BEDGRAPHdir=${OUTDIR}/bedgraphForTdf
TDFdir=${OUTDIR}/tdf
mkdir -p ${BEDGRAPHdir} ${TDFdir}

####################################################################################
################## JOB INFO ########################################################
# A bedgraph condenses a BAM file into a simple histogram of reads for each region in the genome.
# Each line is in the format of chromosome"\t"start"\t"stop"\t"num_reads_overlapping_region, where \t means Tab
# Since our data is paired, that means that there are 2 reads per fragment of RNA captured (one starting from each end of the fragment). Therefore, we need to consider both reads in our analysis.

################################################
##### Extract Reads 1 & 2 from the Bamfile #####

# First, we extract Read1 from the bamfile. 
# Because SAMtools uses a code to save information about reads, this can seem complicated. We tell SAMtools that we want the first read of a pair
#    via this code by giving it the code for a "first read in a pair": 40. We say this via the -f parameter. If you'd like to learn more
#   aboout these codes and how to use them, here's a cheat sheet: https://broadinstitute.github.io/picard/explain-flags.html

# <ANNOTATE WHAT THE -h AND -b STAND FOR HERE> *Hint* module load samtools and run samtools view -help
# <ANNOTATE WHAT THE {INDIR}/${FILENAME}.sorted.bam and ${BEDGRAPHdir}/${FILENAME}.pairfirst.bam NAMES WILL BE>

samtools view \
        -h -b -f 64 \
        ${INDIR}/${FILENAME}.sorted.bam \
        > ${BEDGRAPHdir}/${FILENAME}.pairfirst.bam

echo Finished pulling out the first read with flags using BAM ${INDIR}/${FILENAME}.sorted.bam
echo Saved only first reads of Bam in ${BEDGRAPHdir}/${FILENAME}.pairfirst.bam

# Now, we do the same thing but get the Read2 from the bamfile which has the code 128 (https://broadinstitute.github.io/picard/explain-flags.html).
samtools view \
        -h -b -f 128 \
        ${INDIR}/${FILENAME}.sorted.bam \
        > ${BEDGRAPHdir}/${FILENAME}.pairsecond.bam

echo Finished pulling out the second read with flags using BAM ${INDIR}/${FILENAME}.sorted.bam
echo Saved only second reads of Bam in ${BEDGRAPHdir}/${FILENAME}.pairfirst.bam

####################################################################################
##### Make the bedGraphs with genomeCoverageBed for the First and Second Reads #####

# This library is "reverse stranded" which means that if the bam file says Read1 is mapped to the positive strand, 
#     that fragment actually matches the negative strand. 
# Biologists: This is the case when Read 2 is aligned to the RNA strand. More details on different library preps and how to figure them out can be found here: https://littlebitofdata.com/en/2017/08/strandness_in_rnaseq/ and here: https://www.reneshbedre.com/blog/determine-strandedness-rnaseq.html#google_vignette. 
# Non-biologists: You don't need to worry about the exact biologicial why although it's helpful, but this should remind you why it is imperative to UNDERSTAND the data you are analyzing, otherwise you can analyze it completely incorrectly but still get results.

# <ANNOTATE WHAT THE ${FILENAME} will provide >

############# First Reads In Pair (POS & NEG Strands) #############
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

############# Second Reads In Pair (POS & NEG Strands) #############

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

######################################################################
##### Combine the separate bedgraph files into ONE bedgraph file #####

##### Joining the read1 and read2 bedgraphs ######
# 1. unionBedGraphs combines the bedgraphs into chr"\t"start"\t"stop"\t"num_reads_overlapping_file1"\t"num_reads_overlapping_file2
# 2. The "|" takes this output and feeds into Awk. We tell Awk to print columns 1-3 with tabs in between. Then to add the last two columns together.
# 3. The ">" saves the file.

    unionBedGraphs \
        -i ${BEDGRAPHdir}/${FILENAME}.pairfirst.pos.bed ${BEDGRAPHdir}/${FILENAME}.pairsecond.pos.bed \
        | awk -F '\t' {'print $1"\t"$2"\t"$3"\t"($4+$5)'} \
        > ${BEDGRAPHdir}/${FILENAME}.pos.bedGraph

    unionBedGraphs \
        -i ${BEDGRAPHdir}/${FILENAME}.pairfirst.neg.bed ${BEDGRAPHdir}/${FILENAME}.pairsecond.neg.bed \
        | awk -F '\t' {'print $1"\t"$2"\t"$3"\t"($4+$5)'} \
        > ${BEDGRAPHdir}/${FILENAME}.neg.bedGraph

##### Joining the Positive and negative bedgraphs ######
cat \
        ${BEDGRAPHdir}/${FILENAME}.pos.bedGraph \
        ${BEDGRAPHdir}/${FILENAME}.neg.bedGraph \
        > ${BEDGRAPHdir}/${FILENAME}.bed

##### Sorting the final bedgraph ######
sortBed \
        -i ${BEDGRAPHdir}/${FILENAME}.bed \
        > ${BEDGRAPHdir}/${FILENAME}.BedGraph

echo Finished making the strand corrected sorted bedgraph.

######################################################################
##### Get the compressed version of the Bedgraph as TDF #####

# <ANNOTATE what ${genomeFasta} is referring to (hint look at the very top of the file>

# Compressing the bedGraph into a binary TDF file for easy viewing in IGV
igvtools toTDF \
        ${BEDGRAPHdir}/${FILENAME}.BedGraph \
        ${TDFdir}/${FILENAME}.tdf \
        ${genomeFasta}

echo Finished making TDF file.

# Just some helpful info to print out at the end of a job
echo Job finished at `date +"%T %a %d %b %Y"`
