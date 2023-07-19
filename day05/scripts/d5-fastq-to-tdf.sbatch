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

##########################################################################################
################## RUN JOB ###############################################################

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

##########################################################################################
################## PRINT JOB INFO ########################################################

# Obtaining file stats from BAM input file.
samtools flagstat \
        -@ 8 ${BAM} \
        > ${STATSdir}/${BAMroot}.bam.flagstat \
        2>${STATSdir}/${BAMroot}.bam.flagstat.err

echo Finished obtaining BAM file stats.

# Obtaining flag from first read pair.
# http://davetang.org/wiki/tiki-index.php?page=SAMTools#Extracting_only_the_first_read_from_paired_end_BAM_files
# 0x0040 is hexadecimal for 64 (i.e. 16 * 4), which is binary for 1000000, corresponding to the read in the first read pair.
samtools view \
        -h -b -f 0x0040 \
        ${BAM} \
        > ${BEDGRAPHdir}/${BAMroot}.pairfirst.bam

echo Finished pulling out the first read with flags.

# Obtaining flag from second read pair.
# https://broadinstitute.github.io/picard/explain-flags.html
# 128 means second in pair, 128 in hexadecimal is 0x0080
samtools view \
        -h -b -f 0x0080 \
        ${BAM} \
        > ${BEDGRAPHdir}/${BAMroot}.pairsecond.bam

echo Finished pulling out the second read with flags.

# Using genomecoverage and swap the strand info on the second pair.
# Using genomecoverage on each of them.

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

echo Finished putting the read2 negative strand reads on the negative strand in the bedgraph file.

# Sorting the Bedgraphs
sortBed \
        -i ${BEDGRAPHdir}/${BAMroot}.pairfirst.pos.bed \
        > ${BEDGRAPHdir}/${BAMroot}.pairfirst.pos.BedGraph.sort

sortBed \
        -i ${BEDGRAPHdir}/${BAMroot}.pairfirst.neg.bed \
        > ${BEDGRAPHdir}/${BAMroot}.pairfirst.neg.BedGraph.sort

sortBed \
        -i ${BEDGRAPHdir}/${BAMroot}.pairsecond.pos.bed \
        > ${BEDGRAPHdir}/${BAMroot}.pairsecond.pos.BedGraph.sort

sortBed \
        -i ${BEDGRAPHdir}/${BAMroot}.pairsecond.neg.bed \
        > ${BEDGRAPHdir}/${BAMroot}.pairsecond.neg.BedGraph.sort

echo Finished sorting the bedgraphs.

# Adding the two Bedgraphs, putting the values in columns 4 and 5.
unionBedGraphs \
        -i ${BEDGRAPHdir}/${BAMroot}.pairfirst.pos.BedGraph.sort \
        ${BEDGRAPHdir}/${BAMroot}.pairsecond.pos.BedGraph.sort \
        > ${BEDGRAPHdir}/${BAMroot}.pos.Bedgraphcol

unionBedGraphs \
        -i ${BEDGRAPHdir}/${BAMroot}.pairfirst.neg.BedGraph.sort \
        ${BEDGRAPHdir}/${BAMroot}.pairsecond.neg.BedGraph.sort \
        > ${BEDGRAPHdir}/${BAMroot}.neg.Bedgraphcol

# Summing columnss 4 and 5.
awk -F '\t' '{OFS="\t"; print $1,$2,$3,$4+$5;}' \
        ${BEDGRAPHdir}/${BAMroot}.pos.Bedgraphcol \
        > ${BEDGRAPHdir}/${BAMroot}.pos.Bedgraph

awk -F '\t' '{OFS="\t"; print $1,$2,$3,$4+$5;}' \
        ${BEDGRAPHdir}/${BAMroot}.neg.Bedgraphcol \
        > ${BEDGRAPHdir}/${BAMroot}.neg.Bedgraph

echo Finished adding the postive and negative reads from both strands back together.

# Concatenating the two Bedgraphs files.
cat \
        ${BEDGRAPHdir}/${BAMroot}.pos.Bedgraph \
        ${BEDGRAPHdir}/${BAMroot}.neg.Bedgraph \
        > ${BEDGRAPHdir}/${BAMroot}.bed

# Sorting the final Bedgraph so it can be divided by millions mapped and converted into TDF.
sortBed \
        -i ${BEDGRAPHdir}/${BAMroot}.bed \
        > ${BEDGRAPHdir}/${BAMroot}.BedGraph

echo Finished making the strand corrected sorted bedgraph.

# Converting millionsmapped corrected Bedgraph into TDF.
igvtools toTDF \
        ${BEDGRAPHdir}/${BAMroot}.BedGraph \
        ${TDFdir}/${BAMroot}.tdf \
        ${genomeFasta}

echo Finished making TDF file.
echo Job finished at `date +"%T %a %d %b %Y"`
