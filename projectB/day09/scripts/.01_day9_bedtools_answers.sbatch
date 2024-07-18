#!/bin/bash
#SBATCH --job-name=d9_bedtools                          # Job name
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<EMAIL>        # Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:10:00                                 # Time limit hrs:min:sec
#SBATCH --mem=2gb                                       # Memory limit
#SBATCH -p short
#SBATCH --output=/scratch/Users/<username>/day9/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/<username>/day9/eofiles/%x_%j.err


################## USER INPUT ######################################

chrom_sizes='/scratch/Shares/public/genomes/hisatfiles/hg38/hg38.chrom.sizes'
ANNOT='/scratch/Shares/public/sread2024/cookingShow/day9b/annotation/hg38_refseq_cleaned.sorted.bed'
CHIP_DIR='/scratch/Shares/public/sread2024/cookingShow/day8b/whole_genome/macs/'
HCT116_DMSO_PEAKS=${CHIP_DIR}'/HCT116_DMSO_peaks_clean.narrowPeak'
HCT116_NUTLIN_PEAKS=${CHIP_DIR}'/HCT116_NUTLIN_peaks_clean.narrowPeak'
results='/scratch/Users/<username>/day9/bedtools_results/'
TSS_ANNOT=${results}/hg38_refseq_tssregion.bed

mkdir -p ${results}

################## JOB INFO ########################################

#printf "\nBAM Directory: $BAMS"
#printf "\nOutput Directory: $OUT"
#printf "\nOutput Rootname: $FILENAME"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

################### LOAD MODULES ###################################

module load bedtools/2.28.0

################## SET VARIABLES ###################################



################## SCRIPT ##########################################

# Write your bedtools/file manipulation commands here

# --- Get only Chr 21 from annotations and from peaks 
#grep '^chr21' ${ANNOT} > ${results}/chr21_annotations.bed

# --- Which p53 peaks are associated with genes (not intergenic)
bedtools intersect -wo -a ${HCT116_DMSO_PEAKS} -b ${ANNOT} > ${results}/p53_peaks_in_genes_hct_dmso.bed
cat ${results}/p53_peaks_in_genes_hct_dmso.bed | wc -l

# --- Which p53 peeks are entirely within genes 
bedtools intersect -f 1 -wo -a ${HCT116_DMSO_PEAKS} -b ${ANNOT} > ${results}/p53_peaks_completelyin_genes_hct_dmso.bed
cat ${results}/p53_peaks_completelyin_genes_hct_dmso.bed | wc -l

# --- Which p53 peaks are intergenic 
bedtools intersect -v -wo -a ${HCT116_DMSO_PEAKS} -b ${ANNOT} > ${results}/p53_peaks_intergenic_hct_dmso.bed
cat ${results}/p53_peaks_intergenic_hct_dmso.bed | wc -l

# --- Which genes overlap with p53 peaks? 
bedtools intersect -wo -a ${ANNOT} -b ${HCT116_DMSO_PEAKS} > ${results}/genes_with_p53_peaks_hct_dmso.bed
cat ${results}/genes_with_p53_peaks_hct_dmso.bed | wc -l


### --- At this point, have them download these to their computer --- ### 
### --- to look at in IGV.                                        --- ### 


# --- Which p53 peaks overlap promoters. 
# --- Promoters here will be defined as 1kb upstream of the TSS and 100bp downstream of TSS

cat ${ANNOT} \
| awk '{FS=OFS="\t"}; {if ($6 == "+") print $1,$2,($2+1),$4,".",$6}' \
> ${results}/hg38_refseq_tss_pos.bed

cat ${ANNOT} \
| awk '{FS=OFS="\t"}; {if ($6 == "-") print $1,($3-1),$3,$4,".",$6}' \
> ${results}/hg38_refseq_tss_neg.bed

cat ${results}/hg38_refseq_tss_pos.bed ${results}/hg38_refseq_tss_neg.bed \
> ${results}/hg38_refseq_tss.bed

bedtools slop -l 1000 -r 100 -s -i ${results}/hg38_refseq_tss.bed -g ${chrom_sizes} > ${results}/hg38_refseq_tssregion.bed

# --- Which genes have p53 ChIP sites in their promoters
bedtools intersect -wo -a ${HCT116_DMSO_PEAKS} -b ${TSS_ANNOT} > ${results}/p53_peaks_over_promoters.bed
cat ${results}/p53_peaks_over_promoters.bed | wc -l


### ------------------------------------------------------------------------------- ###
### ------------------------ Now for NUTLIN treated HCT116 ------------------------ ###
### ------------------------------------------------------------------------------- ###

# --- Which p53 peaks are associated with genes (not intergenic)
bedtools intersect -wo -a ${HCT116_NUTLIN_PEAKS} -b ${ANNOT} > ${results}/p53_peaks_in_genes_hct_nutlin.bed
cat ${results}/p53_peaks_in_genes_hct_nutlin.bed | wc -l

# --- Which p53 peeks are entirely within genes 
bedtools intersect -f 1 -wo -a ${HCT116_NUTLIN_PEAKS} -b ${ANNOT} > ${results}/p53_peaks_completelyin_genes_hct_nutlin.bed
cat ${results}/p53_peaks_completelyin_genes_hct_nutlin.bed | wc -l

# --- Which p53 peaks are intergenic 
bedtools intersect -v -wo -a ${HCT116_NUTLIN_PEAKS} -b ${ANNOT} > ${results}/p53_peaks_intergenic_hct_nutlin.bed
cat ${results}/p53_peaks_intergenic_hct_nutlin.bed | wc -l

# --- Which genes overlap with p53 peaks? 
bedtools intersect -wo -a ${ANNOT} -b ${HCT116_NUTLIN_PEAKS} > ${results}/genes_with_p53_peaks_hct_nutlin.bed
cat ${results}/genes_with_p53_peaks_hct_nutlin.bed | wc -l

# --- Which genes have p53 ChIP sites in their promoters
bedtools intersect -wo -a ${HCT116_NUTLIN_PEAKS} -b ${TSS_ANNOT} > ${results}/p53_peaks_over_promoters_hct_nutlin.bed
cat ${results}/p53_peaks_over_promoters_hct_nutlin.bed | wc -l



### ------------------------------------------------------------------------------- ###
### ------------------ Jaccard Calculation for HCT116 with Nutlin ----------------- ###
### ------------------------------------------------------------------------------- ###

# --- First need to sort for jaccard

# --- Jaccard used to calculate similarity between p53 peaks in DMSO and genes
bedtools jaccard -a ${HCT116_DMSO_PEAKS} -b ${ANNOT} > ${results}/jaccard_p53_peaks_in_genes_hct_dmso.txt

# --- Jaccard used to calculate similarity between p53 peaks in Nutlin and genes
bedtools jaccard -a ${HCT116_NUTLIN_PEAKS} -b ${ANNOT} > ${results}/jaccard_p53_peaks_in_genes_hct_nutlin.txt

# --- Shuffle p53 locations in the nutlin samples, then sort 
bedtools shuffle -chrom -i ${HCT116_NUTLIN_PEAKS} -g ${chrom_sizes} > ${results}/p53_nutlin_peaks_shuffled.bed
sort -k1,1 -k2,2n ${results}/p53_nutlin_peaks_shuffled.bed > ${results}/p53_nutlin_peaks_shuffled_sorted.bed

# --- Similarity between shuffled p53 peaks in Nutlin and genes
bedtools jaccard -a ${results}/p53_nutlin_peaks_shuffled_sorted.bed -b ${ANNOT} > ${results}/p53_nutlin_shuff_genes_jaccard.txt

