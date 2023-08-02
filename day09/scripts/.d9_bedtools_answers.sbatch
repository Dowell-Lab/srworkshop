#!/bin/bash
#SBATCH --job-name=d9_bedtools                              # Job name
#SBATCH --mail-type=NONE                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lynn.sanford@colorado.edu                        # Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:10:00                                 # Time limit hrs:min:sec
#SBATCH --mem=2gb                                       # Memory limit
#SBATCH --output=/scratch/Users/lynn-sanford/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/lynn-sanford/eofiles/%x_%j.err

################## JOB INFO ########################################

printf "\nBAM Directory: $BAMS"
printf "\nOutput Directory: $OUT"
printf "\nOutput Rootname: $FILENAME"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

################### LOAD MODULES ###################################

module load bedtools/2.25.0

################## SET VARIABLES ###################################

chrom_sizes='/scratch/Shares/public/genomes/hisatfiles/hg38/hg38.chrom.sizes'
chr21_genes='/scratch/Users/lynn-sanford/day9/bedtools/hg38_refseq_chr21.bed'
BACH1_peaks='/scratch/Users/lynn-sanford/day9/bedtools/BACH1_peaks_chr21.bed'
results='/scratch/Users/lynn-sanford/day9/bedtools'

mkdir -p "$results"

################## SCRIPT ##########################################

# Write your bedtools/file manipulation commands here

# Bedtools intersect
bedtools intersect -wa -a "$BACH1_peaks" -b "$chr21_genes" \
  > "$results"/BACH1_peaks_over_genes.bed

bedtools intersect -f 1 -wa -a "$BACH1_peaks" -b "$chr21_genes" \
  > "$results"/BACH1_peaks_completely_over_genes.bed

bedtools intersect -v -wa -a "$BACH1_peaks" -b "$chr21_genes" \
  > "$results"/BACH1_peaks_intergenic.bed

bedtools intersect -wa -a "$chr21_genes" -b "$BACH1_peaks" \
  > "$results"/genes_over_BACH1_peaks.bed

# Awk manipulation
cat "$chr21_genes" \
  | awk '{FS=OFS="\t"}; {if ($6 == "+") print $1,$2,($2+1),$4,".",$6}' \
  > "$results"/hg38_refseq_chr21_tss_pos.bed

cat "$chr21_genes" \
  | awk '{FS=OFS="\t"}; {if ($6 == "-") print $1,($3-1),$3,$4,".",$6}' \
  > "$results"/hg38_refseq_chr21_tss_neg.bed

cat "$results"/hg38_refseq_chr21_tss_pos.bed \
  "$results"/hg38_refseq_chr21_tss_neg.bed \
  > "$results"/hg38_refseq_chr21_tss.bed

# Bedtools slop
bedtools slop -l 500 -r 100 -s -i "$results"/hg38_refseq_chr21_tss.bed -g "$chrom_sizes" \
  > "$results"/hg38_refseq_chr21_promoter.bed

# Bedtools intersect
bedtools intersect -wa -a "$BACH1_peaks" -b "$results"/hg38_refseq_chr21_promoter.bed \
  > "$results"/BACH1_peaks_over_promoters.bed

# Bedtools jaccard and shuffle
bedtools jaccard -a "$BACH1_peaks" -b "$chr21_genes" > "$results"/BACH1_genes_jaccard.txt

bedtools shuffle -chrom -i "$BACH1_peaks" -g "$chrom_sizes" > "$results"/BACH1_peaks_shuffled.bed

sort -k1,1 -k2,2n "$results"/BACH1_peaks_shuffled.bed > "$results"/BACH1_peaks_shuffled_sorted.bed

bedtools jaccard -a "$results"/BACH1_peaks_shuffled_sorted.bed -b "$chr21_genes" > "$results"/BACH1shuff_genes_jaccard.txt







