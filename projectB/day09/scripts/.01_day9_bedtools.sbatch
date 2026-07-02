#!/bin/bash
#SBATCH --job-name=d9_bedtools                          # Job name
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=meco9877@colorado.edu        # Where to send mail
#SBATCH --nodes=1                                       # Numbers of nodes
#SBATCH --ntasks=1                                      # Number of CPU (tasks)
#SBATCH --time=00:10:00                                 # Time limit hrs:min:sec
#SBATCH --mem=2gb                                       # Memory limit
#SBATCH -p short
#SBATCH --output=/scratch/Users/meaghansullivan/day9/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/meaghansullivan/day9/eofiles/%x_%j.err


################## USER INPUT ######################################

chrom_sizes='/scratch/Shares/public/genomes/hisatfiles/hg38/hg38.chrom.sizes'
ANNOT='/scratch/Shares/public/sread/cookingShow/day9b/annotation/hg38_refseq_cleaned.sorted.bed'
CHIP_DIR='/scratch/Shares/public/sread/cookingShow/day8b/whole_genome/macs/'
HCT116_DMSO_PEAKS=${CHIP_DIR}/HCT116_DMSO_peaks_clean.narrowPeak
HCT116_NUTLIN_PEAKS=${CHIP_DIR}/HCT116_NUTLIN_peaks_clean.narrowPeak
results='/scratch/Users/meaghansullivan/day9/bedtools_results/'
TSS_ANNOT=${results}/hg38_refseq_tssregion.bed

mkdir -p ${results}
################## JOB INFO ########################################

printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

################### LOAD MODULES ###################################

module load bedtools/2.31.1
echo "Bedtools is loaded"
echo "Starting bedtools commands"

################## SCRIPT ##########################################

###-------------------------------------------------###
# Write your bedtools/file manipulation commands here #
###-------------------------------------------------###

# --- Which p53 peaks are associated with genes (not intergenic)

bedtools intersect -wa -a ${HCT116_NUTLIN_PEAKS} -b ${ANNOT} > ${results}/hct116_nutlin_p53_overlap_genes.bed

# --- Which p53 peeks are entirely within genes 

bedtools intersect -wa -f 1 -a ${HCT116_NUTLIN_PEAKS} -b ${ANNOT} > ${results}/hct116_nutlin_p53_in_genes.bed

# --- Which p53 peaks are intergenic 

bedtools intersect -wa -v -a ${HCT116_NUTLIN_PEAKS} -b ${ANNOT} > ${results}hct116_nutlin_p53_intergenic.bed

# --- Which genes overlap with p53 peaks? 

bedtools intersect -wb -a ${ANNOT} -b ${HCT116_NUTLIN_PEAKS} > ${results}/genes_overlap_p53peaks_hct116_nutlin.bed

# --- How many peaks or genes are in each of these files? 

cat ${results}/hct116_nutlin_p53_overlap_genes.bed | wc -l
cat ${results}/hct116_nutlin_p53_in_genes.bed | wc -l 
cat ${results}hct116_nutlin_p53_intergenic.bed | wc -l
cat ${results}/genes_overlap_p53peaks_hct116_nutlin.bed | wc -l

# --- Get TSS location 

cat "$ANNOT" \
| awk '{FS=OFS="\t"}; {if ($6 == "+") print $1,$2,($2+1),$4,".",$6}' \
> "$results"/hg38_refseq_tss_pos.bed
cat "$ANNOT" \
| awk '{FS=OFS="\t"}; {if ($6 == "-") print $1,($3-1),$3,$4,".",$6}' \
> "$results"/hg38_refseq_tss_neg.bed
cat "$results"/hg38_refseq_tss_pos.bed \
"$results"/hg38_refseq_tss_neg.bed \
> "$results"/hg38_refseq_tss.bed

# --- Create promoter annotation file 

bedtools slop -l 1000 -r 100 -s -i "$results"/hg38_refseq_tss.bed -g ${chrom_sizes} > ${results}/hg38_refseq_tssregion.bed

# --- Which genes have p53 ChIP sites in their promoters

bedtools intersect -wa ${HCT116_NUTLIN_PEAKS} -b ${TSS_ANNOT} > ${results}/hct116_nutlin_p53_overlap_promoters.bed

###------------------------------------------------###
###---- Now do for DMSO treated HCT116 Samples ----###
###------------------------------------------------###


# --- Which p53 peaks are associated with genes (not intergenic)

bedtools intersect -wa -a ${HCT116_DMSO_PEAKS} -b ${ANNOT} > ${results}/hct116_dmso_p53_overlap_genes.bed

# --- Which p53 peeks are entirely within genes

bedtools intersect -wa -f 1 -a ${HCT116_DMSO_PEAKS} -b ${ANNOT} > ${results}/hct116_dmso_p53_in_genes.bed

# --- Which p53 peaks are intergenic

bedtools intersect -wa -v -a ${HCT116_DMSO_PEAKS} -b ${ANNOT} > ${results}/hct116_dmso_p53_intergenic.bed

# --- Which genes overlap with p53 peaks?

bedtools intersect -wb -a ${ANNOT} -b ${HCT116_DMSO_PEAKS} > ${results}/genes_overlap_p53peaks_hct116_dmso.bed




###---- Jacquard Calculations ----###

# --- Jaccard used to calculate similarity between p53 peaks in DMSO and genes
bedtools jaccard -a ${HCT116_DMSO_PEAKS} -b ${ANNOT} > ${results}/jaccard_p53_peaks_in_genes_hct_dmso.txt

# --- Jaccard used to calculate similarity between p53 peaks in Nutlin and genes
bedtools jaccard -a ${HCT116_NUTLIN_PEAKS} -b ${ANNOT} > ${results}/jaccard_p53_peaks_in_genes_hct_nutlin.txt

# --- Shuffle p53 locations in the nutlin samples, then sort
bedtools shuffle -chrom -i ${HCT116_NUTLIN_PEAKS} -g ${chrom_sizes} > ${results}/hct116_nutlin_p53peaks_shuffled.bed
bedtools sort -k 1,1 -k2,2n -i ${results}/hct116_nutlin_p53peaks_shuffled.bed > ${results}/hct116_nutlin_p53peaks_shuffled.sorted.bed

# --- Similarity between shuffled p53 peaks in Nutlin and genes
bedtools jaccard -a ${results}/hct116_nutlin_p53peaks_shuffled.sorted.bed -b ${ANNOT} > ${results}/p53_nutlin_shuff_genes_jaccard.txt



