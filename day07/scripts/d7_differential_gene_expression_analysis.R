# Differential Gene Expression Analysis

## Authors: Jacob Stanley (edited by Daniel Ramirez & Rutendo Sigauke)

##########################################
# Loading libraries.                     #
##########################################
library(DESeq2)


##########################################
# Loading the 2 input files.             #
##########################################

# 1. Load the metadata for the analysis

conditionsTableFile <- "/PATH/TO/day07/data/Andrysik2017_samples.tsv"
conditionsTable <- read.table(conditionsTableFile,
                              sep = "\t",
                              header = TRUE)
rownames(conditionsTable) <- conditionsTable$sample

# 2. Load the raw counts
geneCountsTableFile <- "/PATH/TO/day07/data/Andrysik2017_counts.tsv"
geneCountsTable <- read.table(geneCountsTableFile, 
                              header = TRUE, sep = "\t", fill = TRUE, 
                              stringsAsFactors = FALSE, na.strings = "")
#colnames(geneCountsTable) <- conditionsTable$sample

head(geneCountsTable)
conditionsTable

##########################################
# Loading files onto DESeq2              #
##########################################

dds <- DESeqDataSetFromMatrix(countData = geneCountsTable, 
                              colData = conditionsTable,
                              design = ~ condition)
dds
##########################################
# Removing lowly expressed genes.        #
##########################################

dds <- dds[rowSums(counts(dds)) > 1,]
dds

##########################################
# Running DESeq2.                        #  
##########################################

DEdds <- DESeq(dds)

##########################################
# Display size factors.                  #
##########################################

sizeFactors(DEdds)
colSums(counts(DEdds, normalized = FALSE))
colSums(counts(DEdds, normalized = TRUE))

################### ######################
# Plots dispersion estimates.            #
##########################################

plotDispEsts(DEdds, main = "Dispersion Estimates")


##########################################
# Set alpha and conditions.              #
##########################################

alphaValue <- 0.05
contrast <- c("condition", "nutlin", "dmso")

##########################################
# Extract and shrink results.            #
##########################################

results <- results(DEdds, alpha = alphaValue, contrast = contrast)
results_shrunk <- lfcShrink(DEdds, contrast = contrast, res = results)

##########################################
# Plots MA plot.                         #
##########################################

DESeq2::plotMA(results_shrunk,
               alpha = alphaValue,
               main = "RNA-seq\nHCT116\nDMSO vs Nutlin",
               xlab = "mean of normalized counts",
               ylab = "log fold change",
               ylim = c(-5,5))

##########################################
# Plots a given gene normalized counts.  #
##########################################

gene <- "CDKN1A"
plotCounts(DEdds, gene, 
           intgroup = "condition", 
           normalized = TRUE)


#########################################
# Sorts results by padj.                #
#########################################

results_shrunk <- results_shrunk[order(results_shrunk$padj), ]
                                 
                                 
#########################################
# Filters and stores results.           #
#########################################
results_shrunk_sig <- subset(results_shrunk, padj < alphaValue)

write.table(results_shrunk_sig, 
            sep = "\t",
            quote = FALSE, 
            row.names = TRUE,
            col.names = TRUE,
            "/PATH/TO/day07/data/Andrysik2017_RNAseq_Nutlin_results.tsv")

res_shrink_expressed <- as.data.frame(results_shrunk)
res_shrink_expressed <- res_shrink_expressed[!is.na(res_shrink_expressed$padj),]

outdir="/PATH/TO/day07/data/"

write.csv(rownames(res_shrink_expressed), file = paste0(outdir,"backgroundgenes.csv"),row.names = FALSE, col.names = FALSE, quote = FALSE)
write.csv(rownames(results_shrunk_sig), file = paste0(outdir,"siggenes.csv"),row.names = FALSE, col.names = FALSE, quote = FALSE)

rnkdf <- tibble(gene = rownames(results_shrunk),
				rnk = -log(res$pvalue) * sign(res$log2FoldChange)) %>%
	arrange(desc(rnk)) %>% drop_na()

## Write out the table without any additional information
write.table(rnkdf, file = paste0(outdir,"deseq_res_for_gsea.rnk"),
			append = FALSE, col.names = FALSE, row.names = FALSE,
			quote = FALSE, sep = "\t")
