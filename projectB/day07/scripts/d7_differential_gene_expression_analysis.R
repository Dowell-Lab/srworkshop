# Differential Gene Expression Analysis

## Authors: Jacob Stanley (edited by Daniel Ramirez, Rutendo Sigauke, & Samuel Hunter)

##########################################
# Loading libraries.                     #
##########################################
library(DESeq2)
library(pheatmap)
library(ggplot2)
##########################################
# Loading the 2 input files.             #
##########################################

# 1. Load the metadata for the analysis

conditionsTableFile <- "/PATH/TO/day7/Andrysik2017_samples.tsv"
conditionsTable <- read.table(conditionsTableFile,
                              sep = "\t",
                              header= TRUE)

conditionsTable

# 2. Load the raw counts

geneCountsTableFile <- "/PATH/TO/day7/Andrysik2017_counts.csv"
geneCountsTable <- read.table(geneCountsTableFile,
                              header=TRUE,
                              row.names = "GeneID",
                              sep = "\t",
                              stringsAsFactors = FALSE)
head(geneCountsTable)

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
# Plots for EDA                          #
##########################################

normcounts <- log2(as.data.frame(counts(DEdds, normalized = TRUE)) + 1)
hist(normcounts)
boxplot(normcounts)

rld <- rlog(DEdds)
DESeq2::plotPCA(rld, intgroup=c("replicate"))

select <- order(rowMeans(counts(DEdds,normalized=TRUE)),
                decreasing=TRUE)[1:60]

columns_for_heatmap <- as.data.frame(colData(DEdds)[,c("condition","replicate")])
pheatmap(normcounts[select,], 
         cluster_rows=TRUE, 
         show_rownames=TRUE,
         cluster_cols=FALSE, 
         annotation_col=columns_for_heatmap)

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
results_shrunk <- lfcShrink(DEdds, contrast = contrast, res = results, type = "normal")

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

##########################################
# Plots Volcano plot.                    #
##########################################
results_shrunk$threshold <- ifelse(results_shrunk$padj < alphaValue, "significant", "not_significant")
ggplot(results_shrunk) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  ggtitle("Volcano Plot") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme_classic()


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
