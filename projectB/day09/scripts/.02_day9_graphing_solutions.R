


# Set environment ---------------------------------------------------------

  library(tidyverse) # load ggplot2 
  setwd('your/path/to/<your short read dir>/day9/') # set your local working directory 
  outdir <- 'results/'
  

# Read in data ------------------------------------------------------------

  ### --- READ IN CHIP DATA 
  dmso_peaks <- read_tsv('data/bedtools_results/hct116_dmso_p53_overlap_genes.bed', col_names = FALSE) 
  nutlin_peaks <- read_tsv('data/bedtools_results/hct116_nutlin_p53_overlap_genes.bed', col_names = FALSE)
  
  ### --- GIVE CHIP DATA USEFUL COLUMN NAMES
  set_colnames <- c('chip_chr', 'chip_start', 'chip_end', 'peak_id', 'score', 
                    'chip_strand', 'signalValue', 'log_pval', 'log_qval', 
                    'summit', 'gene_chr', 'gene_start', 'gene_end', 
                    'gene', '.', 'gene_strand', 'overlap') 
  colnames(dmso_peaks) <- set_colnames
  colnames(nutlin_peaks) <- set_colnames
  
  ### --- READ IN DE DATA 
  de_import <- read_tsv('data/deseq_output/hct116_deres.txt')# read in data here 
  
  ### --- EXTRACT SIGNIFICANTLY CHANGING GENES
  de <- de_import[de_import$padj < 0.05, ]

  
  
    
  # Examining ChIP Peaks ----------------------------------------------------
  
  ### --- Peaks in EITHER DMSO or Nutlin
  all_p53_bound <- union(dmso_peaks$gene, nutlin_peaks$gene)# get a list of gene names
  all_p53_bound_num <- length(all_p53_bound) 
  print(all_p53_bound_num)
 
  
  ### --- Peaks that are in genes that were measured
  bound_genes <- intersect(all_p53_bound, de_import$GeneID) # get a list of gene names 
  bound_genes_num <- length(bound_genes)
  print(bound_genes_num)

  print(length(bound_genes)) # print the number of genes in bound_de_genes
  print(bound_genes)
  
  ### --- Peaks that are in DE genes 
  bound_de_genes <- intersect(all_p53_bound, de$GeneID) # get a list of gene names 
  bound_de_genes_num <- length(bound_de_genes) # print the number of genes in bound_de_genes
  print(bound_de_genes_num)


### --- Create a venn diagram
  #install.packages("ggvenn") #install on this line. Once you install, comment this line out 
  library(ggvenn) # don't forget to load your new library on this line 
  venn_list <- list(
    dmso_peaks = dmso_peaks$gene, 
    nutlin_peaks = nutlin_peaks$gene,
    de_genes = de$GeneID
  ) # add everything that you want to compare into this list of lists 

  ggvenn(venn_list) #create venn diagram on this line!


  
  ### --- Export results for GO analysis   
  write.table(bound_de_genes, 
              paste0(outdir, 'genes_sig_p53peak.txt'), 
              col.names=FALSE, row.names=FALSE, 
              quote = FALSE)
  write.table(bound_genes, 
              paste0(outdir, 'genes_p53peak.txt'), 
              col.names=FALSE, row.names=FALSE, 
              quote = FALSE)
  

  
# Working with differential expression data -------------------------------
  
  ### --- Load differential expression analysis results
  
  hct_de <- read_tsv('data/deseq_output/hct116_deres.txt')
  hct_p53ko_de <- read_tsv('data/deseq_output/hct116p53ko_deres.txt')
  sjsa_de <- read_tsv('data/deseq_output/sjsa_deres.txt')
  mcf7_de <- read_tsv('data/deseq_output/mcf7_deres.txt')
  

# Make a venn diagram of just what is differential in each cell line?  

  hct_filt <- hct_de[(abs(hct_de$log2FoldChange) > 3 & hct_de$padj <= 0.05), ] 
  hct_p53ko_filt <- hct_p53ko_de[(abs(hct_p53ko_de$log2FoldChange) > 3 & hct_p53ko_de$padj <= 0.05), ] 
  sjsa_filt <- sjsa_de[(abs(sjsa_de$log2FoldChange) > 3 & sjsa_de$padj <= 0.05), ] 
  mcf7_filt <- mcf7_de[(abs(mcf7_de$log2FoldChange) > 3 & mcf7_de$padj <= 0.05), ] 


  
  de_venn <- list(
    HCT116 = hct_filt$GeneID, 
    `HCT116 p53-KO` = hct_p53ko_de$GeneID, 
    SJSA = sjsa_filt$GeneID, 
    MCF7 = mcf7_filt$GeneID
  )
  
  ggvenn(de_venn)
  

# Making a heatmap  -------------------------------------------------------

  #install.packages('pheatmap') # install heatmap on this line
  library(pheatmap) # load pheatmap on this line 
  
  ### --- FILTER datasets so that they are a reasonable size for a heatmap
  hct_filt <- hct_de[hct_de$padj <= 0.05, ]
  hct_p53ko_filt <- hct_p53ko_de[hct_p53ko_de$padj <= 0.05, ]
  sjsa_filt <- sjsa_de[sjsa_de$padj <= 0.05, ]
  mcf7_filt <- mcf7_de[mcf7_de$padj <= 0.05, ]
    
  ### --- Print dimensions below
  print('hct116')
  dim(hct_filt)
  print('hct116 p53ko')
  dim(hct_p53ko_filt)
  print('sjsa')
  dim(sjsa_filt)
  print('mcf7')
  dim(mcf7_filt)
    
  ### --- Rename columns so that the dataframes can be combined.
  colnames(hct_filt)[-1] <- paste('hct', colnames(hct_filt)[-1], sep = '_')
  colnames(hct_p53ko_filt)[-1] <- paste('hct_p53ko', colnames(hct_p53ko_filt)[-1], sep = '_')
  colnames(sjsa_filt)[-1] <- paste('sjsa', colnames(sjsa_filt)[-1], sep = '_')
  colnames(mcf7_filt)[-1] <- paste('mcf7', colnames(mcf7_filt)[-1], sep = '_')
  
  

# Join all of the DE data into one dataframe ------------------------------

  # Use one of the following functions 
  #     ~ full_join, inner_join, left_join, right_join   
  
  df <- full_join(hct_filt[, c(1,3)], hct_p53ko_filt[, c(1,3)])
  df <- full_join(df, sjsa_filt[, c(1,3)])
  df <- full_join(df, mcf7_filt[, c(1,3)])

  # Move GeneID to be the row names of the df and remove the GeneID column
  rownames(df) <- df$GeneID
  df <- df[,-c(1)]
  

# Filter to keep only rows with data in at least 2 cell lines -------------

  df$zero_count <- rowSums(is.na(df))
  df <- df[df$zero_count < 3, ]
  dim(df)
  
  df <- df[1:4]
  df[is.na(df)] <- 0
   

# Graph using pheatmap ----------------------------------------------------

  pheatmap(df) 
  
  breaks = seq(min(df), max(df), length.out = 52)
  pheatmap(df, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(52), 
           breaks = breaks,
           border_color = NA
  )
  

######### --------------------------------------------------------#########
######### --------------------------OPTIONAL----------------------#########
######### --------------------------------------------------------#########


# Challenge question 1. ---------------------------------------------------
nutlin_peaks_unique <- nutlin_peaks[, c('chip_start', 'chip_end', 'peak_id', 'gene')]
nutlin_peaks_unique <- nutlin_peaks_unique[!duplicated(nutlin_peaks_unique), ]
nutlin_peaks_unique$width <- nutlin_peaks_unique$chip_end - nutlin_peaks_unique$chip_start

dmso_peaks_unique <- dmso_peaks[, c('chip_start', 'chip_end', 'peak_id', 'gene')]
dmso_peaks_unique <- dmso_peaks_unique[!duplicated(dmso_peaks_unique), ]
dmso_peaks_unique$width <- dmso_peaks_unique$chip_end - dmso_peaks_unique$chip_start

peaks_widths <- rbind(nutlin_peaks_unique, dmso_peaks_unique)
peaks_widths$sample <- sapply(strsplit(peaks_widths$peak_id, '_'), `[`, 2)

#library(ggridges)
ggplot(peaks_widths, aes(y = sample, x = width, fill = sample)) + 
  #geom_violin()
  #geom_density_ridges(alpha = 0.5)
  geom_boxplot()

