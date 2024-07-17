


# Set environment ---------------------------------------------------------

  library(ggplot2) # load ggplot2 
  setwd('your/path/to/2024_shortread/day9/') # set your local working directory 
  outdir <- 'results/'
  

# Read in data ------------------------------------------------------------

  ### --- READ IN CHIP DATA 
  dmso_peaks <- read.table('bedtools_results/p53_peaks_in_genes_hct_dmso.bed', 
                           header = FALSE) 
  nutlin_peaks <- read.table('bedtools_results/p53_peaks_in_genes_hct_nutlin.bed', 
                             header = FALSE)
  
  ### --- GIVE CHIP DATA USEFUL COLUMN NAMES
  set_colnames <- c('chip_chr', 'chip_start', 'chip_end', 'peak_id', 'score', 
                    'chip_strand', 'signalValue', 'log_pval', 'log_qval', 
                    'summit', 'promoter_chr', 'promoter_start', 'promoter_end', 
                    'gene', '.', 'gene_strand', 'overlap') 
  colnames(dmso_peaks) <- set_colnames
  colnames(nutlin_peaks) <- set_colnames
  
  ### --- READ IN DE DATA 
  de.import <- read.table('data/deseq_output/hct116_deres.txt')# read in data here 
  
  ### --- EXTRACT SIGNIFICANTLY CHANGING GENES
  de <- de.import[abs(de.import$log2FoldChange) > 2 & de.import$padj < 0.05, ]

  
  
    
  # Examining ChIP Peaks ----------------------------------------------------
  
  ### --- Peaks in EITHER DMSO or Nutlin
  all_peaks <- union(dmso_peaks$gene, nutlin_peaks$gene)# get a list of gene names 
  print(length(all_peaks)) # print the number of genes in all_peaks
  
  ### --- Peaks that are in genes that were measured
  bound_genes <- intersect(all_peaks, de.import$GeneID) # get a list of gene names 
  print(length(bound_genes)) # print the number of genes in bound_de_genes
  print(bound_genes)
  
  ### --- Peaks that are in DE genes 
  bound_de_genes <- intersect(all_peaks, de$GeneID) # get a list of gene names 
  print(length(bound_de_genes)) # print the number of genes in bound_de_genes
  print(bound_de_genes)
  

  ### --- Create a venn diagram

  #install.packages("ggvenn") #install on this line. Once you install, comment this line out 
  library(ggvenn) # don't forget to load your new library on this line 
  venn_list <- list(
    dmso_peaks = dmso_peaks$gene, 
    nutlin_peaks = nutlin_peaks$gene,
    de_genes = de$GeneID
  ) # add everything that you want to compare into this list of lists 

  ggvenn(venn_list) #create venn diagram on this line!


  ### --- Calculate Percentage
  percent <- (19 + 64)/length(all_peaks) * 100
  print(percent)
  
  
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
  
  hct.de <- read.table('data/deseq_output/hct116_deres.txt')
  hct_p53ko.de <- read.table('data/deseq_output/hct116p53ko_deres.txt')
  sjsa.de <- read.table('data/deseq_output/sjsa_deres.txt')
  mcf7.de <- read.table('data/deseq_output/mcf7_deres.txt')
  
  
  ### --- Create a venn diagram of the cell lines 
  
  de_venn <- list(
    HCT116 = hct.de$GeneID, 
    `HCT116 p53-KO` = hct_p53ko.de$GeneID, 
    SJSA = sjsa.de$GeneID, 
    MCF7 = mcf7.de$GeneID
  )
  
  ggvenn(de_venn)
  
  # consider modifying the venn diagram below 
  # more helpful version of the venn diagram 
  hct.filt <- hct.de[(hct.de$padj <= 0.05), ] #abs(hct.de$log2FoldChange) > 3 & 
  hct_p53ko.de.filt <- hct_p53ko.de[(hct_p53ko.de$padj <= 0.05), ] #abs(hctko.de$log2FoldChange) > 3 & 
  sjsa.filt <- sjsa.de[(sjsa.de$padj <= 0.05), ] #abs(sjsa.de$log2FoldChange) > 3 & 
  mcf7.filt <- mcf7.de[(mcf7.de$padj <= 0.05), ] #abs(mcf7.de$log2FoldChange) > 3 & 
  
  de_venn <- list(
    HCT116 = hct.filt$GeneID, 
    #`HCT116 p53-KO` = hct_p53ko.de$GeneID, 
    SJSA = sjsa.filt$GeneID, 
    MCF7 = mcf7.filt$GeneID
  )
  
  ggvenn(de_venn)
  

# Making a heatmap  -------------------------------------------------------

  #install.packages('pheatmap') # install heatmap on this line
  library(pheatmap) # load pheatmap on this line 
  
  ### --- FILTER datasets so that they are a reasonable size for a heatmap
  hct.filt <- hct.de[(abs(hct.de$log2FoldChange) > 3 & hct.de$padj <= 0.05), ]
  hct_p53ko.filt <- hct_p53ko.de[(abs(hct_p53ko.de$log2FoldChange) > 3 & hct_p53ko.de$padj <= 0.05), ]
  sjsa.filt <- sjsa.de[(abs(sjsa.de$log2FoldChange) > 3 & sjsa.de$padj <= 0.05), ]
  mcf7.filt <- mcf7.de[(abs(mcf7.de$log2FoldChange) > 3 & mcf7.de$padj <= 0.05), ]
    
  ### --- Print dimensions below
  print('hct116')
  dim(hct.filt)
  print('hct116 p53ko')
  dim(hct_p53ko.filt)
  print('sjsa')
  dim(sjsa.filt)
  print('mcf7')
  dim(mcf7.filt)
    
  ### --- Rename columns so that the dataframes can be combined.
  colnames(hct.filt)[-1] <- paste('hct', colnames(hct.filt)[-1], sep = '_')
  colnames(hct_p53ko.filt)[-1] <- paste('hct_p53ko', colnames(hct_p53ko.filt)[-1], sep = '_')
  colnames(sjsa.filt)[-1] <- paste('sjsa', colnames(sjsa.filt)[-1], sep = '_')
  colnames(mcf7.filt)[-1] <- paste('mcf7', colnames(mcf7.filt)[-1], sep = '_')
  
  

# Join all of the DE data into one dataframe ------------------------------

  # Use one of the following functions 
  #     ~ full_join, inner_join, left_join, right_join   
  
  df <- full_join(hct.filt[, c(1,3)], hct_p53ko.filt[, c(1,3)])
  df <- full_join(df, sjsa.filt[, c(1,3)])
  df <- full_join(df, mcf7.filt[, c(1,3)])

  df <- df %>% column_to_rownames('GeneID')
  

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

