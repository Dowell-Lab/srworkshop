


# Set environment ---------------------------------------------------------

  library() # load ggplot2 
  setwd('your/path/to/<your short read dir>/day9/') # set your local working directory 
  outdir <- 'results/'
  

# Read in data ------------------------------------------------------------

  ### --- READ IN CHIP DATA 
  dmso_peaks <- read_tsv('path', <options>) 
  nutlin_peaks <- read_tsv('path', <options>)
  
  ### --- GIVE CHIP DATA USEFUL COLUMN NAMES
  set_colnames <- c('chip_chr', 'chip_start', 'chip_end', 'peak_id', 'score', 
                    'chip_strand', 'signalValue', 'log_pval', 'log_qval', 
                    'summit', 'gene_chr', 'gene_start', 'gene_end', 
                    'gene', '.', 'gene_strand', 'overlap') 
  colnames(dmso_peaks) <- set_colnames
  colnames(nutlin_peaks) <- set_colnames
  
  ### --- READ IN DE DATA 
  de_import <- read_tsv('path')# read in data here 
  
  ### --- EXTRACT SIGNIFICANTLY CHANGING GENES
  de <- de_import[de_import$"<column_name>" < "<filter_conditions>", ]

  
  
    
  # Examining ChIP Peaks ----------------------------------------------------
  
  ### --- Peaks in EITHER DMSO or Nutlin
  all_p53_bound <- 
  all_p53_bound_num <- 
  print(all_p53_bound_num)
 
  
  ### --- Peaks that are in genes that were measured
  bound_genes <- 
  bound_genes_num <- 
  print(bound_genes_num)

  
  ### --- Peaks that are in DE genes 
  bound_de_genes <- 
  bound_de_genes_num <- 
  print(bound_de_genes_num)


### --- Create a venn diagram
  #install.packages("ggvenn") #install on this line. Once you install, comment this line out 
  library(ggvenn) # don't forget to load your new library on this line 
  venn_list <- list( ) # add everything that you want to compare into this list of lists 

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
  
  hct_de <- read_tsv('<path>')
  hct_p53ko_de <- read_tsv('<path>')
  sjsa_de <- read_tsv('<path>')
  mcf7_de <- read_tsv('<path>')
  

# Make a venn diagram of just what is differential in each cell line?  

  hct_filt <- hct_de[<filter_conditions>, ] 
  hct_p53ko_filt <- hct_p53ko_de[<filter_conditions> , ] 
  sjsa_filt <- sjsa_de[<filter_conditions>, ] 
  mcf7_filt <- mcf7_de[<filter_conditions>, ] 


  
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
  
  df <- # join first two df here
  df <- # join df with the third df here
  df <- # join df with the fourth df here

  # Move GeneID to be the row names of the df and remove the GeneID column
  df <- as.data.frame(df)
  rownames(df) <- df$GeneID
  df <- df[,-c(1)]
  

# Filter to keep only rows with data in at least 2 cell lines -------------

  df$zero_count <- rowSums(is.na(df))
  df <- df[df$zero_count "<operator>" "<put what filtering by here>", ]
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
nutlin_peaks_unique <- nutlin_peaks[]
nutlin_peaks_unique <- 
nutlin_peaks_unique$width <- 

dmso_peaks_unique <- dmso_peaks[]
dmso_peaks_unique <- 
dmso_peaks_unique$width <- 

peaks_widths <- rbind()
peaks_widths$sample <- sapply(strsplit(peaks_widths$peak_id, '_'), `[`, 2)


ggplot(peaks_widths, aes(y = sample, x = width, fill = sample)) + 
  geom_boxplot()

