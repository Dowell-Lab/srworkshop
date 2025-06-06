


# Set environment ---------------------------------------------------------

  library() # load ggplot2
  setwd() # set your local working directory
  outdir <- # define an output directory

# Read in data ------------------------------------------------------------

  ### --- READ IN CHIP DATA 
  dmso_peaks <- # read in data here, remember that bed files don't have column names     
  nutlin_peaks <- # read in data here 
  
  ### --- GIVE CHIP DATA USEFUL COLUMN NAMES
  
  # The "set_colnames" line is already written correctly to rename you columns
  set_colnames <- c('chip_chr', 'chip_start', 'chip_end', 'peak_id', 'score', 
                    'chip_strand', 'signalValue', 'log_pval', 'log_qval', 
                    'summit', 'promoter_chr', 'promoter_start', 'promoter_end', 
                    'gene', '.', 'gene_strand', 'overlap') 
  colnames(dmso_peaks) <- set_colnames
  colnames(nutlin_peaks) <- set_colnames
  
  ### --- READ IN DIFFERENTIAL EXPRESSION DATA 
  de.import <- # read in data here 
  
  ### --- EXTRACT SIGNIFICANTLY CHANGING GENES
  de <- # code to filter for significant genes 
  
    
    
  
# Examining ChIP Peaks ----------------------------------------------------

  ### --- Peaks in EITHER DMSO or Nutlin
  all_p53_bound <- # get a list of gene names 
  print() # print the number of genes in all_p53_bound
  
  
  ### --- Peaks that are in genes that were measured by DESeq2 at all
  bound_genes <- intersect() # get a list of gene names 
  print() # print the number of genes in bound_genes
  
  
  ### --- Peaks that are in DE genes 
  bound_de_genes <- intersect() # get a list of gene names 
  print() # print the number of genes in bound_de_genes
  print() # print the list of gene IDs
    
  ### --- Create a venn diagram
  #install on this line. Once you install, comment this line out 
  library() # don't forget to load your new library on this line 
  venn_list <- list() # add everything that you want to compare into this list of lists 
  #create venn diagram on this line!

  ### --- Calculate Percentage


  ### --- Export results for GO analysis  
  write.table(bound_de_genes, 
              paste0(outdir, 'genes_sig_p53peak.txt'), 
              )
  write.table(bound_genes, 
              paste0(outdir, 'genes_p53peak.txt'), 
              )
  
  

# Working with differential expression data -------------------------------

  ### --- Load differential expression analysis results

  hct.de <- read.table('')
  hct_p53ko.de <- read.table('')
  sjsa.de <- read.table('')
  mcf7.de <- read.table('')
  
  
  
  ### --- Create a venn diagram of the cell lines 

  de_venn <- list(  )
  ggvenn(de_venn)
  
  # consider modifying the venn diagram below 
  

# Making a heatmap  -------------------------------------------------------

  # install heatmap on this line 
  # load pheatmap on this line 
  
  ### --- FILTER datasets so that they are a reasonable size for a heatmap
  hct.filt <- 
  hct_p53ko.filt <- 
  sjsa.filt <- 
  mcf7.filt <- 
    
  ### --- Print dimensions below
    
    
  ### --- Rename columns so that the dataframes can be combined.
  colnames(hct.filt)[-1] <- paste('hct', colnames(hct.filt)[-1], sep = '_')
  colnames(hct_p53ko.filt)[-1] <- paste('hct_p53ko', colnames(hct_p53ko.filt)[-1], sep = '_')
  colnames(sjsa.filt)[-1] <- paste('sjsa', colnames(sjsa.filt)[-1], sep = '_')
  colnames(mcf7.filt)[-1] <- paste('mcf7', colnames(mcf7.filt)[-1], sep = '_')
  
  

# Join all of the DE data into one dataframe ------------------------------

  # Use one of the following functions 
  #     ~ full_join, inner_join, left_join, right_join   
  
  df < - # name the combined data frame just df
  # Move GeneID to be the row names of the df and remove the GeneID column

# Filter to keep only rows with data in at least 2 cell lines -------------
  # Hint, use rowsums() and is.na() 
  # look up how to replace "NA" with 0 with Google 
  

# Graph using pheatmap ----------------------------------------------------

  pheatmap()
  

######### --------------------------------------------------------#########
######### --------------------------OPTIONAL----------------------#########
######### --------------------------------------------------------#########


# Challenge question 1. ---------------------------------------------------

  nutlin_peaks_unique <- # select columns mentioned in instructions
  nutlin_peaks_unique <- # remove duplicated rows 
  nutlin_peaks_unique$width <- # calculate width of the ChIP peak 
  
  dmso_peaks_unique <- # select columns mentioned in instructions
  dmso_peaks_unique <- # remove duplicated rows 
  dmso_peaks_unique$width <- # calculate width of the ChIP peak 
  
  peaks_widths <- # combine the two dataframes you made above
  # BELOW is the code to create a new column with just the sample ID in it 
  # You NEED this column to plot (that's a hint ;) ) 
  peaks_widths$sample <- sapply(strsplit(peaks_widths$peak_id, '_'), `[`, 2)
  
  
  # Use ggplot to create a boxplot, violin plot, or ridge plot
  