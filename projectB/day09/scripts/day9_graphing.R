


# Set environment ---------------------------------------------------------

  library() # load tidyverse 
  setwd() # set your local working directory 

# Read in data ------------------------------------------------------------

  ### --- READ IN CHIP DATA 
  dmso_peaks # read in data here, remember that bed files don't have column names     
  nutlin_peaks # read in data here 
  
  ### --- GIVE CHIP DATA USEFUL COLUMN NAMES
  
  # The "set_colnames" line is already written correctly to rename you columns
  set_colnames <- c('chip_chr', 'chip_start', 'chip_end', 'peak_id', 'score', 
                    'chip_strand', 'signalValue', 'log_pval', 'log_qval', 
                    'summit', 'promoter_chr', 'promoter_start', 'promoter_end', 
                    'gene', '.', 'gene_strand', 'overlap') 
  colnames(dmso_peaks) <- set_colnames
  colnames(nutlin_peaks) <- set_colnames
  
  ### --- READ IN DIFFERENTIAL EXPRESSION DATA 
  de.import # read in data here 
  
  ### --- EXTRACT SIGNIFICANTLY CHANGING GENES
  de <- # code to filter for significant genes 
  
    
    
  
# Examining ChIP Peaks ----------------------------------------------------

  ### --- Which peaks are in EITHER DMSO or Nutlin?
  all_peaks <- # get a list of gene names 
  print() # print the number of genes in all_peaks
  
  ### --- Which peaks are in DE genes?
  bound_de_genes <- # get a list of gene names 
  print() # print the number of genes in bound_de_genes
  

  ### --- Create a venn diagram
  #install on this line. Once you install, comment this line out 
  library() # don't forget to load your new library on this line 
  venn_list <- list() # add everything that you want to compare into this list of lists 
  #create venn diagram on this line!


  ### --- Export results for GO analysis  
  background <- # code here
  all_de <- # code here
  up_de <- # code here
  down_de <- # code here
    
  write_tsv(background, )
  write_tsv(all_de, )
  write_tsv(up_de, )
  write_tsv(down_de, )
  

  

# Working with differential expression data -------------------------------

  ### --- Load differential expression analysis results

  hct.de <- read_tsv('')
  hctko.de <- read_tsv('')
  sjsa.de <- read_tsv('')
  mcf7.de <- read_tsv('')
  
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

# Filter to keep only rows with data in at least 2 cell lines -------------
  # Hint, use rowsums() and is.na() 
  

# Graph using pheatmap ----------------------------------------------------

  pheatmap()
  
  