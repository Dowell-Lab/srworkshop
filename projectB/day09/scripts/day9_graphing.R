


# Set environment ---------------------------------------------------------

  library() # load tidyverse 
  setwd() # set your local working directory 
  
  

# Read in data ------------------------------------------------------------

  ### --- READ IN CHIP DATA 
  dmso_peaks # read in data here     
  nutlin_peaks # read in data here 
  
  ### --- GIVE CHIP DATA USEFUL COLUMN NAMES
  set_colnames <- c('chip_chr', 'chip_start', 'chip_end', 'peak_id', 'score', 
                    'chip_strand', 'signalValue', 'log_pval', 'log_qval', 
                    'summit', 'promoter_chr', 'promoter_start', 'promoter_end', 
                    'gene', '.', 'gene_strand', 'overlap') 
  colnames(dmso_peaks) <- set_colnames
  colnames(nutlin_peaks) <- set_colnames
  
  ### --- READ IN DE DATA 
  de.import # read in data here 
  
  ### --- EXTRACT SIGNIFICANTLY CHANGING GENES
  
  
# First, the set of genes bound by p53 in EITHER DMSO or Nutlin -----------

  all_peaks <- # get a list of gene names 
  print() # print the number of genes in all_peaks
  

# Second, the intersection of p53 peaks and DE genes ----------------------

  bound_de_genes <- # get a list of gene names 
  print() # print the number of genes in bound_de_genes
  
  

# Create a venn diagram  --------------------------------------------------

  #install on this line. Once you install, comment this line out 
  library() # don't forget to load your new library on this line 
  venn_list <- list() # add everything that you want to compare into this list of lists 
  #create venn diagram on this line!

  

# Export results for GO analysis  -----------------------------------------

  background <- # code here
  all_de <- # code here
  up_de <- # code here
  down_de <- # code here
    
  write_tsv(background, )
  write_tsv(all_de, )
  write_tsv(up_de, )
  write_tsv(down_de, )
  
  
  
  