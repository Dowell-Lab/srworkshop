library(tidyverse)
library(data.table)              

indir = "/Users/allenma/simulatedata/"

gtf_path <- paste0(indir, "gencode.v27.chr_patch_hapl_scaff.annotation.gtf.gz")

# Read GTF (skip comment lines)
gtf <- fread(
  gtf_path,
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  col.names = c("seqname","source","feature","start","end","score","strand","frame","attr")
)

# Extract gene_id and gene_name from column 9
gene_info <- gtf %>%
  mutate(
    gene_id   = str_match(attr, 'gene_id "([^"]+)"')[,2],
    gene_name = str_match(attr, 'gene_name "([^"]+)"')[,2]
  ) %>%
  filter(!is.na(gene_id)) %>%
  distinct(gene_id, gene_name, seqname)

head(gene_info)


fn=paste0(indir, "genotype.csv")
genotype <- read.csv(fn, row.names = 1)
head(genotype)

fn=paste0(indir, "expression.csv")
RNAlevel <- read.csv(fn, row.names = 1)
head(RNAlevel)

#what are the steps to making a simulated person with T21 from the D21 invidials?











#here is one path you could use
#step 1, make a gene_info_chr21 dataframe from gene_info but it only has chr21 genes
#step 2, a genotype_D21 dataframe that is genotype but only has people whose genotype is D21
#step 3, make a RNAlevel_D21 that is RNA_level but only has D21 invidials (hint, uses the Dataframe from step 2)
#step 4, make a RNAlevel_simulated_T21 that starts as just a copy of RNAlevel_D21
#then usse the function mutate to by multiple expression by 1.5 if the gene is on chr 21
#this will use the dataframe from step 1

#How will you test if it worked? Plot one gene from both D21 indivals and simualted T21 individals. 

