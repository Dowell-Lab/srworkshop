## Learning R : example for submitting in sbatch
## Author: Taylor Jones (2022), Rutendo Sigauke (2023)

###############################################################################
## In this script we plan to submit the R script as an SBATCH job            ##
## This is useful for when you have a computationally intensive job that     ##
## requires more resources than on your personal computer.                   ##
## For example, counting reads in a bam file over gene features for RNA-seq. ##
###############################################################################

##############################################################################
# Set your working directory first
workdir <- '/Users/maallen3/srworkshop/day06/scripts/'
outdir <- workdir
setwd(workdir)
getwd()

##############################################################################
## In the Learning_R.R script you processed data(iris) and data(mtcars)     ##
## NOW, we are processing data(mtcars) data in a script called              ##
## Learning_R_inclass.R                                                     ##
##############################################################################

## NOTE: Submitting R job as a script is NOT interactive
## Therefore, you have to save all outputs in order to view them

args = commandArgs(trailingOnly=TRUE)

# 0: load packages

library(tidyverse)
library(ggplot2)
library(ggsave)

# 1: load in the data for all genes

pathtoindata<-"/scratch/Shares/public/sread2024/cookingShow/day6/populationRNAseq/"
genesymbolchr <-paste0(pathtoindata,"genesymbolchr.csv")
ncdf_long<-"D21_T21_RNA_long.csv"

# 2: find this gene

myfavoritegene<-"RUNX1"

# 3: Create a simple violin plot of this gene in D21 and T21 individuals

thisgeneEnsmbl <- genesymbolchr[genesymbolchr$SYMBOL==myfavoritegene,] %>% drop_na()
print(thisgeneEnsmbl)

plot_sample_types <- c("D21", "complete_T21", "mosaic_trisomy_21")

for (ENSEMBLname in thisgeneEnsmbl$ENSEMBL){
  print(ENSEMBLname)
  thisgene_ncdf_long <- ncdf_long %>% dplyr::filter(ENSEMBL==as.character(ENSEMBLname))
  
  thisgene_ncdf_long <- thisgene_ncdf_long %>% filter(sample_type %in% plot_sample_types)
  
  thisgene_ncdf_long <- merge(thisgene_ncdf_long, acomorbidf, on=Patient)
  
  titleoffileandgraph = paste(myfavoritegene,"_", ENSEMBLname,"_",whichoutput)
  
  p <- ggplot(thisgene_ncdf_long, aes(x=sample_type, y=normalized_counts, group=interaction(sample_type, comorbidity), color=comorbidity))+
    geom_violin() + 
    geom_boxplot(outlier.shape = NA)+ 
    geom_point(position = position_jitterdodge(), alpha=0.3) +
    ggtitle(paste0(myfavoritegene, " and ", acomorbid))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()
  
  print(p)
  
  fn = paste(outdir,myfavoritegene,"_",ENSEMBLname,"_and_", acomorbid, ".png", sep="")
  print(fn)
  ggsave(filename = fn, plot=p)
}


# 4: Optional- create and save another pretty plot for mtcars data below

