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
workdir <- "/scratch/Users/maallen3/workshop-day6/"
indir <- '/scratch/Shares/public/sread2025/cookingShow/day6/populationRNAseq/'
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

#args = commandArgs(trailingOnly=TRUE) 

# 0: load packages

library(ggplot2)

# 1: load in the data for all genes

genesymbolchr <-paste0(indir,"genesymbolchr.csv")
ncdf_long<-paste0(indir,"D21_T21_RNA_long.csv")
genesymbolchrdf <- read.csv(genesymbolchr)
ncdf_longdf <- read.csv(ncdf_long)



# 2: find this gene

#myfavoritegene<-"RUNX1"
args = commandArgs(trailingOnly=TRUE) 

# 3: Create a simple violin plot of this gene in D21 and T21 individuals


thisgeneEnsmbl <- genesymbolchrdf[genesymbolchrdf$SYMBOL==myfavoritegene,] 

plot_sample_types <- c("D21", "complete_T21", "mosaic_trisomy_21")

for (ENSEMBLname in thisgeneEnsmbl$ENSEMBL){
  print(ENSEMBLname)
  ncdf_longdf_thisgene <- ncdf_longdf[ncdf_longdf$ENSEMBL==ENSEMBLname,] 
  titleoffileandgraph = paste(myfavoritegene,"_", ENSEMBLname)
  
  p <- ggplot(ncdf_longdf_thisgene, aes(x=sample_type, y=normalized_counts))+
    geom_violin() + 
    geom_boxplot(outlier.shape = NA)+ 
    geom_point(aes(colour="name"),position = position_jitterdodge(), alpha=0.3) +
    ggtitle(titleoffileandgraph)+
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()
  print(p)
  
  fn = paste(outdir,myfavoritegene,"_",ENSEMBLname,".png", sep="")
  print(fn)
  ggsave(fn, plot = last_plot(), width = 6, height = 4)
}



