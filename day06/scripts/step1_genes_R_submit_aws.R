## Learning R : example for submitting in sbatch
## Author: Taylor Jones (2022), Rutendo Sigauke (2023), Mary Allen (2025)

# Set your working directory first
workdir <- "/scratch/Users/maallen3/workshop-day6/"
indir <- '/scratch/Shares/public/sread2025/cookingShow/day6/populationRNAseq/'
outdir <- workdir
setwd(workdir)
getwd()

# 0: load packages

library(ggplot2)

# 1: load in the data for all genes

genesymbolchr <-paste0(indir,"genesymbolchr.csv")
ncdf_long<-paste0(indir,"D21_T21_RNA_long.csv")
genesymbolchrdf <- read.csv(genesymbolchr)
ncdf_longdf <- read.csv(ncdf_long)


# 2: find this gene

#myfavoritegene<-"RUNX1" #this was my orginal test of the script
#then I moved to using args so I can pull in the gene name from the sbatch script
args = commandArgs(trailingOnly=TRUE)
myfavoritegene = args[1]

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
