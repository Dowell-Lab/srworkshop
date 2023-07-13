#' ---
#' title: "Counting Reads with featureCounts"
#' author: "Taylor Jones (edited by Rutendo Sigauke)"
#' date: "`r Sys.Date()`"
#' output: 
#'   html_document: 
#'     toc: yes
#'     toc_float:
#'       collapsed: yes
#'       smooth_scroll: yes
#' ---

#' # Introduction
#' 
#' Here we will learn how to download a package, what metadata table is (and why it is important), and run featureCounts, which counts reads over genes.
#' 
#' # We will want to start fresh and clear our environment.
#' 
#' Start by clearing your console. To do this hit `Ctrl+l` or go to Edit --> Clear Console. Clear your environment and plots by hitting the broom icon in both those cells reset our working directory.
#' (NB: You can also open the `d6_featureCounts.R` in vim to make edits)
## ----setwd, eval=TRUE--------------------------------------------------------------------------------------------

workdir <- '/PATH/TO/WORKING/DIRECTORY'
setwd(workdir)
getwd()

#' # Install packages 
#' 
#' Like before: most packages you can install with this syntax: `install.packages("PACKAGE")`. Such as this package: `install.packages("tidyverse")` --great package!
#' HOWEVER, some of the bioinformatic software is located in Bioconductor, including `RSubread`(R wrapper for [**Subread**](http://subread.sourceforge.net/)) which contains `featureCounts`.

#' ## Install `BiocManager`
#' 
#' Note: `eval=FALSE` means the code in the chunk is NOT evaluated. To evaluate the code chunk, change `eval=FALSE` to `eval=TRUE`.
#' The installations only needs to be run once. Turn back `eval=FALSE` after installation.

## ----installBiocmanager, eval=FALSE------------------------------------------------------------------------------
## ##This only needs to be run once. Turn back eval=FALSE after installation.
## 
## if (!requireNamespace("BiocManager", quietly = TRUE))  ## Install BiocManager
##   install.packages("BiocManager")  #this will take a moment to run


#' ## Install `RSubread`.
#' 
## ----installRsubread, eval=FALSE---------------------------------------------------------------------------------
## BiocManager::install("Rsubread")  #This actually installs RSubread

#' ## Load the `Rsubread` package
#' This actually loads the library into the environment
#' 
## ----loadRsubread, eval=TRUE-------------------------------------------------------------------------------------

library("Rsubread")

#' Here, this calls for help of a package. This is useful if you have no idea where to start.
#' 
## ----rsubreadHelp, eval=FALSE------------------------------------------------------------------------------------
## 
## ??Rsubread # click on the Help pages search results

## ----featureCountsHelp, eval=FALSE-------------------------------------------------------------------------------
## 
## ??Rsubread::featureCounts # we can look specifically at featureCounts and the featureCount flags
 
#' # Counting reads with `Rsubread`
#' 
#' NOTE: We run `featureCounts` on a cluster because bam files are LARGE.
#' We want to run this on the cluster but using command line `R` is not always great.
#' An alternative is to generate a script like this then submit this script to the cluster in an sbatch script.
#' 
#' ## Inputs for `featureCounts`
#' 
#' ### Bam files
#' 
#' Read in bam file list
## ----bamdir, eval=TRUE-------------------------------------------------------------------------------------------
bamdir <- '/scratch/Shares/public/sread/data_files/day6/bam'

filelist <- list.files(path=bamdir,
                       pattern="sorted.bam$",
                       full.names=TRUE)

outdir <- paste(workdir,'/', 'counts', '/', sep='') ##naming our outdir
dir.create(outdir) ###creating the directory

#' ### Annotation file in GTF format 
#' 
#' We use a GTF -- this is necessary for `featureCounts` and more useful for RNA-seq data as it is cleaner/faster for counting over exons
#' 
## ----loadGTF, eval=TRUE------------------------------------------------------------------------------------------

hg38gtf <- "/scratch/Shares/sread/data_files/day7/annotations/hg38_ucsc_genes_chr21.gtf"

#' ## Read counting using `featureCounts`
#' 
#' These parameters below are for specifying how to count over the gtf file.
#' 
## ----featurecounts, eval=TRUE------------------------------------------------------------------------------------
# Sink will save out stdout -- check all of these settings (e.g. isPairedEnd)!

# Define root name for output files
fileroot <- "chr21_Ethan_Eric"

sink(paste0(outdir, fileroot, "_featureCounts_gene_rnaseq.txt")) #we want to save the output as a file
fc <- featureCounts(files=filelist,
                    annot.ext=hg38gtf,
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="exon",
                    GTF.attrType="gene_id",
                    useMetaFeatures=TRUE,
                    allowMultiOverlap=TRUE,
                    largestOverlap=TRUE,
                    countMultiMappingReads=TRUE,
                    isPairedEnd=TRUE,
		    strandSpecific=2,
                    nthreads=4)  #when you move to a bigger machine change to 8

#Flag info here: ??Rsubread::featureCounts for specifics on each flag

sink()

#' ## Write results 
#' 
#' We'll keep our accession number (GeneID) and length for filtering genes in DESeq2. 
#' We will also be writing different files to separate files (listed below).
#' - `featureCounts_gene_rnaseq.txt` : GeneID, Length, Counts
#' - `coverage.csv` : Counts
#' - `.stat.csv` : Coverage Statistics
#' - `.annotation.csv` : GeneID, Chromosome, Start, End, Strand, Length
#' - `.targets.csv` : file names for input `bam` files
#' 
## ----writecounts, eval=TRUE--------------------------------------------------------------------------------------
# Write output file with write.table

write.table(x=data.frame(fc$annotation[,c("GeneID","Length")],
                         fc$counts,stringsAsFactors=FALSE),
            paste0(outdir, fileroot, file="_featureCounts_gene_rnaseq.txt"),
            quote=FALSE,sep="\t",
            row.names=FALSE)

# Write out separate files for counts and other count statistics
# These files are written as comma separated files with write.csv

write.csv(fc$counts, paste(outdir, fileroot,".coverage.csv", sep=""))
write.csv(fc$stat, paste(outdir, fileroot,".stat.csv", sep=""))
write.csv(fc$annotation, paste(outdir, fileroot,".annotation.csv", sep=""))
write.csv(fc$targets, paste(outdir, fileroot,".targets.csv", sep=""))
