#' ---
#' title: "Counting Reads with featureCounts"
#' author: "Taylor Jones (edited by Rutendo Sigauke, Samuel Hunter, Lynn Sanford, Malia Fredrickson)"
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
#' Here we will learn how to install and download a package, what metadata table is (and why it is important), and run featureCounts, which counts reads over genes.

# Set working directory
workdir <- '/PATH/TO/WORKING/DIRECTORY'
setwd(workdir)
getwd()

# install Rsubread if needed and load
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager",
                     repos = "https://cloud.r-project.org")
}

if (!requireNamespace("Rsubread", quietly = TRUE)) {
    BiocManager::install("Rsubread", ask = FALSE, update = FALSE)
}

library(Rsubread)

#' The following commands call for help pages of a package. This is useful if you have no idea where to start.

## ??Rsubread # click on the Help pages search results
## ??Rsubread::featureCounts # we can look specifically at featureCounts and the featureCount flags
 
#' Counting reads with `Rsubread`
#' 
#' Inputs for `featureCounts` are BAM files. We run `featureCounts` on a cluster because bam files are LARGE.

# Read in bam file list
bamdir <- '/scratch/Shares/public/sread/data_files/day7/bam'

filelist <- list.files(path=bamdir,
                       pattern="sorted.bam$",
                       full.names=TRUE)

outdir <- paste(workdir,'/', 'counts', '/', sep='') ## naming our outdir
dir.create(outdir) ## creating the directory

#' Annotation file in GTF format 
#' 
#' We use a GTF -- this is necessary for `featureCounts` and more useful for RNA-seq data as it is fast and clean for counting over exons

hg38gtf <- "/scratch/Shares/public/sread/data_files/day7/annotations/hg38_ucsc_genes_chr21.gtf"

#' Read counting using `featureCounts`

fc <- featureCounts(
    files=filelist,
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
    nthreads=2 ## Make sure this number matches the sbatch script ntasks
)  

#' Again, look at ??Rsubread::featureCounts for specifics on each flag

#' Write results 
#' 
#' We'll keep our accession number (GeneID) and length for filtering genes in DESeq2. 
#' We will also be writing different files to separate files (listed below).
#' - `featureCounts_gene_rnaseq.txt` : GeneID, Length, Counts
#' - `coverage.csv` : Counts
#' - `.stat.csv` : Coverage Statistics
#' - `.annotation.csv` : GeneID, Chromosome, Start, End, Strand, Length
#' - `.targets.csv` : file names for input `bam` files

# Define root name for output files
fileroot <- "chr21_Ethan_Eric"

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
