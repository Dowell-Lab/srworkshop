# Day 6 | Introduction to R

This workshop gives an introduction to R! Briefly, R is a statistical programming language with a wide array of applications. In this workshop you will learn how to run R code in a few environments. You will learn basic R commands, how to load tables, save tables, plot publication ready figures and save all these files.

## Inclass Worksheets

We will cover running R three ways:

1. In an R console

- `Day6_worksheet1_Introduction_to_R.md`  

2. In RStudio 

- `Day6_worksheet2_R_in_Rstudio.md` and `Learning_R.R` will guide you

3. Submitting and R script as an sbatch job

- `Learning_R_submit_aws.R` R script
- `Submit_Rscript.sbatch` sbatch script calling the `Learning_R_submit_aws.R` R script

4. R challenge for really advanced students:

-   Find some data table online or your own data table.
-   Copy the data to the supercomputer.
-   Could you create an R script to graph the data for one subsection (i.e. graph data for one person, or one gene or one deer, or one sports team)?
-   Create a sbatch script that runs the R script.
-   Graph every subsection via the R script on separate CPUs.
-   An example which will not work because you don't have tidyverse installed is under scripts and called step1. 
   

## Homework

1. Homework can be done in R Studio

- `Learning_R_Additional_Practice.R`

2. Additional homework is to install R packages we will use later in the week. Both packages can be found on bioconductor.

Project B : Multi-omics (RNA-seq & ChIP-seq)

More details in **Project B Day07 Worksheets** folder 

- [rsubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html) : Install on AWS R
- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) : Install on personal computer R

