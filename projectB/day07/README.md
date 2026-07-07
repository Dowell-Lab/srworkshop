# Project B Day 7 | RNA-seq: Counting Reads and Differential Expression

Today, you will learn to count reads using *featureCounts*. Then, you will learn about exploratory data analysis and QC checks associated with differential analysis. Finally, you will use the differential analysis pipeline *DESeq2* to identify a list of genes which are differentially expressed in our dataset. This list will be used later in the week when we integrate it together with our ChIP-seq results.

## Before Day 7
- Please watch the following videos:
  - <a href="https://www.youtube.com/watch?v=nQDpoM2vc8k" target="_blank">B7.1 | Counting Reads</a>
  - <a href="https://www.youtube.com/watch?v=Tk4Q91Kvr7o" target="_blank">B7.2 | Differential Expression Analysis</a>
  - <a href="https://www.youtube.com/watch?v=JfYo1eXSbtg" target="_blank">B7.3 | DESeq2</a>
  - <a href="https://www.youtube.com/watch?v=X6p3E-QTcUc" target="_blank">B7.4 | Multifactor Designs in DESeq2 (Optional)</a>

- Install *DESeq2* on your personal computer if you have not done so.

## BEFORE WE START

We will go over an M and M example. 

Please count your M and Ms and put them here: https://tinyurl.com/MnMstats

## In-class Worksheets

1. We will start with counting reads using *featureCounts*. Since counting requires more resources, we will be coutning reads on the AWS server. 

- `Day7_featurecounts_worksheet.pdf` worksheet with details on running *featureCounts*
- `d7_featureCounts.R` R script with code to count reads using *featureCounts* from *Rsubread*
- `d7_featureCounts.sbatch` sbatch script that calls the aboive R scripts allowing us to count reads on the server

2. Second, we will perform differential gene expression analysis with *DESeq2*. This process will take in as input read counts from featureCounts.

- `Day7_differential_expression_worksheet.pdf` A worksheet with R code we will run in Rstudio to perform differential gene experession analysis

## Homework

The homework is to perform differential gene expession analysis on a new dataset. In the homework folder, the counts for the dataset are already generated. You will need to generate your own conditions table.

- `Day7_homework.pdf`
