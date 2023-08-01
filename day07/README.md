# Day 7 | RNA-seq: Counting Reads and Differential Expression

Today we will count RNA-seq reads and used the counts to run DEseq2.

DESeq2 the most commonly used program for differential expression assessment. We'll talk about how to interpret results and build quality designs.


###

BEFORE WE START
We will do the M and M example. 
Please count your M and Ms and put them here. 

https://docs.google.com/spreadsheets/d/1xaJAa6lmNU_PdwQs_-4M2thAbBNSiyuwXmRyp6NpP2w/edit#gid=1181236382

## Inclass Worksheets

1. We will start with counting reads using *featureCounts*. Since counting requires more resources, we will be coutning reads on the AWS server. 

- `Day7_featurecounts_worksheet.docx` worksheet with details on running *featureCounts*
- `d7_featureCounts.R` R script with code to count reads using *featureCounts* from *Rsubread*
- `d7_featureCounts.sbatch` sbatch script that calls the aboive R scripts allowing us to count reads on the server

2. Second, we will perform differential gene expression analysis with *DESeq2*. This process will take in as input read counts from featureCounts.

- `Day7_differential_expression_worksheet.docx` A worksheet with R code we will run in Rstudio to perform differential gene experession analysis

## Homework

The homework is to perform differential gene expession analysis on a new dataset. In the homework folder, the counts for the dataset are already generated. You will need to generate your own conditions table.

- `Day7_homework.pdf`
