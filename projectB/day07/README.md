# Project B Day 7 | RNA-seq: Counting Reads and Differential Expression

Today we will count RNA-seq reads and use the counts to run *DESeq2*.

First, you will learn to count reads using *featureCounts*. Then, you will learn about exploratory data analysis and QC checks associated with differential analysis. Finally, you will use the differential analysis pipeline *DESeq2* to identify a list of genes which are differentially expressed in our dataset. This list will be used later in the week when we integrate it together with our ChIP-seq results.


## BEFORE WE START

We will do the M and M example. 

Please count your M and Ms and put them here: https://tinyurl.com/MnMstats

- If you are online, I'm sorry we can't give you candy. 

i. First login to tha AWS and `cd` into `srworkshop` and run `git pull`. 
ii. Simulate getting candy by running the following R script.

- You can change the size of your handfull to anything between 50 and 150. Here I have it set to 100.

To sample from the _Green Bowl_, run:
```Rscript /Users/<yourusername>/srworkshop/projectB/day07/mnm_activity/grabahandful.R 100 greenbowl```

To sample from the _Red Bowl_, run:
```Rscript /Users/<yourusername>/srworkshop/projectB/day07/mnm_activity/grabahandful.R 100 redbowl```


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
