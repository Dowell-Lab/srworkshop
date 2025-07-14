# First, you will need to log into the AWS and follow:
 
1. Instructions for installing RSubread on AWS are here `Day7_installing_Rsubread.pdf`
2. Submitting read count script with `Day7_featurecounts_worksheet.pdf`

To save on time and compute resources, you will first use a subsampled bam file to generate counts. 
We are on the AWS because Day7_featurecounts_worksheet.pdf requires the bam files.

3. Then you will need to log off the AWS to do `Day7_differential_expression_worksheet.pdf`

- For this worksheet you will instead use the full counts data. The path is provided for you in the worksheet.

- It is easier to use Deseq2 on your own computer with Rstudio and the counts file you need for input is very small.

- At the end of today, you will be left with a list of differentially expressed genes. You will use this gene list later in the week.
