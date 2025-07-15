# First, you will need to log into the AWS and follow:


1. Instructions for installing RSubread on AWS here `Day7_installing_Rsubread.pdf`

2. If you are online, you will be getting a simulated set of M&M candy. 

- First login to tha AWS and `cd` into `srworkshop` and run `git pull`. 
- Simulate getting candy by running the R code below. You can change the size of your handfull to anything between 50 and 150. Here I have it set to 100.
   `Rscript /Users/<yourusername>/srworkshop/projectB/day07/mnm_activity/grabahandful.R 100 greenbowl`
or
`Rscript /Users/<yourusername>/srworkshop/projectB/day07/mnm_activity/grabahandful.R 100 redbowl`

- Enter your candy counts to the following spreadsheet: https://tinyurl.com/MnMstats

3. Submitting read count script in `Day7_featurecounts_worksheet.pdf`

- To save on time and compute resources, you will first use a subsampled bam file to generate counts. 
 
- We are on the AWS because `Day7_featurecounts_worksheet.pdf` requires the bam files which are very BIG to load and process on your personal computer.

4. Then you will need to log off the AWS to do `Day7_differential_expression_worksheet.pdf`

- For this worksheet you will instead use the full counts data. The path is provided for you in the worksheet.

- It is easier to use Deseq2 on your own computer with Rstudio and the counts file you need for input is very small.

- At the end of today, you will be left with a list of differentially expressed genes. You will use this gene list later in the week.
