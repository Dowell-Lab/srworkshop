## How do you test a tool or code to make sure its doing what you want?

Answer: Make simulated data. 


Go to this link and download the three files into a directory on your computer. 
https://drive.google.com/drive/u/2/folders/1Y30ygVpg-XF3iZAYIIAET23b9jb0vpcw

Then open `Day10/scripts/simuateT21.R` in Rstudio. 

I have already set it up to read in three files
1. `expression.csv`
2. `genotype.csv`
3. The gene gtf

Your job is to make simulated T21 individals by finding everyone with the genotype "D21" and every gene on chromosome 21 and muliplying the expression values for D21 invidials on chr21 genes by 1.5. 

Hints: people's psudonyms are in the column `Random_name` in both `genotype.csv` and `expression.csv`. All other columns in `expression.csv` are genes and the values in the `expression.csv` are the expression of that gene. 

Chromosome number is in the dataframe `gene_info` in a column called `seqname`.
