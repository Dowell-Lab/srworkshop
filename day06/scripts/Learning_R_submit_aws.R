## Learning R : example for submitting in sbatch
## Author: Taylor Jones (2022), Rutendo Sigauke (2023)

###############################################################################
## In this script we plan to submit the R script as an SBATCH job            ##
## This is useful for when you have a computationally intensive job that     ##
## requires more resources than on your personal computer.                   ##
## For example, counting reads in a bam file over gene features for RNA-seq. ##
###############################################################################

##############################################################################
# Set your working directory first
workdir <- '/path/to/working/directory/'
setwd(workdir)
getwd()

##############################################################################
## In the Learning_R.R script you processed data(iris) and data(mtcars)     ##
## NOW, we are processing data(mtcars) data in a script called              ##
## Learning_R_inclass.R                                                     ##
##############################################################################

## NOTE: Submitting R job as a script is NOT interactive
## Therefore, you have to save all outputs in order to view them

# 1: First, load mtcars data
data(mtcars)

# Question: Where do you view the outputs from the following commands?
head(mtcars)
dim(mtcars)
summary(mtcars)

# 2: Save the data as a csv file
write.csv(mtcars,
          paste0(workdir,"mtcars.csv"))

# 3: Create a simple scatter plot of mpg vs. wt using base R plot() function
png(paste0(workdir,'mtcars_mpg_wt_scatterplot.png'))
plot(mtcars$mpg, mtcars$wt,
     col='steelblue',
     main='Scatterplot MPG vs WT',
     xlab='mpg',
     ylab='wt',
     pch=19)
dev.off()

# 4: Optional- create and save another pretty plot for mtcars data below

