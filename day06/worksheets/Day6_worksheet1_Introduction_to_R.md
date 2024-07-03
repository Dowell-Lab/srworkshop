# Short Read Day 6 Worksheet | Introdution to R and RStudio (R you ready to learn R?!)

Authors: Taylor Jones (2022), Rutendo Sigauke (2024)

## Before class

- Make sure you have [RStudio](https://www.rstudio.com/products/rstudio/download/) and the latest version of R (R 4.4.1) installed on your computer!

- You want to get RStudio Desktop - Free

- R is awesome for plot generation and any statistical analysis you want to do- even those that aren’t related to next gen sequencing!

## Useful resources

- [Introduction to R](https://cran.r-project.org/doc/manuals/r-release/R-intro.pdf): This is an introduction to the ins and outs of R syntax.

- [R for Data Science (2e)](https://r4ds.hadley.nz/): This discusses applying R to data science with tips on exploratory data analysis.

- [The Pirate’s Guide to R](https://bookdown.org/ndphillips/YaRrr/): This is an excellent (and famous) comprehensive guide to R. I highly recommend this for learning as it starts at the basics and helps you build up skills for data analysis. It is also an easy read.

- [Beautiful Plotting Tips](http://zevross.com/blog/2014/08/04/beautiful-plotting-in-r-a-ggplot2-cheatsheet-3/): Cheat sheet for making beautiful, publication worthy plots in ggplot2!

## Useful packages

- Data isn’t always beautiful. To tidy data up there is a very useful package called tidyvrse with instructions in R for Data Science. Info can be found here: https://www.tidyverse.org/

- To plot data beautifully (ggplot2 is also a part of tidyverse): https://ggplot2.tidyverse.org/

## Section A: Getting familiar with R

You can access R through the command line or through RStudio. We will start on the command line. 

We use this on a big cluster when a viz instance isn’t available. You can write an R script and submit it with a sbatch script. 

1. Go ahead and launch R in the command line on the AWS server. Different versions of R are already installed. We can check which R is installed by running `module spider R`

Load the newest version of R listed using `module load R/3.6.1`

![Module load](images/module_load_R.png)

Lauch R in the terminal! Do this by typing R and Entering. Notice when you do this your username will be replaced by `>` indicating you are in R.

![Launch R](images/launch_r.png)

Now we can look at some basic R syntax.

2. Now that you are in R do a basic calculation of your choosing to test it out!

```
19 + 22
```

![Add in R](images/adding_in_R.png)

Here I typed 19+22(Enter) and it returned the answer to the summation which is 41.

3.  Set some variables. In R you can set variables that represent numbers, lists, dataframes etc. Let’s start
with numbers. To create a variable, you use the following syntax: `VARIABLE_NAME <- CONTENT` For example:

![Set variable](images/set_variables_R.png)

4. We can perform arithmetic on these variables now that they are set. We get the same answer, 41.

![Add variable](images/add_numeric_variables_R.png)

5. Now let’s set a slightly more complex variable. Let’s make a vector. To create a vector, we used the `c()` function for combine or concatenate. 

```
a_vector <- c(2, 10, 27, 100, 36)
```

![A vector](images/create_cvector.png)

You can see it returns the vector we just created.

Here I made a vector. Let’s call the vector to check its contents.

![call a vector](images/call_a_vector.png)

6. We can check data classes with `class()`. This can be applied to vectors, lists, data.frames etc.

```
class(a_vector)
```

6. Let’s test some basic functions on this list such as sum, mean, and standard deviation.

![operations on vector](images/operations_on_vector.png)

You invoke R functions by calling the function name (such as `sum`) followed by your variable in parentheses (such as `a_vector`). Some functions require more than one input, with the inputs separated by comma. Let’s exit command line R and switch to RStudio to try it out!

7. Exit R in the command line by typing `q()`, saving the workplace image saves the work we did.

![Save workspace](images/save_workspace.png)

## Section B: Exploring R in RStudio

Section B will be in RStudio. I have provided 2 different scripts that we will investigate. We will start with Learning_R.R in RStudio.
Open up RStudio on the viz instance. Once it is open go to Fil -> Open File -> Select Learning_R.R 

Further instructions are in this script!

## Section C: Submitting an R script as an sbatch job

The other script is Learning_R_submit_aws.R! This is an example R script where we submit R code as an sbatch job. Therefore, we save all outputs in order to view them after the process is complete. Submitting R scripts as an sbatch job is very useful when you have compute intensive processes (examples are counting reads).
