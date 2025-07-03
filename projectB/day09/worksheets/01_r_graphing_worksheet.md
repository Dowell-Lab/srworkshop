# Day9: Making Figures in R

Meaghan Courvan
24-07-18

## Learning objectives:

The goal of this section is to become familiar with working with dataframes in R, using R functions to answer biological questions, and making figures.

One common question people have when running multimodal experiments (for example with RNA-seq and ChIP-seq) is:  

* What genes are differentially expressed and have a p53 ChIP peak in their promoter?
* Do these genes fall into one functional classification, molecular pathway, or biological function?

We’ve gotten almost to answering this question throughout this project, and we’ll finish here. First, let’s recap what we’ve done.

1.  We calculated differential gene expression on Day 7.
2.  We called p53 ChIP peaks on Day 8.
3.  Earlier today, we used bedtools to get a list of p53 ChIP peaks overlapping genes in control (DSMO) conditions.
4.  We got an equivalent list of p53 ChIP peaks overlapping genes in experimental (Nutlin) conditions.

All that’s left is to do some basic set manipulation here in R and create some figures.

## Download the necessary data from the AWS
1. On your **local computer**, pick a location that's convenient for you and make a "day9" folder.
2. Inside your "day9" folder, make a "data" folder and a "results" folder.
3. From the AWS, you'll transfer your bedtools results: *</scratch/Users/YOUR_USERNAME/day9/bedtools_results>* to the "day9/data" folder on your local computer
4. From the AWS, transfer the DESeq2 results I am providing for you: */scratch/Shares/public/sread2025/cookingShow/day9b/deseq_output* to your "day9/data" folder.

## Now work in the *day9_graphing.R* script
### First, set up the environment.

At the top of your script, set your working directory and the results directory. 

### Read in the data

* Load in bed files you created that contain p53 peaks that overlap genes, in both DMSO and Nutlin treated conditions. In case you’ve forgotten, the function you need is `read.table()`. We just downloaded these above, and they should be in something like "*your/path*/day9/data/bedtools_results" folder.  
* Since we’re working with a bed file, it doesn’t have column names. In the `read.table()` function, make sure to set header=FALSE.
* Now that you’ve loaded the bed files, give them usable column names. I’ve provided those for you in the block of code below, you just need to assign them to the dfs.
* Load the differential expression data for HCT116 cells, which you downloaded above into something like "*your/path*/day9/data/deseq_output". Note that this file DOES have a header with the column names.
* There are a lot of genes in this data frame, so we’re going to filter down to a smaller number to create the venn diagram. Make another list that has only genes with an adjusted p-value \< 0.05 and an absolute value log<sub>2</sub>(Fold Change) \> 1.25.

``` r
### --- READ IN CHIP DATA 
  dmso_peaks # read in data here     
  nutlin_peaks # read in data here 
  
### --- GIVE CHIP DATA USEFUL COLUMN NAMES
  set_colnames <- c('chip_chr', 'chip_start', 'chip_end', 'peak_id', 'score', 
                    'chip_strand', 'signalValue', 'log_pval', 'log_qval', 
                    'summit', 'promoter_chr', 'promoter_start', 'promoter_end', 
                    'gene', '.', 'gene_strand', 'overlap') 
  colnames(dmso_peaks) <- set_colnames
  colnames(nutlin_peaks) <- set_colnames
  
### --- READ IN DE DATA 
  de.import # read in data here 
  
### --- EXTRACT SIGNIFICANTLY CHANGING GENES
```

<br>

### Examine lists for overlap

#### Working with lists of genes

These functions will do most of the work when we’re comparing lists of genes.

1.  **union()** ~ This function combines two lists and returns the
    unique elements.
2.  **intersection()** ~ This returns only elements that are common to
    both lists. <br> <br>

**First:** What is the set of genes bound by p53 in either DMSO or Nutlin treated samples? *Hint: we’re comparing lists, not dataframes.*

``` r
### --- Use either of the functions above to get the set of genes bound in either condition.
    # Hint: we’re comparing lists, not dataframes.

all_p53_bound <- 
```

**Second:** Look at the intersection of this list with all genes.

``` r
### --- Use either of the functions above to get the set of genes with p53 bound that are present at all in the DESeq2 results. 

bound_genes <- 
```

**Third:** Look at the intersection of this list with <u>differentially expressed</u> genes.

``` r
### --- Use either of the function above to get the set of genes with p53 bound that are differentially expressed. 

bound_de_genes <- 
```

Print out the number of genes that are differentially expressed and have a p53 ChIP peak in any condition. Also, print out the list of genes.

### Creating figures in R.

**Next** we’ll make a Venn diagram showing the overlap of genes with ChIP peaks from the Nutlin samples, genes with ChIP peaks from the DSMO samples, and genes wtih a ChIP peak that are differentially expressed.

To make figures in R, we can call a specific package or function that makes the type of graph we’re looking for. In this case, we will install a package to make venn diagrams, [ggvenn](https://github.com/yanlinlin82/ggvenn). 

* Visit that website and follow the installation instructions. 
* Load the library at the top of this next code block. 
* Look at the Quick Start section of the documentation. Do you see that the function expects a lists of lists? 
* Make a list of lists that contains genes with p53 peaks under DMSO treatment, genes with p53 peaks under Nutlin treatment, and genes with a p53 ChIP peak that are differentially expressed.

```r
#install on this line. Once you install, comment this line out 
library() # don't forget to load your new library on this line 
venn_list <- list() # add everything that you want to compare into this list of lists 
#create venn diagram on this line!
```
**Finally**, report the percentage of any ChIP peaks that overlap with a differentially expressed gene.

```r
# Just use R to do basic math!
# You can get the numbers you need from the venn diagram you made. 
```

### Export our work for GO analysis

We will do GO analysis on genes with a p53 ChIP peak that are also differentially expressed. For this, we want to have a list of our genes of interest and a list of background. **Background** means the set of genes we actually measured in this experiment. For us, this is the list of genes that are measured in our DESeq data that have any ChIP peak. Our **genes of interest** are those genes with a ChIP peak that are differentially expressed in the DESeq data. 

- Save each list in a place that is descriptive and convenient for you.

*How to use write.table()*:\
`write.table()` has a lot of options. These are the ones you need to set: 

* col.names=FALSE
* row.names=FALSE
* quote = FALSE

```r
write.table(  , file =  paste0(outdir, 'genes_p53peak.txt'),  )
write.table(  , file = paste0(outdir, 'genes_sig_p53peak.txt'),  )

```
> NOTE: the `paste0` command concatenates two strings.

### Comparing DE genes between cell lines.

This study looked at p53 responses across 4 cell lines: SJSA, MCF7, HCT116, and HCT116 with p53 KO. Let’s see how consistent the p53 response is across these cells lines. We’ll do this by making two types of graphs.

1.  A venn diagram, as we did above.
2.  A heatmap of differential expression.

**First**, load in the data from *"day9/data/deseq_output/"* folder.

``` r
# ---- READ IN DATA
de_hct116 <- #
de_hct116_p53ko <- #
de_mcf7 <- #
de_sjsa <- #
```

**Second**, make a venn digram like we did above with all analyzed genes. Is this a helpful visualization? What do you like about it? What are the limitations of the venn diagram? Can you change the input data to create a more helpful version of this same venn diagram?

``` r
# First venn 
venn_list_de <- list()
ggvenn()

# Do you think you can make a more helpful version of this?
# If so, put your code below 
```

**Third**, we’ll experiment with heatmaps as an alternative to venn diagrams. I like heatmaps for a couple reasons. They make it easier to visualize many samples, and also they show more of the quantitative data. We’ll use the package [pheatmap](https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/) (aka pretty heatmap).

You may have already installed this package during Day 7. If not, click on the link above, and skim the tutorial. The first lines tell you how to install pheatmap. Install it below. Then comment out the install line and load the package on the line below.

**What do we want the heatmap to display?**

- There are two quantitative variables you get out of DESeq2, the log<sub>2</sub>(Fold Change) and the adjusted p-value. Which do you want to plot?
- If we were to plot the data frames we imported, we’d be plotting tens of thousands of genes. That’s far too many to look at! Especially since the large majority of them aren’t significantly changing and therefore aren’t interesting. You need to pick a subset of genes to plot. What criteria do you use to select a smaller list of genes? How many genes do you want in your final list?

1.  Filter your lists of genes for those that are differentially expressed with the cutoffs of your choice. Print out the dimensions of each filtered data frame at the end. 
2.  If you’re happy with your filtering results, great! If you think you need more or fewer genes, adjust your filtering criteria.
3.  For the next step, you dataframes will need unique column names for values calculated by DESeq2, and they all need to have their gene column named the same. Right now they have the gene column named “GeneID”. We’ll leave that, but tag all the other column names with their cell type. I’ve provided that code, you just have to run it.

```r
### --- FILTER 
hct.filt <- 
hct_p53ko.filt <- 
sjsa.filt <- 
mcf7.filt <- 

### --- Print dimensions below


### --- Rename 

colnames(hct.filt)[-1] <- paste('hct', colnames(hct.filt)[-1], sep = '_')
colnames(hct_p53ko.filt)[-1] <- paste('hct_p53ko', colnames(hct_p53ko.filt)[-1], sep = '_')
colnames(sjsa.filt)[-1] <- paste('sjsa', colnames(sjsa.filt)[-1], sep = '_')
colnames(mcf7.filt)[-1] <- paste('mcf7', colnames(mcf7.filt)[-1], sep = '_')
```

These need to be combined into a common data frame that can be used to create a heatmap. There are several options for creating merged dataframes, the most popular being **full_join()**, **inner_join()**, **left_join()**, and **right_join()**. Today, we want to keep all of the entries in all of the dataframes. Choose which of the four functions to use.
> NOTE: We only really need the first three columns 
```r
  # Use one of the following functions 
  #     ~ full_join, inner_join, left_join, right_join   
  
  df <- # name the combined data frame just df

  # Move GeneID to be the row names of the df and remove the GeneID column
```
We want to graph genes that are differentially expressed in at least two out of three cell lines. Right now, if a gene was measured in one sample but not the others, when it was combined during full_join, NA was filled in for that gene. This is helpful for us because it’s easy to count the number of NAs per row using the **rowSums()** function and the **is.na()** function.

- Count the number of NA’s per row and save it as a new column titled “zero_count”.
- Filter the data frame so that you only keep genes with values in at least 2 cell lines.
- Print out the dimensions of this data frame.

There are two last things to adjust the data before we can make a heatmap.

* We’ve added an extra column that is unrelated to our data, the “zero_count” column. Remove this column.
* Pheatmap() will only take in numerical values, so change all your NA values to 0. *Hint* use the
**is.na()** function.

Lastly, graph your data using **pheatmap()**!

I am a picky and opinionated person. I dislike the default pheatmap esthetics. I think the middle of the range should be white, and that the highest values should be red and the lowest values be blue. I also don’t like that each box is outlined in gray, I find it distracting. 

In the Console, type “?pheatmap” and look at the documentation. 

* Find how to specify color scale and outline of the box.
* Rewrite the pheatmap line of code below making the esthetic changes.
* Add any other changes you’d like to see (i.e. text color, row labels, or font sizes).
* I've made plenty of asthetic demands, so I'll just leave a question/hint here for you to ignore if you want to. Is everything in this graph useful? Are there aspects of it that don't add to our understanding of the data?


## Optional Challenge Exercises (if you finish early)

**Question 1.** \
There are more ChIP peaks in the Nutlin samples than the DMSO. Are there other differences between the samples? Let's look at the width of the ChIP peak calls in both samples. This is an opportunity for you to practice doing some more dataframe manipulation and try new types of graphs. I will walk you through how to do compare the distributions of ChIP peak sizes.

These following 3 directions you need to do twice, once for nutlin_peaks, and once for dmso_peaks. 
* Make a new data frame that has only the columns you need for plotting. From the {sample}_peaks data frame, select the columns contain where the peak starts, peak stops, the peak ID, and the gene.
* Look at the dataframe. We can't use it as is. Why is that? Fix that problem using the function duplicated(), which identifies duplicated rows.
* Create a new column in this dataframe that contains the width of the ChIP peak. 

Now that you have these two new dataframes, we can create 1 dataframe and get it into a format that's good for plotting. 
* Combine nutlin_peaks_unique and dmso_peaks_unique into one dataframe using the rbind() function. This function pastes one dataframe on top of the other. It requires that the columns have the same names. This should already be the case (unless you changed anything).
* This dataframe is in LONG format. Read this webpage on the [difference between long and wide formats](https://www.statology.org/long-vs-wide-data/).
* Make sure to read the sections on <u>when</u> to use wide or long format.
* Since this is long format, we need a column that specifies which row is "Nutlin" and which is "DMSO". This is a bit of an odd command to write, so I've written it for you. You can just run it. If you want to understand it, ask me or ChatGPT to explain it to you.

Let's plot!!! \
Now you're ready to compare the distributions of ChIP widths. What sorts of graphs compare distributions? Well, there are [boxplots](https://ggplot2.tidyverse.org/reference/geom_boxplot.html), [violin plots](https://ggplot2.tidyverse.org/reference/geom_violin.html), and [ridgeplots](https://r-graph-gallery.com/294-basic-ridgeline-plot.html). I've listed them in order of complexity. The first two are built into ggplot. Ridgeplot requires installation of another package. Decide how much work you want to do and make one of these plots. 

**Question 2.** \
Have you heard of *GO analysis*? GO stands for gene ontology. I bet many of you have read this is papers! It's a way to see if your statistically significant results all belong to the same biological pathway, molecular process, or compartment of the cell. Move on to the next worksheet to learn more about this. 
