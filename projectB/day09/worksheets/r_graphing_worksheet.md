Day9: Making Figures in R
================
Meaghan Courvan
24-07-18

### Learning objectives:

The goal of this section is to become familiar with working with dataframes in R, using R functions to answer biological questions, and making figures.

One common question people have when running multimodal experiments, for example with RNA-seq and ChIP-seq is:  

* What genes are differentially expressed and have a p53 ChIP peak in their promoter? 
* Do these fall into one functional classification, molecular pathway, or biological funciton.

We’ve gotten almost to answering this question throughout this project, and we’ll finish here. First, let’s recap what we’ve done.

1.  We calculated differential gene expression on Day 7.
2.  We called p53 ChIP peaks on Day 8.
3.  Earlier today, we used bedtools to get a list of p53 ChIP peaks
    overlapping genes in control (DSMO) conditions.
4.  We got an equivalent list of p53 ChIP peaks overlapping genes in
    experimental (Nutlin) conditions.

All that’s left is to do some basic set manipulation here in R and create some figures.

### First, set up the environment.

At the top of your script, load any R libraries you want to use. 

### Read in the data

We will be working with bed files we generated earlier today, as well as differential expression data generated on Day 7. We’ll load that all in here.

* Load in bed files you created that contain p53 peaks that overlap genes, in both DMSO and Nutlin treated conditions. In case you’ve forgotten, the function you need is “read_tsv()”. 
* Since we’re working with a bed file, it doesn’t have column names. In the “read_tsv()” function, make sure to set col_names=FALSE.
* Now that you’ve loaded the bed files, give them usable column names. I’ve provided those for you in the block of code below, you just need to assign them to the dfs.
* Load the differential expression data you generated on Day 7.
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

<br>

#### Working with lists of genes

These functions will do most of the work when we’re comparing lists of genes.

1.  **union()** ~ This function combines two lists and returns the
    unique elements.
2.  **intersection()** ~ This returns only elements that are common to
    both lists. <br> <br>

**First:** What is the set of genes bound by p53 in either DMSO or Nutlin treated samples? *Hint: we’re comparing lists, not dataframes.*

``` r
### --- Use either of the function above to get the set of genes.
    # Hint: we’re comparing lists, not dataframes.

all_p53_bound <- 
```

**Second:** Look at the intersection of this list with DE genes

``` r
### --- Use either of the function above to genes with p53 bound that are differentially expressed. 

all_p53_bound <- 
```

Print out the number of genes that are differentially expressed and have a p53 ChIP peak in any condition. Also, print out the list of genes.

### Creating figures in R.

**Next** we’ll make a Venn diagram showing the overlap of these categories.

To make figures in R, we usually call a specific package or function that makes the type of graph we’re looking for. In this case, we will install a package to make venn diagrams, [ggvenn](https://github.com/yanlinlin82/ggvenn). 

* Visit that website and follow the installation instructions. 
* Load the library at the top of this next code block. 
* Look at the Quick Start section. Do you see that the function expects a lists of lists? 
* Make a list of lists that contains genes with p53 peaks under DMSO treatment, genes wtih p53 peaks under Nutlin treatment, and differentially expressed genes.

```r

#install on this line. Once you install, comment this line out 
library() # don't forget to load your new library on this line 
venn_list <- list() # add everything that you want to compare into this list of lists 
#create venn diagram on this line!

```
**Finally**, report the percentage of ChIP peaks that overlap with a differentially expressed gene.

```r
# Just use R to do basic math!
# You can get the numbers you need from the venn diagram you made. 
```

### Export our work for GO analysis

For GO analysis, we want to have a list of our genes of interest and a list of background. **Background** means the set of genes we actually measured in this experiments. We’ll create the background list, and then export *(1)* this background list, along with *(2)* all differential expressed and bound genes, *(3)* overexpresed and bound genes, *(4)* underexpressed and bound genes.

- Make each list
- Now you need to turn each list into a dataframe. This is because the
  function used to write, write_tsv(), only write dataframes and not
  lists. Use the function as.data.frame() to convert each list.
- Now save each list in a place that is descriptive and convenient for
  you.

### Comparing DE genes between cell lines.

This study looked at p53 responses across 4 cell lines: SJSA, MCF7,
HCT116, and HCT116 with p53 KO. Let’s see how consistent the p53
response is across these cells lines. We’ll do this by making two types
of graphs.

1.  A venn diagram, as we did above.
2.  A heatmap of differential expression.

**First**, load in the data from \_\_\_\_\_\_\_\_\_\_\_.

**Second**, make a venn digram like we did above. Is this a helpful
visualization? What do you like about it? What are the limitations of
the venn diagram? Can you change the input data to create a more helpful
version of this same venn diagram?

``` r
# more helpful version of the venn diagram 
```

**Third**, we’ll experiment with heatmaps as an alternative to venn
diagrams. I like heatmaps for a couple reasons. They make it easier to
visualize many samples, and also they show more of the quantitative
data. We’ll use the package
[pheatmap](https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/)
(aka pretty heatmap).

Click on the link above, and skim the tutorial. The first lines tell you
how to install pheatmap. Install it below. Then comment out the install
line and load the package on the line below.

**What do we want the heatmap to display?**

- There are two quantitative variables you get out of DESeq2, the
  log<sub>2</sub>(Fold Change) and the adjusted p-value. Which do you
  want to plot?
- If we were to plot the data frames we imported, we’d be plotting tens
  of thousands of genes. That’s far too many to look at! Especially
  since large majority of them aren’t significantly changing and
  therefore aren’t interesting. You need to pick a subset of genes to
  plot. What criteria do you use to select a smaller list of genes? How
  many genes do you want in your final list?

1.  Filter your lists of differentially expressed genes below. Print out
    the dimensions of each filtered data frame at the end. Are all the
    data frames the same size? Does anything surprise you?
2.  If you’re happy with your filtering results, great! If you think you
    need more or fewer genes, adjust your filtering criteria.
3.  For the next step, you df’s will need unique column names for values
    calculated by DESeq2, and they all need to have their gene column
    named the same. Right now they have the gene column named “GeneID”.
    We’ll leave that, but tag all the other column names with their cell
    type. I’ve provided that code, you just have to run it.

Now I need to teach you how to combine these into a common data frame
that can be used to create a heatmap. We’ll be doing this today using
the functions **full_join()**, **left_join()**, and **right_join()**.

We want to graph genes that are differentially expressed in at least two
out of three cell lines. We’ll filter the data frame for these genes,
then replace NAs with 0.

Right now, if a gene was measured in one sample but not the others, when
it was combined during full_join, NA was filled in for that gene. This
helpful for us because it’s easy to count the number of NAs per row
using the **rowSums()** function and the **is.na()** function.

- Count the number of NA’s per row and save it as a new column titled
  “zero_count”.
- Filter the data frame so that you only keep genes with values in at
  least 2 cell lines.
- Print out the dimensions of this data frame.

There are two last things to adjust the data before we can make a
heatmap. \* We’ve added an extra column that is unrelated to our data,
the “zero_count” column. Remove this column \* Pheatmap() will only take
in numerical values, so change all your NA values to 0. *Hint* use the
**is.na()** function.

Lastly, graph your data using **pheatmap()**!

I dislike the default pheatmap esthetics. I think a log<sub>2</sub>(Fold
Change) = 0 should be white. I also don’t like that each box is outlined
in gray, I find it distracting. In the Console, type “?pheatmap” and
look at the documentation. Find how to specify color scale and outline
of the box. Rewrite the pheatmap line of code below making the esthetic
changes. Add any other changes you’d like to see. For example, text
color, row labels, or font sizes.
