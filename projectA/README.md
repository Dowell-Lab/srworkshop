# Project A | Single-cell RNA-seq

## 0. We have gotten data from fetuses that are either Trisomy 21 or Eupliod. 

Jardine, L., Webb, S., Goh, I., Quiroga Londoño, M., Reynolds, G., Mather, M., ... & Haniffa, M. (2021). Blood and immune development in human fetal bone marrow and Down syndrome. Nature, 598(7880), 327-331.

Supplemental table 1
https://docs.google.com/spreadsheets/d/17JnPENn9CDxJa6igqSnC1O5RCjmZPxQ3bepf7jaS-MM/edit?gid=0#gid=0

Premade HD5 objects
https://developmental.cellatlas.io/fetal-bone-marrow

Link to tutorial for batch correction with Harmony: https://portals.broadinstitute.org/harmony/articles/quickstart.html. We don't go over it in these scripts due to time constraints.

## 1. Installing packages
If you haven't already, download the following packages locally onto R by following the steps below:

1. Open R (type R on the terminal and you should see something like the output below):
```
R version 4.4.0 (2024-04-24) -- "Puppy Cup"
Copyright (C) 2024 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin20

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.



```

2. Install Seurat according to instructions [here](https://satijalab.org/seurat/articles/install.html). TLDR version is below
*  Type this in R: `remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)`. You don't need to install the optional packages. 
* If it says somehting like remotes can't be found, do `install.packages('remotes')` and try again. 
* If it says that a dependency package isn't installing, try installing from source. For example, if I wanted to install the paackage "Rcpp" from source, I'd google cran Rcpp where the first link would get me to [here](https://cran.r-project.org/web/packages/Rcpp/index.html). I'd find a link to the latest release under the section Downloads, copy the link address of the appropriate version (in my case r-release(x86_64): Rcpp_1.0.12.tgz). Then in R I'd type `install.packages("link_I_copied", repo=NULL, type="source")` so in my case this would be `install.packages("https://cran.r-project.org/bin/macosx/big-sur-x86_64/contrib/4.4/Rcpp_1.0.12.tgz", repo=NULL, type="source")`
* If that doesn't work, install using `install.packages('Seurat').

2. Install CellChat according to instructions on the README [here](https://github.com/jinworks/CellChat)
* For troubleshooting dependencies, use the install from source tip from above
**Note** Make sure you install all of the specified dependencies before trying to install CellChat. 

## 2. Getting the data
If you were successfully able to install Seurat (don't worry if not Cellchat), download the Seurat objects from the zenodo link [here](https://zenodo.org/records/12725642).

If you were successfully able to install Seurat **and** CellChat, download the the CellChat objects from the zenodo link above. 

## 3. Open this cheat sheet.
   https://satijalab.org/seurat/articles/essential_commands.html
