# Day 7 Worksheet – Installing Rsubread		

Author: Rutendo Sigauke (2025)

This short worksheet goes over installing Rsubread in R on the AWS. This process can take a while and should be completed before class. RSubread can be found on Bioconductor and instructions to install the package can be found here: [https://bioconductor.org/packages/release/bioc/html/Rsubread.html](https://bioconductor.org/packages/release/bioc/html/Rsubread.html). 

## Install Rsubread from BiocManager

Installation can be done in the R console (shown below). 

1. First load a version of R installed on the AWS

```sh
module load R/4.3.1
```

2. Rsubread can be installed from the **BiocManager**

- Type `R` in the terminal 

> Note: To install Rsubread from BiocManager requires R 4.5.0. 
> Since were are working with R 4.3.1, the command below will not work.

- We have to first install **BiocManager**. **BiocManager** library can be installed from the R Comprehensive R Archive Network (CRAN).

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

```

- The above command will initiate a new personal library for your R packages. Select **yes**. Additionally, CRAN mirrors from where to download the packages. We will use the USA (KS) mirror.

- Now, we can install **Rsubread** from **BiocManager** to our **R** libraries.

```R
BiocManager::install("Rsubread")
```

> Note: This will take a few seconds. If the library is installed successfully, it can be loaded as shown below without any errors.

```R
library("Rsubread")
```

3. If the above does not work, try installing from source:**

In some cases, you may have to download an older version of Rsubread from source ([https://bioconductor.org/packages/3.17/bioc/src/contrib/Archive/Rsubread/](https://bio\
conductor.org/packages/3.17/bioc/src/contrib/Archive/Rsubread/)).

1. Make sure you have loaded R

```sh
module load R/4.3.1
```

2.  Second, download the tar.gz file with `wget` to a location on the AWS.

```sh
wget https://bioconductor.org/packages/3.17/bioc/src/contrib/Archive/Rsubread/Rsubread_2.14.0.tar.gz
```

3. Then open R, in an R console, install the package from the source file as shown below. 

- Type `R` in the terminal.

- Enter the conda below pointing to the path to the `Rsubread_2.14.0.tar.gz` file.

```R
install.packages(“/path_to_file/Rsubread_2.14.0.tar.gz”, repos = NULL, type="source")
```

