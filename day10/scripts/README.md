# Downloading raw data from NCBI

## Summary about the data

The data downloaded in this example is from Andrysik et al. 2017 [doi:10.1101/gr.220533.117](https://genome.cshlp.org/content/27/10/1645). The study looks at the function of TP53 in three cancer cell lines. Authors use different genomic datasets to understand the function of the transcription factor. 

The data was deposited at the NCBI Gene Expression Omnibus (GEO) database [GSE86222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86222). This GEO accession number is usually reported in the research article.

## Getting the data

Raw sequenced data for the paper can be downloaded from NCBI's [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP083188&o=acc_s%3Aa). Since the raw samples are also hosted on the AWS, we can use `wget <url_to_raw_file>` to download the files. 

The following script `00_download_Andrysik2017ChIPseq.sbatch` uses `wget` to download the *sra* files, and `run_download.sh` loops through the paths (below) listed in the file `andrysik2017_chip.srr` and downloads them locally.

```
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR4090097/SRR4090097
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR4090089/SRR4090089
```
**Note:** these are just a subset of the samples in this project.

1. Update the from the GitHub repository

```
cd /Users/your_username/srworkshop

git pull
```

2. Now go to your `scratch` folder and create a `day10` folder, and in that folder create subdirectories `e_and_o`, `scripts` and `sra`

```
cd /scratch/Users/your_username

mkdir day10

cd day10

mkdir e_and_o scripts sra

```

3. Go into the scripts folder and copy the following files ( `00_download_Andrysik2017ChIPseq.sbatch`, `run_download.sh` and `andrysik2017_chip.srr`) to your scratch scripts folder.

4. Edit the `00_download_Andrysik2017ChIPseq.sbatch` and  `run_download.sh` scripts so they point to your `/scratch/Users/your_username/day10` files.

5. After you have edited the files, you will run the `run_download.sh` script.

```
./run_download.sh 
```

This script will submit the `wget` output to you working directory. Once the script is done, you will need to move the files to a dedicated **sra** folder.

```
mv SRR* /scratch/Users/your_username/day10/sra/
```

## Converting the SRA files to FASTQ files

The downloaded files are in *sra* format and we need to use the [`fastq-dump`](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) tool from the [sra-tools](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit) to extract the raw `fastq` files. 

> **Note:** There is also `fasterq-dump`, **BUT** it is memory hungry. If you do decide to use the `fasterq-dump` instead of `fastq-dump`, request more/adequate memory in your script.

To download the *fastq* files we can now use the command below:

```
fastq-dump ${SRR} --gzip
```

or download in parallel using a wrapper for `fastq-dump` called [`parallel-fastq-dump`](https://github.com/rvalieris/parallel-fastq-dump):

```
parallel-fastq-dump --threads <NumCores> --gzip --sra-id ${SRR}
```