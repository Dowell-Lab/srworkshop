# Downloading raw data from NCBI

## Summary about the data

The data downloaded in this example is from Andrysik et al. 2017 [doi:10.1101/gr.220533.117](https://genome.cshlp.org/content/27/10/1645). The study looks at the function of TP53 in three cancer cell lines. Authors use different genomic datasets to understand the function of the transcription factor. 

The data was deposited at the NCBI Gene Expression Omnibus (GEO) database [GSE86222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86222). This GEO accession number is usually reported in the research article.

## Getting the data

Raw sequenced data for the paper can be downloaded from NCBI's [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP083188&o=acc_s%3Aa). Since the raw samples are also hosted on the AWS, we can use `wget <url_to_raw_file>` to download the files. 

The following script `00_download_Andrysik2017ChIPseq.sbatch` uses `wget` to download the *sra* files, and `run_download.sh` loops through the paths (below) listed in the file `andrysik2017_chip.srr` and downloads them locally.

```
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR4090097/SRR4090097
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR4090096/SRR4090096
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR4090095/SRR4090095
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR4090094/SRR4090094
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR4090093/SRR4090093
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR4090092/SRR4090092
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR4090091/SRR4090091
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR4090090/SRR4090090
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR4090089/SRR4090089
```

The downloaded files are in *sra* format and we need to use the [`fastq-dump`](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) tool from the [sra-tools](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit) to extract the raw `fastq` files. 

> **Note:** There is also `fasterq-dump`, **BUT** it is memory hungry. If you do decide to use the `fasterq-dump` instead of `fastq-dump`, request more/adequate memory in your script.