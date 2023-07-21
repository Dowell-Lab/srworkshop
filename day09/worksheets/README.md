# Short Read Day 9: Peak-centric sequencing worksheets
- Author: Jessica Huynh-Westfall

## Part 1: QC of your ChIP-seq data

As we discussed last week, ensuring your data's quality is a pertinent step before you move forward with your data analysis.

1. After you Log into the AWS, you will make a directory for day 9. Make the subdirectories including one for qc and its subdirectories. Create a working directory for day9 in scratch.
```
[<YourUsername>@ip-172-31-29-36 ~]$ mkdir day9
[<YourUsername>@ip-172-31-29-36 ~]$ cd day9
[<YourUsername>@ip-172-31-29-36 ~]$ mkdir scripts qc fastq bam sam
[<YourUsername>@ip-172-31-29-36 ~]$ cd qc
[<YourUsername>@ip-172-31-29-36 ~]$ mkdir fastqc hisat_mapstats multiqc
```

2. Copy scripts from github into your scripts directory
   a. d9_fastqc.sbatch
   b. d9_hisat2.sbatch
   b. d9_multiqc.sbatch

3. Copy the fastq qc files over from the scratch directory to your qc directory
```
[<YourUsername>@ip-172-31-29-36 ~]$ cp /scratch/Shares/public/sread2023/data_files/day9/fastq qc/fastq
```

4. 
