# Short Read Day 9: Peak-centric sequencing worksheets
- Author: Jessica Huynh-Westfall

## Part 1: QC of your ChIP-seq data

As we discussed last week, ensuring your data's quality is a pertinent step before you move forward with your data analysis.

1. After you Log into the AWS, you will make a directory for day 9 in your scratch directory. Make the subdirectories including one for day09, stdout and stderr, scripts, and fastq (which will be our data we will be using). Create a working directory for day9 in scratch.
```
[<username>@<hostname> ~]$ cd /scratch/Users/<YourUsername>
[<username>@<hostname> ~]$ mkdir day9
[<username>@<hostname> ~]$ cd day9
[<username>@<hostname> ~]$ mkdir scripts eofiles fastq 
```

2. Copy scripts from github repo into your scripts directory
   a. d9_preprocess_chipseq.sbatch
   b. d9_qc_chipseq.sbatch
   c. d9_macs.sbatch
   d. d9_bedtools.sbatch

3. Copy the fastq files over from the scratch directory to your fastq directory. We will only be using BACH* in class. The GAPBA* is another dataset you can practice with in the homework.
```
[<username>@<hostname> ~]$ cp <path/to/sr2023 git repo fastq files> /scratch/Users/<YourUsername>/day9/.
```

4. cd into your scripts directory. Edit and run scripts d9_preprocess_chipseq.sbatch and d9_qc_chipseq.sbatch. The preprocessing will run fastQC and HISAT2 alignment. d9_qc_chipseq.sbatch will run some programs to get QC report and multiQC will compile all the QC files for you to review the quality of the data.
- Preseq predicting and estimating the complexity of a genomic sequencing library. 
- RSeQC has multiple module for sequencing QC with read distribution being one of the programs.
- MultiQC will summarize the output from numerous bioinformatic tools and compile a report for you to review the QC logs.

## Part 2: MACS to call peaks
To study DNA enrichment assays such as ChIP-seq and ATAC-seq, we are introducing the analysis method, Model-based Analysis of ChIP-Seq (MACS). This method enables us to identify transcription factor binding sites and significant DNA read coverage through a combination of gene orientation and sequencing tag position.

1. Edit and run the d9_macs.sbatch. This script has multiple section that is detail out in the worksheet Day9_macs_worksheet.pdf. In addition to calling peaks with and without controls, the script also includes intersecting with blacklist regions. 

2. Pull the output peak files from the run into IGV to see the difference between all the different called regions from the run. For more information on what you should consider for parameters see the worksheet Day9_macs_worksheet.pdf. Use the homework for some thoughts to consider when looking at peak-centric data.

## Part 3: Bedtools
Annotation files are critical for downstream sequencing analysis, and text file manipulation and Bedtools are both very useful ways to work with them.

Follow the Day9_Bedtools_worksheet.docx now.
