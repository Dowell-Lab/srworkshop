# Day 4 Homework 
Authors: Jessica Westfall & Rutendo Sigauke\
Edited: Lynn Sanford, 2026

## Introduction
In a future day of the workshop we will go into more details about RNA-seq libraries. This homework will go over the tasks that we did in class and provide more practice and questions to expand your understanding.

## Homework

### 1. Open the homework script
  Copy the `example_process_rnaseq.sbatch` script from the GitHub repo `day04/scripts` folder to your scratch day4 working directory and make the necessary edits to do the following tasks. 
### 2. FastQC
Evaluate the FASTQ files in `/scratch/Shares/public/sread/homework_data_files/day4/`.

How is the quality of these sequence libraries? Things we want to look at are:
- GC content (Is the library contaminated?)
- Adaptor content (Did the sequencer read into our adaptors?)
- Read duplication (Is our sample overamplified? Depends on library type…)
- Sequence quality/N content (How confident was the sequencer in calling each base?)
- Sequence quality based on flow cell location (Was there a sequencing failure?)
- Base identity at each location (Was there any bias in amplification/ligation?)

### 3. Trimmomatic
Trim the sequence library to remove adapters. Save the output to have the suffix `_trim.fastq` to track the trimming.

When trimming, consider if you are trimming single-ended (SE) or paired-ended (PE) reads. Consider other parameters listed for trimming:
- `ILLUMINACLIP`: Cuts adapter and other illumina-specific sequences from the read
- `SLIDINGWINDOW`: Performs a sliding window trimming approach. It starts scanning at the 5' end and clips the read once the average quality within the window falls below a threshold
- `LEADING`: Cuts bases off the start of a read, if below a threshold quality
- `TRAILING`: Cuts bases off the end of a read, if below a threshold quality
- `CROP`: Cuts the read to a specified length by removing bases from the end
- `HEADCROP`: Cuts the specified number of bases from the start of the read
- `MINLEN`: Drops the read if it is below a specified length

### 4. FastQC of trimmed files
Run FastQC on the trimmed FASTQs to reevaluate the quality.  

Take a look at the FastQC report for the trimmed fastq and ask yourself:
- Are these files trimmed well?
- Can you adjust the parameters to make the trimming more stringent to remove adapter content? 

### 5. HISAT2
Edit the sbatch script to map the two corresponding paired-end FASTQ files (`end1` and `end2`).

Mapping has different parameters to change the mapping efficiency. What would happen if you altered the script which currently has `--very-fast` to `--very-sensitive`? (These parameters are best explained in the Bowtie2 documentation. You also may need to increase your wall time to test them)

### 6. SAMtools
Now we have a huge SAM file.  

Convert the SAM file to a BAM file to produce a compressed binary file that takes up less space.
- What is the size difference between SAM versus BAM?
- What is the difference between the two filetypes that contributes to the size difference?

Sort and index your BAM files.

Remove the unnecessary/redundant files to save space.

Back up scripts and final output files to your home directory.

### 7. IGV
Transfer your `.bam`/`.bam.bai` files to your local computer.

Open IGV Desktop OR the IGV Web App and visualize the sorted BAM files.
- Do you know which reference genome your sequence reads were aligned and mapped to?
- Since both hg19 and hg38 are different versions of the human genome, are they interchangeable in IGV?
- How can you tell if you have the right genome?