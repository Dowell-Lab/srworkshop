# Day 4 Worksheet – Trimmomatic 
Author: Jessica Westfall, 2021\
Edited: Lynn Sanford, 2024

## Introduction
Now that we have evaluated our sequence library initially to determine if the libraries are worth analyzing, we will do some “cleaning up” by trimming unwanted sequences such as adapter sequences. This step is necessary for improved alignment and mapping to the reference genome downstream. Once trimming is completed, we will reevaluate our trimmed files with FastQC for quality to decide if we will move forward with mapping.
>Note: The directory and username used in the screenshot will be for my working directory and username and will be different than yours.

## Make working directories and copy script
Yesterday, we made working directories for running fastQC. Repeat the same process, but this time we will make a directory for trimmomatic.

- Use command `pwd` to determine what directory you are in and if necessary, `cd` to the directory that you want to place your new trimmomatic directory in. 
- Make several new directories and confirm the folders are present.

  ![Make trimming directories](md_images/make_trimming_dirs.png)

- Navigate to the srworkshop repo and `git pull`. Then `rsync` the `d4_trim_qc.sbatch` script from `day04/scripts/` into your script directory. Confirm the file is present in the directory with `ls`.
  - Remember, to copy the script, the command syntax is `rsync <input> <output>`

  ![Rsync trimming script](md_images/rsync_trimming_script.png)

## Trimmomatic

Edit the sbatch script by using `vim <scriptname>` to open your sbatch script in Vim.

Similar to the previous exercise you will need to change the job name, user email, and the standard output and error log directories.
- Change the `job-name=<JOB_NAME>` to a name related to the job you will be running, for example `trim_qc`. Additionally you will want to change the `–mail user=<YOUR_EMAIL>` to your email, as well as the path to your output/error directory for the standard output (`--output`) and error log (`--error`). If you include the terms `%x` and `%j` in your filepaths, `%x` will be replace by your job-name and the `%j` will be replace by the job id that will be assigned by Slurm when you run your sbatch script.

    ![Trimming sbatch header](md_images/trimming_sbatch_header.png)

- Note also that we have changed the numbers of processors (`--ntasks`) for the job, as Trimmomatic can use multiple processors per input file. We'll request 1 node, 8 ntasks, 8gb of memory and 1 hr of wall time. 

Assigning path variables will make your scripts easier to read. In addition, this makes it easier to reference a given path and utilize it in your scripts.
- Here we define six variables.
- For the `INDIR=` (input directory), change the path to the directory where the fastq data files are located.
- For the `OUTDIR=` (output directory), point to the appropriate output file directories for our fastQC and trimmed fastq files.
- Note we make sure that the output directories exist by calling `mkdir` within the sbatch script. Look up what the `-p` parameter does.

  ![Trimmomatic variables](md_images/trimmomatic_variables.png)
 

Load the require modules for running this pipeline. We will be using FastQC and the trimming program trimmomatic. If you are not sure which version of a program is available on the cluster you can save and quit out of Vim, then use the command `module spider <program>` to find the available versions.
- Note again, though, that the version of fastqc that is actually available to us is _____.

  ![Trimmomatic modules](md_images/trimmomatic_modules.png)
 
For the meat of the script, we will be running 3 analysis steps. Several analysis steps run together is called a pipeline. Our pipeline allows us to:
1. Run FastQC on raw sample fastq files
2. Trim the fastq files to remove adapters and other sequence elements
3. Reevaluate the quality of the trimmed fastq files with FastQC.

  ![Trimmomatic pipeline](md_images/trimmomatic_pipeline.png)

The FastQC steps should look familiar to you. Now let's look closer at the trimming step.
- In this script we are running paired end reads. Trimmomatic can be used on both single end or paired-end reads. When setting your parameters use the appropriate adapters. 
- Below are the syntaxes needed to run trimmomatic:
  - For single-end reads\
    ```
    java jar /opt/trimmomatic/0.36/trimmomatic-0.36.jar SE [ -threads <n> ]  [ -phred33 | -phred64 ] [ -trimlog <output_trimlog> ] <input_file> <output_file> ILLUMINACLIP
    ```
  - For pair-end reads\
    ```
    java jar /opt/trimmomatic/0.36/trimmomatic-0.36.jar PE [ -threads <n> ]  [ -phred33 | -phred64 ] [ -trimlog <output_trimlog> ] <input_file1> <input_file2> <output_fileP1> <output_fileU1> <output_fileP2> <output_fileU2> ILLUMINACLIP
    ```
  - In the script, you'll see this command is split into multiple lines. Recall that the `\` at the end is used to break the code up for ease of reading. If `\` does not change color as you see above (it may not be the same color as here, but it should change from the normal text color), you may have an extra space after the `\`. Remove that space or your code will not run properly.

- The final parameter shown in the Trimmomatic commands above is the `ILLUMINACLIP` parameter, which should be replaced by a series of strings that define the trimming behavior you need. Below are several possible trimming strategies that can be implemented, singly or in combination:
  - `ILLUMINACLIP:<path_adapters_fasta>:<seed_mismatches>:<palindrome_clip_threshold>:<simple_clip_threshold>`
  - `LEADING:<quality>`
  - `TRAILING:<quality>`
  - `SLIDINGWINDOW:<window_size>:<required_quality>`
  - `MINLEN:<length>`

Each of the trimming behaviors is described in much more detail in the Trimmomatic manual. Find the manual online and read more about the two we're using in this script, `ILLUMINACLIP` and `CROP`.

After saving all of your changes to the script, run it!

- Submit the job to the job manager SLURM (`sbatch <scriptname>`). The job manager will assign a job id to your run.
  - This pipeline has more tasks than the previous worksheet, so you will want to check the status of your job (`squeue -u <username>`) to see if the job is running (`R`) or completed (`C`).
  - If there are any errors (often times these are just typos in your scripts), you will want to access your error file to make necessary corrections. This is the path/filename specified in the header of the sbatch script. You can view this file using `more`, `less`, `cat`, or open it in Vim.

Check that you have your expected outputs and that they are not empty.
- `ls -lh` displays the long-form list, which includes the size of your files. The sizes should be larger than 0.