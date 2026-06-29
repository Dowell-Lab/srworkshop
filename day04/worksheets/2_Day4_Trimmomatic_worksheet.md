# Day 4 Worksheet – Trimmomatic 
Author: Jessica Westfall, 2021\
Edited: Lynn Sanford (2024), Georgia Barone (2025)

## Introduction
The first step in any short read analysis is to evaluate a sequence library to determine if it is worth analyzing, as you did yesterday. Often, the second step is to “clean up” the data by trimming unwanted sequences such as adapter sequences. This step is necessary for improved alignment and mapping to the reference genome downstream. Once trimming is completed, we will reevaluate our trimmed files with FastQC for quality to decide if we will move forward with mapping.

## Make working directories and copy script
- In your **home** directory, make a `workshop-day4` directory. Within that, make `scripts`, `fastqc`, and `hisat2` subdirectories.
- Navigate to your **scratch** directory and create the same directory with `scripts`, `eofiles`, `fastqc`, and `trimmomatic` subdirectories.

  ![Make trimming directories](md_images/make_trimming_dirs.png)

- Navigate to the srworkshop repo and `git pull`. Then `rsync` the `d4_trim_qc.sbatch` script from `day04/scripts/` into your **scratch** day4 scripts directory.

  ![Rsync trimming script](md_images/rsync_trimming_script.png)

## Trimmomatic

In your **scratch** day4 scripts directory, edit the trimming script in Vim.

You will need to change the job name, user email, and the standard output and error log directories. Use a useful name for the job name, such as `trim_qc`.
- Note also that we have changed the numbers of processors (`--ntasks`) for the job, as Trimmomatic can use multiple processors per input file. We'll request 1 node, 2 ntasks, 1gb of memory and 10 minutes of wall time.
- You'll also see just under the header is a section that specifies certain information about the job that might help in troubleshooting or in documenting your work. Since we're just printing this to standard output, this information will be stored in the `.out` file associated with the job.

Assigning path variables will make your scripts easier to read. In addition, this makes it easier to reference a given path and utilize it in your scripts.
- Modify the `outdir=` (output directory) to point to the appropriate output file directories for our fastQC and trimmed fastq files
- The rest of the variables have been filled in for you. Make sure you understand what each variable points to (if you don't, look for where they are called in the script and explore the fastq input directory. The `adapter_file` is explained a little more in the trimming command section of this worksheet)
- Note that although in this case you already manually made these directories, it's useful to make sure that the output directories exist by calling `mkdir` within the sbatch script. Look up what the `-p` parameter does.

  ![Trimmomatic variables](md_images/trimmomatic_variables.png)
 

Add commands to load the required modules for running this pipeline. We will be using FastQC and the trimming program Trimmomatic. If you are not sure which version of a program is available on the cluster you can save and quit out of Vim, then use the command `module spider <program>` to find the available versions.
- Technically Trimmomatic has a module on the AWS, but since it's run through Java, we need to use the path to the executable, which makes loading the module redundant. Here we've included a variable that stores the path to the Trimmomatic executable.
 
For the meat of the script, we will be running 3 analysis steps. Several analysis steps run together is called a pipeline. Our pipeline allows us to:
1. Run FastQC on raw sample fastq files
2. Trim the fastq files to remove adapters and other problematic sequence elements
3. Reevaluate the quality of the trimmed fastq files with FastQC.

  ![Trimmomatic pipeline](md_images/trimmomatic_pipeline.png)

The FastQC steps should look familiar to you from yesterday. Now let's look closer at the trimming step.
- In this script we are running paired end reads. Trimmomatic can be used on both single end or paired end reads. When setting your parameters use the appropriate format.
- Below is the syntax needed to run Trimmomatic:
  - For single-end reads\
    ```
    java jar <trimmomatic_executable_file> SE [ -threads <n> ]  [ -phred33 | -phred64 ] [ -trimlog <output_trimlog> ] <input_file> <output_file> ILLUMINACLIP
    ```
  - For pair-end reads\
    ```
    java jar <trimmomatic_executable_file> PE [ -threads <n> ]  [ -phred33 | -phred64 ] [ -trimlog <output_trimlog> ] <input_file1> <input_file2> <output_fileP1> <output_fileU1> <output_fileP2> <output_fileU2> ILLUMINACLIP
    ```
  - In the script, you'll see this command is split into multiple lines. Recall that the `\` at the end is used to break the code up for ease of reading. Depending on your Vim color settings, if `\` does not change color compared to the rest of the line, you may have an extra space after the `\`. Remove that space or your code will not run properly.

- The final parameter shown in the Trimmomatic commands above is the `ILLUMINACLIP` parameter, which should be replaced by a series of strings that define the trimming behavior you need. Below are several possible trimming strategies that can be implemented, singly or in combination:
  - `ILLUMINACLIP:<path_adapters_fasta>:<seed_mismatches>:<palindrome_clip_threshold>:<simple_clip_threshold>`
  - `LEADING:<quality>`
  - `TRAILING:<quality>`
  - `SLIDINGWINDOW:<window_size>:<required_quality>`
  - `MINLEN:<length>`

Each of the trimming behaviors is described in much more detail in the Trimmomatic manual. Find the manual online and read more about the two we're using in this script, `ILLUMINACLIP` and `CROP`.
- Note: The script provided a variable that stores the path to the adapter FASTA file. Take a look at this file to see its format.

After saving all of your changes to the script, run it!

- Submit the job (`sbatch <scriptname>`). The job manager will assign a job id to your run.
  - This pipeline has more tasks than yesterday's worksheet, so you will want to check the status of your job (`squeue -u <username>`) to see if the job is configuring (`CF`), pending (`PD`), running (`R`), or complete (not present on queue anymore).
  - If you don't get your expected outputs, or if the files are empty, troubleshoot using your `.err` file (located at the path in your SBATCH header). Again, you can view this file by using `more`, `less`, `cat`, or opening it in Vim.
