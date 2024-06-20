# Day 4 Worksheet - Read mapping and visualization
Author: Qing Yang, 2021
Edited: Lynn Sanford, 2024

1. Make sure that you have two trimmed fastq files from the Trimmomatic step, one for read 1 and one for read 2. If you're still struggling to get Trimmomatic to work, retrieve these two files from `/scratch/Shares/public/sread2024/data_files/day4/cookingshow/`.\
![Trimmed fastq file locations](md_images/trimmed_fastq_file_locations.png)

2. Create a new directory within your working day4 directory called `hisat2`.\
![Day4 directories](md_images/Day4_directories.png)

3. Copy the `d4_mapping.sbatch` script from the GitHub repo to your `day4/scripts/` directory, then open it in Vim.

4. Edit the SLURM header configurations:
    - Change the name of the job to something useful like `hisat2_mapping`
    - Replace any user inputs (`<>`) with your specific values. Make sure your output and error directory exists.
    - HISAT2 can use multiple processors per input file, but still only one node. For this run, edit the fields to reserve 1 node, 8 processors (`ntasks`), 8 Gb memory and 1 hr of walltime.

5. Assign path variables. In this case, `DATADIR` is set to your working directory, and the `TRIM` path variable is the default trimmomatic directory within that. The `HISAT` variable points to the directory you made at the beginning of the worksheet and will be your output directory for this script. Make sure all of these paths are accurate for you.

6. Load modules. HISAT2 is our read mapping software and we'll use SAMtools for file conversion afterward.

7. Finally, look at the read mapping and file conversion commands.
    - As before, note that `\` with no trailing spaces provides a displayed line break that doesn't actually break up the command as run by the computer. This makes the code more readable.
    - The genome index, which is specific to HISAT2, is located at the `-x` path. If you look at this directory on the command line, you'll see that the last term in that path (`genome`, here) should correspond to the prefix for the files in that directory.
    - Make sure you understand the HISAT2 and SAMtools commands we're using here. Look at the online documentation for these tools for more information. One exception is that the HISAT2 parameter `--very-fast` is not well-documented in the HISAT2 manual and actually is a holdover from a previous generation mapping software called Bowtie2 - the Bowtie2 documentation has more information for that parameter.

8. Save your edits and exit Vim, then run the job with the `sbatch` command. You can check your command status with `squeue -u <USERNAME>`.

9. You should have 5 output files in your .../hisat2/ directory when your job is complete. If you don't, check your error file to troubleshoot.

10. To visualize the mapped reads using IGV, you will need to transfer the `...sorted.bam` AND `...sorted.bam.bai` files to your local machine.
    - Copy the full path to your `...sorted.bam` file. One way to display this full path is using the command `realpath <file>`. Then copy it to your clipboard.
    - Open a terminal on your local machine (NOT connected to the AWS).
    - Navigate to the local directory where you want your mapped files.
      - Note that on a Windows machine, it may be easiest to navigate to your Desktop, which is either located at `/mnt/c/Users/<username>/Desktop/` OR `/mnt/c/Users/<username>/OneDrive/Desktop/`
    - Use `rsync` to transfer the bamfile to your current directory on your local machine (`rsync <username>@<AWS_ip_address>:<path to file> ./`)