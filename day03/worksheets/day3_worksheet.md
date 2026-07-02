# Short Read Day 3: Working with Supercomputers
- Author: Zach Maas
- Edited by Lynn Sanford, 2024; Malia Fredrickson, 2025

## Part 1: Write an SBATCH script and submit it to Slurm

The goal of this exercise is to take the skills that you practiced yesterday and use them to write a script formatted so that it can run using the slurm workload manager. Slurm is software that allows for large compute jobs to be delegated across a cluster of computers, scheduling things so that resources are used optimally. This makes it so that many people can use the same compute system at the same time.

**Important:** If you write a script for slurm and run it on the login node instead of submitting it to slurm, you will slow down the system for everyone else, potentially rendering it unusable.

**You MUST use the command** `sbatch` **to submit your slurm scripts.**

Things you need to do:

1. Make a directory in your home directory (`/Users/<username>/`) on the AWS called `workshop-day3`.

2. Navigate to your scratch directory (`/scratch/Users/<username>/`).

3. Make another `workshop-day3` directory here. Within that directory, create directories for `data`, `scripts`, `results`, and `eofiles` (error and output files). 

4. In the `scripts` directory, write a slurm script to copy the following file to your `data` directory:  
`/scratch/Shares/public/sread2025/day3/SRR062641.filt.fastq.gz`
    - **Note:** If you look for this file, you may notice that this fastq filepath is problematic, but that's intentional, so follow these instructions for now.
    - You can name this script whatever you want, but by convention we use a `.sbatch` suffix for slurm scripts.
    - Your sbatch script needs directives added to the top of the file so that the slurm system knows how many compute resources you need:
    ```
    #!/bin/bash
    #SBATCH --job-name=<JOB_NAME>                         # Job Name
    #SBATCH --mail-type=ALL                               # Mail (NONE, BEGIN, END, FAIL, ALL)
    #SBATCH --mail-user=<YOU@EMAIL.COM>                   # Your email address
    #SBATCH --nodes=1                                     # Number of nodes requested, should always be 1
    #SBATCH --ntasks=1                                    # Number of CPUs (processor cores/tasks)
    #SBATCH --mem=256mb                                   # Memory limit, formatted 256mb, 8gb, etc
    #SBATCH --time=00:10:00                               # Time limit hrs:min:sec
    #SBATCH --partition=short                             # Partition/queue requested on server
    #SBATCH --output=/scratch/Users/<YOURUSERNAME>/workshop-day3/eofiles/%x_%j.out
    #SBATCH --error=/scratch/Users/<YOURUSERNAME>/workshop-day3/eofiles/%x_%j.err

    <the rest of your script: a command to copy the file>
    ```

5. Submit your slurm script using the `sbatch` command (`sbatch <scriptname>`) and check the queue with `squeue`.
    - You likely won't see your user ID in the queue, since this script will take less than a second to run, but it is always a good habit to check the queue after submitting a script.

6. Check if your file correctly copied to your `data` directory.

7. If the file isn't currently in your `data` directory (spoiler alert - it's not), look at your output and error files to troubleshoot.
    - Navigate to your `eofiles` directory and make sure that you have `.err` and `.out` files in there.
        - If you don't, then you have a problem with the `--error` and `--output` paths in your slurm script and need to make sure that those paths are valid. Fix the paths and submit the script again with `sbatch`.
        - If you do, look at the files by `cat`,`more`,`less`, or by opening them in Vim. Do you see an error? What do you think it means?

8. As it turns out, I 'forgot' to download this fastq file, and it does not exist anywhere on `/scratch/`, but I can tell you where it is available online. Modify your slurm script to download it from here:
`ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/sequence_read/SRR062641.filt.fastq.gz`
    - `wget` is a command that will download a file from a remote location. The `-c` flag stands for (c)ontinue and means that your download can resume if it is interrupted for some reason.

9. Submit your modified script with `sbatch`, check the queue with `squeue` (though this still runs in a few seconds), and make sure that the fastq file successfully downloaded into your `data` directory.
    - If you don't see the downloaded file and your job is not still running, go back into your `eofiles` directory and check your most recent error file. This is the error file with the highest job ID number. Look at the error within to see if you can understand what went wrong.

10. Use the `md5sum` command to make sure that the checksum for the downloaded file matches `412d95d53ac858dbc899494f9db42ed2`. Matching the checksum for a copied/downloaded file with the checksum for the original is a way of verifying that the file is exactly the same. Why would we want to use this?

11. Back up your script and data to the `workshop-day3` folder in your home directory.


## Part 2: Extend your Slurm script to run FastQC on the downloaded file

The goal of this exercise is to extend the script you just wrote to run FastQC, a quality control tool, on the fastq file that you’ve downloaded. This requires that you use the module system available on the server to load fastqc, run the `fastqc` command, and then download the results from the server to your local machine so that you can look at them.

Things you need to do:

1. Write or extend your previous SBATCH script to run FastQC on the fastq file that you downloaded.
    - If you're not still there, navigate back to your `/scratch/Users/<username>/workshop-day3/scripts/` directory
    - Inside your script (either the same one or a new one - up to you), you will need to use the module system to load fastqc, using fastqc version 0.11.5.
        - The `module` command allows you to load different programs for use in your script.
        - `module spider <query>` will find a program if available.
        - `module load <program>` will load a program to be run. You may need to specify the version if there is more than one version on the cluster. Even if there is only one version, specifying the version in your script documents it for later use/methods write-ups. It's essential to keep track of the specific versions of tools that you use so that your code and analysis can be reproduced by others.
        - **IMPORTANT:** In our case, `module spider` or `module avail` will SHOW more than one version of FastQC for legacy reasons, but there is only actually one version installed - `0.11.5`.
    - If you need help with the fastqc command syntax, outside of your script you can load the module on the command line, then access the command documentation (`fastqc --help`).
        - Unfortunately this is one case where the official online documentation doesn't help that much with running it on the command line. Usually you should take a look at the online documentation for most software packages to make sure you're running them correctly for your use case.
    - Once you've run the sbatch script with your fastqc command (again, run it with `sbatch`), make sure you have a zip file and an html file as output. By default the results are written to the same directory as the data. Move the output of this fastqc run into your results folder, then back it up in your home directory (`/Users/<username>/`) day 3 folder.
        - Optional: using the command line fastqc documentation, figure out how to write the results directly to your `results` directory.
    - Remember, you can view the whole queue with `squeue`, and only **your** jobs with `squeue -u <username>`.

2. Use `rsync` to copy the output files from your FastQC run to your local machine and open them in a browser.
    - The `rsync` command allows you to not just move files around on a single machine, but also between multiple machines.
        - The `-P` flag will show you transfer progress.
        - The `-r` flag will copy directories and the files in them recursively.
    - Make sure you’re on the non-server machine to which you want to download. In this case, it's your local laptop. This will be true whether you are transferring FROM the server or TO the server.
        - Use `hostname` and your prompt color as clues to your current machine.
        - It's best to open up a new terminal session for transferring so that you can still copy and paste paths from the server.
        - Why can’t you send a file from the server to your machine from the server while logged into the server? Why do you have to download the file using your local machine instead?
    - On your local machine, create and/or navigate to the folder where you want to transfer your file.
        - Windows users: This isn't as straightforward as it is on Mac. Please read through <a href="https://github.com/Dowell-Lab/srworkshop/blob/main/resources/Windows_file_locations.md">this information</a> to understand the filesystem.
    - The inter-machine syntax of the `rsync` command is as follows:\
    `rsync <username>@<server IP address>:<full path to file> <destination>`
        - Don't forget the colon in between the server address and the filepath
        - If you are currently in your destination directory, `<destination>` can be replaced with `.`
    - How does this library look, according to the quality control? Can it be analyzed further?
