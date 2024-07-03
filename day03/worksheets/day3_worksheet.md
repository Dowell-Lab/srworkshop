# Short Read Day 3: Working with Supercomputers
- Author: Zach Maas
- Edited by Lynn Sanford, 2024

## Part 1: Write an SBATCH script and submit it to Slurm

The goal of this exercise is to take the skills that you practiced yesterday and use them to write a script formatted so that it can run using the slurm workload manager. Slurm is software that allows for large compute jobs to be delegated across a cluster of computers, scheduling things so that resources are used optimally. This makes it so that many people can use the same compute system at the same time.

Important: If you write a script for slurm and run it on the login node instead of submitting it to slurm, you will slow down the system for everyone else, potentially rendering it unusable.

**You MUST use the command** `sbatch` **to submit your slurm scripts.**

Things you need to do:

1. Make a directory in your home directory (`/Users/<username>/`) on the AWS called `workshop-day3`.

2. Navigate to your scratch directory (`/scratch/Users/<username>/`).

3. Make another `workshop-day3` directory here. Within that directory, create directories for `data`, `scripts`, `results`, and `eofiles` (error and output files). 

4. In the `scripts` directory, write a slurm script to copy the following file to your `data` directory: `/scratch/Shares/public/sread2024/day3/SRR062641.filt.fastq.gz`
    - By convention, we use a `.sbatch` suffix for slurm scripts.
    - Your sbatch script needs directives added to the top of the file so that the slurm system knows how many compute resources you need:
    ```
    #!/bin/bash
    #SBATCH --output=<path to your eofiles directory>/%x_%j.out
    #SBATCH --error=<path to your eofiles directory>/%x_%j.err
    #SBATCH --mail-user=<an email address to get notifications at>
    #SBATCH --job-name=<a descriptive name for the job>
    #SBATCH -p <the name of the slurm partition we’re on>
    #SBATCH -N <the number of nodes to use, should always be 1>
    #SBATCH -c <the number of cores to use>
    #SBATCH --mem=<the amount of memory to use, formatted 256mb, 8gb, etc>

    <the rest of your script: a command to copy the file>
    ```
    - Note: If you're paying close attention, you may notice that this fastq filepath is problematic, but that's intentional, so follow these instructions for now.

5. Submit your slurm script using the `sbatch` command (`sbatch <scriptname>`) and check the queue with `squeue`.
    - You likely won't see your user ID in the queue, since this script will take less than a second to run, but it is always a good habit to check the queue after submitting a script.

6. Check if your file correctly copied to your `data` directory.

7. If the file isn't currently in your `data` directory (spoiler alert - it's not), look at your output and error files to troubleshoot.
    - Navigate to your `eofiles` directory and make sure that you have `.err` and `.out` files in there.
        - If you don't, then you have a problem with the `--error` and `--output` paths in your slurm script and need to make sure that those paths are valid. Fix the paths and submit the script again with `sbatch`.
        - If you do, look at the files by `cat`,`more`,`less`, or by opening them in Vim. Do you see an error?

8. As it turns out, I 'forgot' to download this fastq file, and it does not exist anywhere on `/scratch/`, but I can tell you where it is available online. Modify your slurm script to download it from here:
`ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/sequence_read/SRR062641.filt.fastq.gz`
    - `wget` is a command that will download a file from a remote location. The `-c` flag stands for (c)ontinue and means that your download can resume if it is interrupted for some reason.

9. Submit your modified script with `sbatch`, check the queue with `squeue`, and make sure that the fastq file successfully downloaded into your `data` directory.
    - If you don't see the downloaded file and your job is not still running, go back into your `eofiles` directory and check your most recent error file. This is the error file with the highest job ID number. Look at the error within to see if you can understand what went wrong.

10. Use the `md5sum` command to make sure that the checksum for the downloaded file matches `412d95d53ac858dbc899494f9db42ed2`. Matching the checksum for a copied/downloaded file with the checksum for the original is a way of verifying that the file is exactly the same. Why would we want to use this?

11. Back up your script and data to the `workshop-day3` folder in your home directory.


## Part 2: Extend your Slurm script to run FastQC on the downloaded file

The goal of this exercise is to extend the script you just wrote to run FastQC, a quality control tool, on the fastq file that you’ve downloaded. This requires that you use the module system available on the server to load fastqc, run the `fastqc` command, and then to download the results from the server to your local machine so that you can look at them.

Things you need to do:

1. Write or extend your previous SBATCH script to run FastQC on the fastq file that you downloaded.
    - Navigate back to your `/scratch/Users/<username>/workshop-day3/scripts/` directory
    - Inside your script, you will need to use the module system to load fastqc, using fastqc version 0.11.5.
        - The `module` command allows you to load different programs for use in your script.
        - `module spider <query>` will find a program if available.
        - `module load <program>` will load a program to be run. You may need to specify the version if there is more than one version on the cluster. Even if there is only one version, specifying the version in your script documents it for later use/methods write-ups. It's essential to keep track of the specific versions of tools that you use so that your code and analysis can be reproduced by others.
    - If you need help with the fastqc command syntax, outside of your script you can load the module on the command line, then access the command documentation (`fastqc --help`).
    - Copy the output of this fastqc run into your results folder, then back it up in your home directory (`/Users/<username>/`) day 3 folder.
    - Remember, after your SBATCH script is modified, submit it to the queue with `sbatch`. You can view the whole queue with `squeue`, and only **your** jobs with `squeue -u <username>`.

2. Use `rsync` to copy the output files from your FastQC run to your local machine and open them in a browser. How do they look?
    - The `rsync` command allows you to not just move files around on a single machine, but also between multiple machines.
        - The `-e` flag allows you to specify that you’re using the ssh protocol.
        - The `-P` flag will show you transfer progress.
        - The `-r` flag will copy directories and the files in them recursively
    - Make sure you’re on the non-server machine to which you want to download. In this case, it's your local laptop. This will be true whether you are transferring FROM the server or TO the server.
        - Use `hostname` and your prompt color as clues to your current machine.
        - It's best to open up a new session for transferring so that you can still copy and paste paths from the server.
        - Why can’t you send a file from the server to your machine from the server while logged into the server? Why do you have to download the file using your local machine instead?
            - Hint: IP addresses
    - On your local machine, create and/or navigate to the folder where you want to transfer your file.
        - Unlike on Mac, Windows machines have the Lunix virtual machine set up in an obscure part of the file structure.
        - On Windows 11, you can usually navigate to your Linux files in File Explorer.
        - On Windows 10 and some Windows 11 systems, the easiest way to put your files somewhere accessible is to navigate to your Desktop within the terminal.
            - Upon opening a terminal, `pwd` will show you `/home/<linux username>`. You can then `cd` to either `/mnt/c/Users/<computer username>/Desktop/` OR `/mnt/c/Users/<computer username>/OneDrive/Desktop/` (depending on your version of Windows). Running your `rsync` command with this as the destination will transfer files to your Desktop.
    - The inter-machine syntax of the `rsync` command is as follows:\
    `rsync <username>@<server IP address>:<path to file> <destination>`
        - Don't forget the colon in between the server address and the filepath
        - If you are currently in your destination directory, `<destination>` can be replaced with `./`