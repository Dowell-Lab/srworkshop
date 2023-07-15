# Short Read Day 3: Working with Supercomputers
- Author: Zach Maas

## Part 1: Write a script using slurm metadata that can

The goal of this exercise is to take the skills that you practiced yesterday and use them to write a script formatted so that it can run using the slurm workload manager. Slurm is software that allows for large compute jobs to be delegated across a cluster of computers, scheduling things so that resources are used optimally. This makes it so that many people can use the same compute system at the same time. Important: If you write a script for slurm and run it on the login node instead of submitting it to slurm, you will slow down the system for everyone else, potentially rendering it unusable. Double check your commands.

Things you need to do:

1. Make a folder on the short read system called workshop-day3
2. Make folders there for data, scripts, and results
3. Write a slurm script to download the fastq file located here:
`ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/sequence_read/SRR062641.filt.fastq.gz`
4. Submit your slurm script using the sbatch command and check to see that it does what you expect
5. Use the `md5sum` command to make sure that the checksum for the file matches `412d95d53ac858dbc899494f9db42ed2`. A checksum is a way of verifying that a file is the same as the original we wanted to download. Why would we want to use this?
5. Back up your script and results to your home directory

Some useful things to know:

- Your sbatch script needs additional directives added to the top of the file so that the slurm system knows how many compute resources you need:
```
#!/bin/bash
#SBATCH --output=<a path to a directory>/%x_%j.out
#SBATCH --error=<a path to a directory>/%x_%j.err
#SBATCH --mail-user=<an email address to get notifications at>
#SBATCH -p <the name of the slurm partition we’re on>
#SBATCH -N <the number of nodes to use, should be 1>
#SBATCH -c <the number of cores to use>
#SBATCH --mem=<the amount of memory to use, formatted 256mb, 8gb, etc>
<the rest of your script>
```
- wget is a command that will download a file from a remote location. The -c flag stands for (c)ontinue and means that your download can resume if it is interrupted for some reason.

## Part 2: Extend your slurm script to run fastqc on the downloaded file

The goal of this exercise is to extend the script you just wrote to run fastqc, a quality control tool, on the fastq file that you’re downloading. This requires that you use the module system available on the server to load fastqc and then to download the results from the server to your local machine so that you can look at them.

Things you need to do:

1. Write or extend your script to run fastqc on the fastq file that you downloaded. You will need to use the module system to load fastqc first inside of your script, using fastqc version 0.11.5. Put the output of this fastqc run into your results folder.
2. Use rsync to copy the output files from your fastqc run to your local machine and open them in a browser. How do they look?

Some useful things to know:

- The module command allows you to load different programs for use in your script
- module spider <query> will find a program if available
- module load <program> will load a program to be run. You may need to
specify the version
- The rsync command allows you to not just move files around on a single machine, but
between multiple machines
- The -e flag allows you to specify that you’re using the ssh protocol
- The -P flag will show you progress
- The -r flag will copy directories and the files in them recursively
- Make sure you’re on the machine that you want to download to. Why can’t you
send a file from the server to your machine from the server while logged into the
server? Why do you have to download the file using your local machine instead?
- You can view some documentation on fastqc here:
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/. Often project pages like that will have helpful documentation for explaining what a tool does, although sometimes it's better to read the original paper instead
- It's essential to keep track of the specific versions of tools that you use so that your code and analysis can be reproduced by others
