# Workshop Day 2: Navigating the Unix Environment
Author: Christopher Rauchet (2025)<br>
Edited: Lynn Sanford (2026) 

# Part 1: Working with Files and Directories in the Command Line

The objective of this worksheet is to:
- Learn how to navigate the terminal
- Create a hierarchical directory structure
- Copy files to a directory on the server and your local device
- Work with and manage files

You are welcome and encouraged to try this autonomously with the help of cheat sheets. You can also refer to the step by step instructions for Part 1 below, which contain example code and rationale for each step in the terminal. TAs are available to assist.

> The file we will be using for this worksheet is a GTF, or Gene Transfer Format, which contains gene annotations for the human genome. Each annotation line specifies a set of genomic coordinates along with metadata about the gene. We will work with these files more in week 2.

> IMPORTANT: We'll go over this more tomorrow, but you have TWO main directory trees associated with your account. One is at `/Users/<username>` and is referred to as **your home directory**. On a real cluster this would be backed up, but it's slow for computation. The second is at `/scratch/Users/<username>` and is referred to as **your scratch directory**. This would not be backed up but is where you want to do most of your work. For today's exercises you will only work in your scratch directory, but later you will move back and forth.

## Part 1 Tasks (Try on your own - if you need help, guided steps are below)

1. In `/scratch/Users/<username>`, create a directory `workshop-day2` and subdirectories `scripts`, `eofiles`, `gtf`
2. Copy the gene annotation file `/scratch/Shares/public/sread/data_files/day2/bedfiles/hg38.genes.gtf`  to your `gtf` directory
- Note: DON'T open this file in VIM - it's large and you don't want to accidentally edit it. Do the following steps outside of VIM.
3. Determine the number of annotations in the file
4. Determine the number of annotations from chromosome 1 in the file
5. Create a new file in the `gtf` directory called `chr1_exon.gtf` that contains the gene exon annotations from chromosome 1
6. Copy the `chr1_exon.gtf` file to your local computer

## Guided Steps Part 1
### A.  Create the Directories

1. Start by logging into the AWS using `ssh` followed by `<username>@<IP>`
2. Once on the AWS, change directories using `cd` to `/scratch/Users/<username>` 
3. You can check the contents of this directory using the `ls` command
4. While in your scratch directory make a new directory `workshop-day2` using `mkdir`
5. To check if you made the new directory, you can type `ls` while in your scratch and you should see `workshop-day2` appear
6. Change directories to `workshop-day2` using `cd` 
    - Note: after typing `cd` you should start typing the name of the directory and hit `Tab` to tab complete to save you time and typing
7. Once in `workshop-day2` make the following directories using the `mkdir` command: `scripts`, `eofiles`, `gtf`
    - Note: when making multiple directories within the same directory you can use `mkdir <dir1> <dir2>`
    > Challenge: If you use `mkdir` to make a directory that already exists, it will throw an error. Use the `man` command to see what flag you can add to make it so that `mkdir` will only make a directory if it does not exist
8. You can verify the new directories are created using `ls`
9. Practice navigating around between the directories using relative and absolute paths

### B. Copy over the GTF File

1. Navigate to the `/scratch/Users/<username>/workshop-day2/gtf` directory
    - This is the directory into which we will copy the gtf file
2. Using the `rsync` command, copy the `hg38.genes.gtf` file (`/scratch/Shares/public/sread/data_files/day2/bedfiles/hg38.genes.gtf`) to this directory
    - The format of the `rsync` command is as follows: `rsync <source> <destination>`
    - We are currently in the destination directory. What is the symbol we can use to represent that?
    > Challenge: You do not need to be in the destination directory to rsync a file to that location. What is the command to move the file to the destination directory if you are not in that folder? (hint: you will need to specify a path)
3. Check that the file is there using `ls`
    - Do NOT open this file in VIM (see below)

### C. Examine and manipulate the GTF file
Once a data file is in the place you need it to be, there are a few things to note. Firstly, to maintain responsible conduct of research, you should not open this file using VIM, because editing raw data is bad! Also, opening large data files in VIM can be way too memory-intensive. Instead, we are going to use different tools to view and engage with the data: `head`, `tail`, `less`, `more`. 

1. To look at the format of the gtf file, move to the gtf directory and use `head <filename>` to print the first 10 lines of the file
    - Don't forget about tab complete!
    - Notice the format of the file - it is a tab separated variable format where there is a tab between each column
    - Each row is a unique gene annotation, containing information about the chromosome, type of annotation (exon, etc), location and more.
2. Examine the contents of the file with `more <filename>` and `less <filename>`. Both of these commands bring up an interactable interface rather than just printing lines to your screen. You can exit this interface with `q` or the `Esc` key
    - You cannot edit files with these commands
    - These commands also *don't* read the whole file into memory - they only read it in line by line, so you can run them on the head node
3. Now, let's figure out how many annotations are in this file. Because each line is an annotation, we just need to calculate the number of lines to know the number of annotations. To do so we are going to use the word count (`wc`) command, with an additional flag. Use the `man` command on `wc` to figure out what flag tells you how many lines are in the file
    - Once you find the flag use the `wc <flag> <filename>` format to find the number of lines in the file
    > Challenge: Change directories, and try writing the same command except with a relative path and/or absolute path to the file.
4. Warning: Clearly, there are many lines! if you were to print out this entire file, that would take a long time to load on your screen. If you accidentially print out the whole file, cancel the command by hitting `ctrl-C`
    > Now, often times a large data set like this might be redundant for your purposes. For example, maybe you would only need the annotations on a certain chromosome, or just the exons, or both. We are now going to take this file and extract certain information. 
    
    > Stop and Think: How could you write a script to figure out how many annotations are from chromosome 1 only? What are the commands we will need? How about the flags? How do we combine two commands together?

5. To figure out how many annotations are from chromosome 1, we will need to combine `grep` and `wc` with additional flags
    - `grep` is a command that can print lines that contain a specific pattern using `grep <pattern> <file>`.  Intuitively, to extract the lines with chromosome 1, one might write the following line of code `grep chr1 <filename>` but there are serious issues with this line of code. Stop and think what would could go wrong here? 
        - If you run that line of code, first of all, `grep` would print out the output to your screen, and you do not need millions of lines printed in your terminal!
        - Additionally, if you just gave `grep` *chr1* it would print out every line that contains the pattern anywhere, including *chr10*, *chr11*, etc.. Check the `man` page for `grep` to find the flag that will only print the lines that contain the term as a whole word, so that when we say *chr1*, we ONLY print lines with *chr1*
    - To combine `grep` with `wc` you will need to directly pipe (`|`) the output of `grep` to `wc`
6. Pipe `grep` output to `wc` to find out how many annotations are from chromosome 1 using a combined command which takes into account the necessary flags: `grep <flag> <pattern> <file> | wc <flag>`
    > Challenge: How many annotations on chr6 are there?

    > Now, let's say that you don't just want to know how many annotations there are for a chromosome, but instead you want to create a file that has the specific annotations that you are looking for. As an example, create a file that has only the annotations of exons from chromosome 1. Name this file `chr1_exon.gtf`

    > Stop and Think: How can you use `grep` and `|` to get lines that contain only *chr1* and *exon*? How do you output that to a new file?

7. While still in the `gtf` directory, write a line of code to take the annotations only from exons on chr1 to an output file based on this: `grep <flag> <pattern> <filename> | grep <flag> <pattern> > <outputfile>.gtf`
    - Note: The extra `>` above is not a typo. This character writes the output of the second `grep` command to a file
    - Think about why we do not need a filename for the second `grep` in that pipe
    - Be sure to include a `.gtf` at the end of the output file so that you can remember the file type

8. After you create the file, use `ls` in the directory to see if it is there
9. Use `head` to check the first 10 lines and see if they are only chr1 exons, and `tail` to see what the last lines look like

### D. Copy a file from the server to your local computer

Once you have data files, it is important to be able to transfer them to your local computer so you can then use them for things like viewing them on IGV.

Stop and Think: Based on the example in class, how would you copy the `chr1_exon.gtf` file to your local device? Do you want to be logged into the server when you do this or on the local computer? Do you want the file to end up in a particular place?

1. When you are ready to copy the file, open a new window in the terminal or log out of the server. You want to be on your local computer when copying a file from the AWS server - the IP address for the AWS is fixed and that of your laptop is not
2. Change directories on the local computer to the place where you want the file to be, or have the relative or absolute path for the destination ready
3. Use the `rsync` command with the following format: `rsync <username>@<awsip>:</server_absolute_path_to_file> <local_destination_directory_path>`
    - To get the absolute path of the file on the AWS, you can type it all out or while still logged onto the aws you can use `realpath <file>` and copy that output. This path MUST begin with `/` within the `rsync` command
    - You must have the `:` separating the server address from the file path
    - Windows paths are sometimes confusing. You can refer to [this worksheet](https://github.com/Dowell-Lab/srworkshop/blob/main/resources/Windows_file_locations.md) to get more information about Windows paths.
4. Check to see if the file appeared on your local computer using `ls`

# Part 2: File Management and Permissions

Now that you have worked through a series of steps and have a better sense of navigating the terminal, we are going to work on a new skill - using manual (`man`) information to identify how to create commands on your own. We are here to help of course, but the most common way that bioinformaticians learn is by having to dive into manuals yourself!

## Background
In computational work, file management is an important part of sharing and protecting data. A common file management technique is *compression*. Compression removes any redundancy in files to make them smaller. You may have already worked with compressed `.zip` files or archives.

Different algorithms can be used for compression. `zip` is one, but it's not as commonly used on Linux as the `gzip` algorithm. Using `gzip` automatically adds a `.gz` to the file, indicating it is zipped. The command `gunzip` takes the compressed file and decompresses it to the original file.

Another file management technique is *archiving*. Archiving is the process of gathering up many files and building them together into a single large file, often done as a part of backing up a system or making many files easily portable. On Linux, you use `tar` to archive files or extract archived files. A file that ends in `.tar` is an archived file. These are often compressed into `.tar.gz` files.

An additional important concept for file management is modifying *permissions*. On your computer, you probably run into files that are read-only. This means that you don't have permission to write to the file. On compute clusters, permissions are essential to determine who has access to what files and what they can do with those files. You also may want to change permissions on files to make them unable to be changed or deleted. You can only change permissions on files you 'own', i.e. have created yourself from scratch or as a copy of another file.

## Objective
You will read about `.tar.gz` files, copy a compressed archive file into your directory, and decompress and extract the files from it. Then you will change the permissions on one of those files. 

Steps:
1. In the `workshop-day2` directory you made in your scratch, create a subdirectory called `bedfiles`
2. Copy this file to your `bedfiles` directory: `/scratch/Shares/public/sread/data_files/day2/bedfiles/chr1_bedfiles.tar.gz` 
3. Using the resources below, decompress and extract the files from `chr1_bedfiles.tar.gz`
4. Change the permissions on any one of the files so that the owner has read, write and execute permissions; the group has read and execute permissions; and the public has only read permission
5. Copy one of the files to your local computer (for more practice)


**Some useful commands and tools:**
- Basic Utilities:
    - `man` (manual) followed by the name of a command will tell you how to use a  command
        - This might seem overwhelming at first because it includes a lot of information, but getting used to learning how to read these manuals is important
        - Type `q` or hit `Esc` to exit from man mode
    - `ssh` lets you log into a remote server like the AWS
    - `whoami` will print your current username
    - `hostname`: This is useful to make sure you’re actually on the short read server and not your local machine
    - `mkdir`: make a folder
- Navigating Around:
    - `pwd` (print working directory) will tell you which directory you’re in
    - `cd` (change directory) lets you move to a different directory. You can navigate using absolute or relative paths. For example, if you’re starting in `/scratch/Users/<username>` and you want to be in the `workshop-day2` directory inside of it, typing `cd /scratch/Users/<username>/workshop-day2` will move you using the *absolute path* and typing `cd workshop-day2` will move you using the *relative path*
        - Relative paths have some shortcuts: `..` (two dots) means "the directory immediately above the one that I'm in" and `.` (one dot) means "the directory that I'm currently in"
    - `ls` (list) prints out information about the files in your current directory
- Working with Files:
    - `mv` (move) moves a file from a source path to a destination path
    - `cp` (copy) copies a file from a source path to a destination path
    - `rm` (remove) removes a file or directory. **WARNING**: there is no recycle bin on a Unix system. If you remove a file, it is gone forever with no chance of recovery. **DOUBLE CHECK YOUR COMMANDS**
    - `gzip`/`gunzip` allows you to create and extract `.gz` zipped files
    - `tar` (tape archive) allows you to create and extract archive files that have extensions of `.tar`, `.tar.xz`, `.tar.gz`, and more. (see here for more information: https://phoenixnap.com/kb/extract-tar-gz-files-linux-command-line)
    - `less` allows you to view files on the command line without editing them
    - `cat` (concatenate) prints out the contents of a file directly to your terminal
    - `head`/`tail` show you the first and last lines of a file, respectively (10 by default). These are useful for quickly making sure a file is in the right format
    - `vim` is a text editor that allows you to edit files for writing scripts. `nano` is a simpler text editor for writing scripts on the command line (but won't be covered in this class)
    - `chmod` changes the permission on files

 
# Part 3 (optional challenge): Write a script

One of the reasons that we have taught you how to use VIM is because we use VIM to write scripts. A script is a series of commands that the shell can execute sequentially. So rather than typing commands directly into the command line interface, we can create a text document that has the ability to be run. This will be the focus of tomorrow, but you can get a head start now if you'd like.

To indicate something is a bash shell script the first line must be `#!/bin/bash`. On subsequent lines, you can write commands that are executable by the bash shell.

To submit a bash script, you will need to exit out of the script, and on the command line interface use the `bash` or `sh` command. 

*Goal*: Learn how to take a series of steps like you ran before and combine them into a single script that you can run for reproducibility.

### Steps:
1. Make a new file called `sync_files.sh` in your scratch `workshop-day2/scripts` directory using `vim`
2. Write `#!/bin/bash` as the first line of the file
3. Write a series of commands on separate lines in your script to make a subdirectory in `workshop-day2` called `gtf_from_script` and copy the same human annotation gtf file from Part 1 to this new directory. As a tip, absolute paths in scripts make them much more reproducible and reduce errors
4. Run the script on the command line interface using `sh <script_file>` to submit the script
5. Check to see if the directory is made, and if the file copied over















