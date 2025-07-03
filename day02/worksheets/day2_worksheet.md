# Workshop Day 2: Navigating the Unix Enviorment
Author: Christopher Rauchet (2025) 
## Part 1: Working with Files and Directories in the Command Line


The objective of this worksheet is to learn how to navigate the terminal, create a hierarchical directory structure, copy files to a directory in the server and the local device, and work with and manage the files.


The file we will be using is a GTF, or Gene Transfer Format, which contains gene annotations for the human genome. Depending on your comfort level or experience, you are welcome to try this autonomously, or refer to the step by step instructions for Part 1 below, which contain screenshots of the code and rationale for each step in the terminal. The more you can do independently the better, but if you need help you can rely on the steps below, and your TAs.

Part 1 Steps (Try it on your own if you feel comfortable, if not jump to A)

1. In your /scratch/Users/<username> create a directory *day2-workshop* and subdirectories *scripts*    *eofiles*    *gtf*
2. Copy the gene annotation file */scratch/Shares/public/sread2025/data_files/day2/bedfiles/hg38.genes.gtf*  to your *gtf* directory
3. Determine the number of annotations in the file
4. Determine the number of annotations from chromosome 1 in the file
5. Create a new file in the *gtf* directory, the contains the gene annotation only exons from chromosome 1, called *chr1_exon.gtf*
6. Copy the *chr1_exon.gtf* file to your local device

Guided Steps Part 1
### A.  Create the Directories

1. Start about by logging into the AWS using the `ssh` followed by `<username>@<IP>`
2. Once in the AWS, change directories using `cd` to `/scratch/Users/<username>` 
3. You can check the contents of this directory using the `ls` command
4. While in your scratch make a new directory *day2-workshop* using `mkdir`
5. To check if you made the new directory, you can type `ls` while in your scratch and you should see *day2-workshop* appear
6. Change directories to *day-2-workshop* using `cd` 
- Note: after typing `cd` you can start typing the name of the directoy and hit the tab key to tab complete, it will save you time.
7. Once in *day2-workshop* make the following directories using the `mkdir` command: *scripts*   *eofiles*   *gtf*
- Note: when making multiple directories within the same direcotry you can use `mkdir <dir1> <dir2>`
- Challenge: If you use `mkdir` to make a direcotory that already exists, it will override that directory and everything in it. Use the `man` command to see what flag you can add to make it so that `mkdir` will only make a directory if it does not exist
8. You can verify the new directories are created using `ls`
9. Practice navigating around between the directories using relative and absolute paths


### B. Copy over the GTF File

1. Navigate to the */scratch/Users/<username>/day2-workshop/gtf* directory
- This directory is the directory we will copy the gtf file into
2. Using the `rsync` command, copy the *hg38.genes.gtf file: */scratch/Shares/public/sread2025/data_files/day2/bedfiles/hg38.genes.gtf*  to your gtf directory
- The format of the `rsync` command is as follows: `rsync <source> <destination>`
- Recall, because we are currntly in the desintation directory what is the symbol we can use to represent that?
- Challenge: You do not need to be in the destination directory to rsync a file to there, what would the script be for moving the file if your current working directory were simply: */scratch/Users/<username>* ?
3. Check that the file is there using `ls`


### C. Working with the GTF file

Once a data file is in the place you need it to be, there are a few things to note. Firstly, to maintain responsible conduct of research, you should not open the file using vim, because editing raw data is bad! Instead, we are going to use different tools to view and engage with the data: head, tail, less, more. 

1. To look at the format of the gtf file, move to the gtf directory and use `head <filename>` to print the first 10 lines of the file.
- Don't forget about tab complete!
- Notice the format of the file, it is a tab seperated variable format where there is a space between each column
- Each row is a unique gene annotation, containing information about the chromosome, type of annotation (exon, etc), location and more.
3. Now, let's figure out how many annotation are in this file. Because each line is an annotation, we just need to calculate the number of lines to know the number of annotations. To do so we are going to use the word count command, with a flag. Use the `man` command on `wc` to figure out what flag tells you how many lines are in the file.
- Once you find the flag use the `wc -flag <filename>` format to find the number of lines in the file
- Challenge: Change directories, and try writing the same script excpet with a relative path as well as the absolute path to the file.
4. Warning: Clearly, there are many lines! if you were to print out this entire file, that would take a long time to load on your screen. If you accidentially print out the whole file, cancel the command by hitting the control key and c at the same time.

Now, often times a large data set like this might be redundant for your purposes. For example, maybe you would only need the annotations on a certain chromosome, or just the exons, or both. We are now going to take this file and extract certain information out of it. 

Stop and Think: How could we write a script to figure out how many annotations are from chromosome 1 only? What are the commands we will need? How about the flags? How do we combine two commands together?


4. To figure out how many annotations are from chromosome 1, we will need to combine `grep` and `wc` with the use of flags.
- `grep` is tool that can print lines that contain a certain line or symbol using `grep <term> <file>` .  Intuitively, to extract the lines with chromosome 1, one might write the following line of code `grep chr1 <filename>` but there are serious problems with this line of code. Stop and think what would would go wrong here? 
- If you typed that line of code, first of all, grep would print out the output to your screen, kind of like `cat` and we do not need thousands of lines of code printed to your computer!
- Addtionally, if we just gave grep *chr1* it would print out every line that contains that, so *chr10* *chr11* etc.. To combat this, use the `man` command on `grep` to find the flag that will only print the lines that contain the term as a whole word, so that when we say chr1, we only print lines with chr1 and not chr13 for example
- To combine it with `wc` we will need to directly pipe command with `wc`
6. Pipe `grep` and `wc` to find out how many annotations are from chromosome 1 using a script which takes into account the necessary flags: `grep -flag <term> <file> | wc -flag`
- Challenge: How many annotations of exons on chr6 are there?


Now, let's say that we don't just want to know how many annotations there are for a chromosome, but instead we want to create a file that has just the specific annotations that we are looking for. As an example, let's create a file that has only the annotations of exons from chromosome 1, to a file that we will name *chr1_exon.gtf*

Stop and Think: How can we use `grep` and `|` to get lines that contain only chr1 and exons? How do we output that to a new file?

6. While still in the *gtf* directory, write a line of code to take the annotations only from chr1 and exons to an output file based on this: `grep -flag <term> <filename> | grep -flag <term> > <outputfile>.gtf
- Note: Think about why we do not need a filename in the second `grep` in that pipe
- The key output indicator is '>'
- Be sure to include a .gtf at the end of the output file so that you can remember the file type.

8. After you create the file, use `ls` in the directory to see if it is there
9. Use `head` to check the first 10 lines and see if there is just chr1 and exon, and `tail` to see what the last lines look like. 

### D. Copying a file to from the server to the local device

Once you have data files, it is important to be able to transfer them to your local device so you can then use them for things like viewing them on IGV

Stop and ThinK: Based on the example in class, how would you copy the chr1_exon.gtf file to your local device? Do you want to be logged into the server when you do this or on the local computer? Do you want the file to end up in a particualr place?

1. When you are ready to copy the file, open a new window on in the terminal or log out of the server. You want to be on the local device when copying a file from the server.
2. Change directories on the local device to the place where you want the file to copy over, or have the realtive or absolute path for the destination ready.
3. Use the `rsync` command with the following format: `rsync <user>@<awsip>:</file_absolute_path_server> <destination directory path local>`
To get the absolute path of the file, you can type it all out or while still logged onto the aws you can use `realpath` and copy that to get the absolute path
4. Check to see if the file appeared on the local device using `ls` in that destination directory, and you can even try and view the file on IGV if that makes you happy!


## Part 2: File Management and Permissions

Now that you have worked through a series of steps and have a better sense of navigating the terminal, we are going to work on a new skill - using the manual to identify how to create commands on your own. Rather than have each step laid out one by one, you will have to create a series of steps independetly to try a new task. We are here to help of course, but the most common way that bioinformaticians learn is by having to dive into the manual yourself! For this task, you will independtenly navigate file management techniques in the terminal.

Background: In computational work, file management is an important part of sharing and protecting data. A common file management technique is *compressing files*. File compression gets removes any redunancy in a file to make the file smaller and take up less storage. A file can be compressed, or zipped, using `gzip`. Using `gzip` adds a .gz to the file, indicating it is zipped. The command `gunzip` takes the commpressed file and decpmresses it to the original file. An additional file managment technique is *archiving*. Archiving is the process of gathering up many files and building them together into a single large file, often done as a part of backing up a system or long term storage. To archive files into one large file, and extract the files - meaning to make the archived file go back to the multiple original files - you use `tar`. A file that ends in *.tar* is an archived file. To make things really interesting, we can compress an archived file to save even more space, those files end in *.tar.gz*

 
Objective: For this, you will read about .tar.gz files, copy an archived compressed file into your directory decompress and extract the files from it. Then you will change the permission on one of those files indepndently. 

Steps:
1. In the *day2-workshop* directory you made in your scratch, create a subdirectory called *bedfiles*
2. Copy this file to your *bedfiles* directory */scratch/Shares/public/sread2025/data_files/day2/bedfiles/chr1_bedfiles.tar.gz* 
3. Using the resources below, figure out how to decompress and extract the files from *chr1_bedfiles.tar.gz*, and do so
4. Pick one of the files and using the information from the video last night, change the permission on one of the files so that the owner has read, write and execute; the group has read and execute, and the others has read only permissions
5. Copy one of the files to your local device


**Some useful commands and tools:**
- Basic Utilities:
-  **man** (manual) followed by the name of a command will tell you how to use a  command
                - This might seem overwhelming at first because it includes everything but getting used to learning how to read this is really helpful
                - type q or click "esc" to escape from man mode
- **ssh** lets you log into a remote server like the one we’re using during the short  read workshop
- **whoami** will print your current active user
- **uname-n** This is useful to  make sure you’re actually on the short read server and not your local machine
- **mkdir** make a folder
- Navigating Around:
- **pwd**  (print working directory) will tell you which directory you’re in
- **cd** (change directory) lets you move to a different directory. You can move folders using absolute or relative paths. For example, if  you’re starting in /scratch/ and you want to be in the Workshop directory inside of it, typing cd /scratch/Workshop will move you using  the absolute path and typing cd Workshop will move you using the *relative path*
                - Relative paths have some shortcuts .. (two dots) means "the directory immediately above the one that I'm in" and . (one dot) means "the directory that I'm in"
- **ls** (list) prints out information about the directory that you’re in
- Working with Files:
- **mv** (move) moves a file from a source path to a destination path
- **cp** (copy) copies a file from a source path to a destination path
- **rm** (remove) removes a file or directory. **WARNING** there is no recycle bin on a  unix system. If you remove a file it is gone forever, with no chance of recovery. **DOUBLE CHECK YOUR COMMANDS**
- **zip/unzip** allows you to create and extract .zip archive files
- **tar** (tape archive) allows you to create and extract archive files that have the  extension of .tar, .tar.xz, .tar.gz, and more. (see here for more information: https://phoenixnap.com/kb/extract-tar-gz-files-linux-command-line)
- **less** allows you to view files on the command line without editing them
- **cat** (concatenate) prints out the contents of a file directly to your terminal
- **head/tail** show you the first and last 10 lines of a file, respectively. These are  useful for quickly making sure a file is in the right format
- **vim** is an advanced text editor that allows you to edit files for writing scripts. **nano** is a simpler text editor for writing scripts on the command line
- **chmod** changes the permission on files

 
















