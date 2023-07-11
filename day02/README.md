# Day 2! Working in the Unix Environment

## Part 1: Copy a bedGraph file to our home directory and inspect it
Author: Zach Maas

*Goal*: Help you familiarize yourself with the unix computing environment by performing some basic tasks and troubleshooting things on your own.

### What you get to do:
1. Log onto the Short Read Workshop server using SSH
2. Move to your user directory on */scratch/Users/<username>*
3. Make a folder called *workshop-day2* and a folder inside of that called *bedfiles*
4. Make folders inside *workshop-day2* called *results*, *scripts*, *bin*, and *data*
5. Copy the folder in */scratch/Shares/public/sread2022/data_files/day2/bedfiles/chr1_bedfiles.tar.gz* to your *bedfiles* directory
6. Decompress and extract the *chr1_bedfiles.tar.gz* folder in your *bedfiles* directory
7. Look at the format of the extracted files. What kind of information do these files give you?

### Some help/encouragement:
If you don’t know how to do something, don’t be discouraged! Look at the documentation  we’ve provided below or search the internet for how to do what you’re trying to accomplish. 

**Some useful commands and tools:**
- Basic Utilities:
	- **man** (manual) followed by the name of a command will tell you how to use a  command 
	- **ssh** lets you log into a remote server like the one we’re using during the short  read workshop 
	- **whoami** will print your current active user
	- **uname-n** This is useful to  make sure you’re actually on the short read server and not your local machine
- Navigating Around:
	- **pwd**  (print working directory) will tell you which directory you’re in 
	- **cd** (change directory) lets you move to a different directory. You can move folders using absolute or relative paths. For example, if  you’re starting in /scratch/ and you want to be in the Workshop directory inside of it, typing cd /scratch/Workshop will move you using  the absolute path and typing cd Workshop will move you using the *relative path* 
		- Reelative paths have some shortcuts .. (two dots) means "the directory immediately above the one that I'm in" and . (one dot) means "the directory that I'm in"
	- **ls** (list) prints out information about the directory that you’re in 
- Working with Files:
	- **mv** (move) moves a file from a source path to a destination path 
	- **cp** (copy) copies a file from a source path to a destination path 
	- **rm** (remove) removes a file or directory. **WARNING** there is no recycle bin on a  unix system. If you remove a file it is gone forever, with no chance of recovery. **DOUBLE CHECK YOUR COMMANDS**
	- **zip/unzip** allows you to create and extract .zip archive files 
	- **tar** (tape archive) allows you to create and extract archive files that have the  extension of .tar, .tar.xz, .tar.gz, and more.  
	- **less** allows you to view files on the command line without editing them
	- **cat** (concatenate) prints out the contents of a file directly to your terminal
	- **head/tail** show you the first and last 10 lines of a file, respectively. These are  useful for quickly making sure a file is in the right format 
	- **vim** is an advanced text editor that allows you to edit files for writing scripts. **nano** is a simpler text editor for writing scripts on the command line 
