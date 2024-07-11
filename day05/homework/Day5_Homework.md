# Homework Day 5
Authors: Daniel Ramírez (2022), Samuel Hunter (2023), Hope Townsend (2024)

## Everyone
1.	Watch Day 6 videos
2.	Install R (https://rweb.crmda.ku.edu/cran/)) and R-studio (https://www.rstudio.com/products/rstudio/download/)

## If you took more than two hours to do the Assessment:
0.	Create a directory to hold your bash script

1.	Create the startings of a bash script (Hint: use Vim and add the one-line heading needed to tell your computer that it is a bash script)

2.	In the bash script print out the words `Hello world!`. Run the script with `sh` to check it worked!

3.	Now let’s practice using variables. Let’s say we want the script to say `Hello` to any person we want. Using bash variables, have the script print `Hello Ethan!` and `Hello Eric!` by using a bash variable called `NAME`.
    - Hint: define `NAME` as one name, then use `echo` to print out a string including that variable. Run the script, then change the value of `NAME` to the other name and run it again.

4.	Now edit this script to become a SBATCH script 
    - Add the SBATCH required headings to the bash script to make it SBATCH. Use `nodes=1`, `ntasks=1`, `time` is 1 minute, and `mem=1gb`. Don’t forget to change the paths to the output and error!

5.	Now annotate each step of `day05/scripts/d5-fastq-to-tdf.sbatch` with what the different parameter options mean for the following commands: (Examples of what a good annotation might look like is on lines 192-194)
    - `hisat2` (line 52)
    - `samtools view` (line 63 & 69)
    - `samtools index` (line 77)
    - `genomeCoverageBed` (line 153)

6.	If you’re still struggling with Bash Variables (like `${}`), go through `day05/scripts/d5-fastq-to-tdf.sbatch` and write out what each of the bash variables equals on each line.

7.	What does `|` do? Find an example of where this is used.

8. If you want to challenge yourself:
    - How do you have to edit the script `d5-fastq-to-tdf.sbatch` and `run_d5-fastq-to-tdf.sh` for the scripts to run on datasets sequenced as single-end reads? 

