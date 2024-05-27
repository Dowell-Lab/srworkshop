# Part 1 - Bedgraphs and TDFs
Authors: Daniel Ramírez (2022), Samuel Hunter (2023), Hope Townsend (2024)
## 0. Staging your work area
1. Navigate to **your** github repo clone and **git pull** to "pull" any updates that someone had "pushed" to the repository to your own computer. 
```
cd /Users/your_username/sr2024
git pull
```
2. Make a directory for error and output files for today (usually we call the directory e_and_o OR eofiles), if you don't already have one. (Hint: if you forget, google how to do this on the command line)

## 1. Making Bedgraphs and TDFs
BAM files that we made yesterday are an intermediate file (mapped reads to the genome). We want to extrapolate more info fsrom them, in this case, the number of reads that are at each region. This will allow downstream analysis like determining how transcribed a region is or for easier visualization.

* Difference between bedgraphs and TDFs: TDFs are the compressed version of bedgraphs (meaning they take up less space on your computer and are faster to work with).
* If you didn't make it through all the scripts yesterday, we've provided the files you need at the following folder:
  
```
/scratch/Shares/public/sread2023/cookingShow/day5/bam/
	chr21Eric_repA.RNA.sorted.bam
	chr21Eric_repA.RNA.sorted.bam.bai
```

1. First, we need to make a script that will run the code to get the bedGraphs and TDFs from the bam files.
	a. Navigate to your github repository under day05/scripts (*Hint* cd). Then open the script/file **d5-bam-to-tdf.sbatch** (*Hint* vim or vi). 
	
 	b. Edit the script to match YOUR information (Remember to delete the <> brackets!):

	*  JOB_NAME to the name bam_to_tdf
	* output and error files (--output and --error) should use your directory made in Step 0. **Important**: keep the %x_%j.out and %x_%j.err parts as those are the file names where the % symbol means that the job id (j) will populate the names. This makes it really easy to know which error/output files connect to which jobs you ran.
	* OUTDIR to point to your day5 directory (change <YOUR_USERNAME> to your username)

3. Second, before you run a script, you should make sure you have a basic understanding of what is happening. A good way of doing this is by adding comments to the script with "#". You see examples of this in the first X.
	
 	a. Read through the comments that already exist to understand what is happening. Fill out the parts that say "<ANNOTATE ...>" with your own annotations according to the prompt. If you have any key questions, raise your red sticky note. UNDERSTANDING how the variables are being used and the basic steps is key for the rest of the assessment.
	
 	b. Now you can exit the script (make sure to save the file when you do). Then run it using the code `sbatch d5-bam-to-tdf.sbatch`.
    
5. Check that everything ran smoothly!

	a. Look at the error and output files you generated and see if everything ran correctly.

	b. Download the TDF to your local computer and open it up inthe IGV web app (or Desktop version if you have it installed).

    
# Part 2 BONUS: For Loops and Pipelines
**If you finish Part 1 before we start the assessment**, start working here. Otherwise, go straight to the assessment (Assessment.md).

In the past couple days, we've run samples one-by-one using individual scripts for each step. This is fine but becomes tedious when we're using dozens/hundreds of samples. Instead we can use a second script to submit jobs automatically using a for loop:

1. Open the script loop_example.sh and read it. What do you expect to be printed out?
2. Now exit vim and run the script with `bash example_for_loop.sh`
	* How many times did we go through the loop?
	* How many times did DONE! get printed out? Why?
	* What would happen if you changed the 3 to 5?
3. Now let's see how we might apply this to our code. In our code, we want to do the same process for multiple files. We can do this by pointing to the folder with the files we want to keep running. 

	a. So open loop_example_with_files.sh and edit the path to point to your fastq files directory. Run it same as last time.

	b. Now instead of just printing th efile name, let's add the command `head -n 1 ${FOLDER}/${index}` under the line do. In this case, we should also be seeing the first line of the files printed out.
4. But we are ultimately trying to run many samples at once on the same pipeline. We could do a for loop within the script with the code, but this means that a sample couldn't be run on the code till a previous sample is finished. So how can we do it so that all the samples can run at the SAME TIME. 

	* Think of it as we have an obstacle course that 5 kids need to run on. We could have each kid take turns, OR we could just have 5 obstacle courses and have the kids run at the same time. Computers make the latter much easier.

	a. First, we need the code that the samples will be run on (this doesn't change between samples and is like the obstacle course). In our example, this would be the **d5-fastq-to-tdf.sbatch**. This code runs the same code as done in d5-bam-to-tdf.sbatch but starts with the fastq files instead. Briefly edit the following portions to allow it to apply to you (Remember to delete the <> brackets!):
	* YOUREMAIL to your email
 	* output and error files (--output and --error) should use your directory made in Step 0. Important: keep the %x_%j.out and %x_%j.err parts as those are the file names where the % symbol means that the job id (j) will populate the names. This makes it really easy to know which error/output files connect to which jobs you ran.

	b. Second, we need our samples/input-files that need to go through the code (like the kids who need to run the obstacle courses). In this case, this would be the fastqs in your fastq directory.

	c. Finally, we need to have a script that will run tell the samples to be run through the code (think of it as the person who says Ready Set Go). Naming this "run_script_using.sbatch" usually helps keep track of it. An example can be found with **run_d5-fastq-to-tdf.sh**. 
	* This is where we actually use the for loop
4. Now edit the script run_d5-fastq-to-tdf.sh to run. 

	a. Edit the OUTDIR variable to point to wherever you want (I'd recommend your day5).

	b. Read the annotations and make sure you understand what's happening. Check out the variables to which you're assigning values in the script to see what you're doing. If you don't understand what's happening, hold up a red sticky note.

	c. Now run the script with `bash run_d5-fastq-to-tdf.sh`

   	* You should see output like this:
	    ```
	    sample1_day5_igv.RNA
	    Submitted batch job (number here starting with 9)
	    sample2_day5_igv.RNA
	    Submitted batch job (number here starting with 9)
	    ```
    
Now you were able to get the TDF files with a single script!

