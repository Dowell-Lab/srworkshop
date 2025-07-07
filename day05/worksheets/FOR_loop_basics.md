# `for` loop basics
Author: Lynn Sanford, 2025

## This worksheet covers:
<ul>
<li>`for` loop syntax in Unix</li>
<li>Serial `for` loops (most basic)</li>
<li>How to use `for` loops with SLURM to parallelize scripts</li>
</ul>

## `for` loop syntax in Unix
`for` loops in any coding language are blocks of code that loop for a given condition. In this class, our most common use of `for` loops is to loop through code for each file of a certain type in a directory (such as all FASTQ files).

A simple example in Unix looks like this:
```
for filepath in ./*.fastq
do

# The actual operations of the loop
echo "$filepath"

done
```
In this loop, Unix takes in `./*.fastq` as a list of paths to fastq files in the current directory (For the purposes of this demonstration, I've moved into the directory to the fastq files we're using today, `/scratch/Shares/public/sread2025/data_files/day5/fastq/for_loops_fastq`).

Then, in each iteration of the loop, the next path in that list is stored as the variable `filepath`. `do` and `done` bracket the code block that runs during each iteration.

## Serial `for` loops
A `for` loop will loop its code serially, so it finishes the whole loop of code and then begins the next iteration. When the input list has no more elements, the loop finishes.

If you run a `for` loop on the head node, it will run serially on the head node. For example, the following loop will pause for 5 seconds per file.

```
for filepath in ./*.fastq
do
echo $filepath
sleep 5
done
```

Running this loop on the head node will use head node resources for the entire length of the loop (in this case, 4 files = 20 seconds). As you now know, this is NOT a good idea.

You can also run a loop in a SBATCH script (Day5 `for_loop_sleep_serial.sbatch`). This will run the loop serially on a child node.

## How to use `for` loops with SLURM to parallelize scripts

What if your code block takes an hour to run on each file? What if it takes two days to run on each file? What if you have forty files?

In order to speed up analysis, we can parallelize running scripts. This involves using a `for` loop to read in all the input files, as in the previous examples, and then submit a separate SBATCH script for each file.

>**Note:** You'll notice that this actually does involve running a `for` loop on the head node. But submitting jobs is not at all resource-intensive and is very fast, so this is an acceptable use of the head node.

To do this, we use two scripts:
<ul>
<li>The SBATCH script that runs on each file (`for_loop_sleep_parallel.sbatch`)</li>
<li>A coordinator BASH script (`for_loop_sleep_coordinator.sh`)</li>
</ul>

Look at the coordinator script. It takes each value of the input list as in the previous examples, but now it feeds that path into an SBATCH submission with the parameter `--export`. Each submission then goes as a separate job to a child node.

Now look at the SBATCH script. It does NOT contain a loop, since it's only working on one file. But it does use the variable exported by the coordinator script.

When running the coordinator script, you will immediately see the jobs submitted.