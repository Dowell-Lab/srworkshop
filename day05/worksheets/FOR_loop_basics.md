# <code>for</code> loop basics
Author: Lynn Sanford, 2025

## This worksheet covers:
<ul>
<li><code>for</code> loop syntax in Unix</li>
<li>Serial <code>for</code> loops (most basic)</li>
<li>How to use <code>for</code> loops with SLURM to parallelize scripts</li>
</ul>

## <code>for</code> loop syntax in Unix
<code>for</code> loops in any coding language are blocks of code that loop for a given condition. In this class, our most common use of <code>for</code> loops is to loop through code for each file of a certain type in a directory (such as all FASTQ files).

A simple example in Unix looks like this:
```
for filepath in ./*.fastq
do

# The actual operations of the loop
echo "$filepath"

done
```
In this loop, Unix takes in <code>./*.fastq</code> as a list of paths to fastq files in the current directory (I'm using the directory to the fastq files we're using today, <code>/scratch/Shares/public/sread2025/data_files/day5/fastq/for_loops_fastq</code>).

Then, in each iteration of the loop, the next path in that list is stored as the variable <code>filepath</code>. <code>do</code> and <code>done</code> then bracket the code block that runs during each iteration.

## Serial <code>for</code> loops
A <code>for</code> loop will loop its code serially, so it finishes the whole loop of code and then begins the next iteration. When the input list has no more elements, the loop finishes.

If you run a for loop on the head node, it will run serially on the head node. For example, the following loop will pause for 5 seconds per file.

```
for filepath in ./*.fastq
do
echo $filepath
sleep 5
done
```

Running this loop on the head node will use head node resources for the entire length of the loop (in this case, 4 files = 20 seconds).

You can also run this loop in a SBATCH script (Day5 <code>for_loop_sleep_serial.sbatch</code>). This will run the loop serially on a child node.

## How to use <code>for</code> loops with SLURM to parallelize scripts

What if your code block takes an hour to run on each file? What if it takes two days to run on each file? What if you have forty files?

In order to speed up analysis, we can parallelize running scripts. This involves using a <code>for</code> loop to read in all the input files, as in the previous examples, and then submit a separate SBATCH script for each file.

>**Note:** This involves running a <code>for</code> loop on the head node. But submitting jobs is not at all resource-intensive, so this is an acceptable use of the head node.

To do this, we use two scripts
<ul>
<li>The SBATCH script that runs on each file (<code>for_loop_sleep_parallel.sbatch</code>)</li>
<li>A coordinator BASH script (<code>for_loop_sleep_coordinator.sh</code>)</li>
</ul>

The coordinator script takes each value of the input list and feeds it into the SBATCH submission with the parameter <code>--export</code>. Each submission then goes as a separate job to a child node.