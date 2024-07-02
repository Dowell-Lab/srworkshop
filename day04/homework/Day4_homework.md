# Day 4 Homework 
Authors: Jessica Westfall & Rutendo Sigauke\
Edited: Lynn Sanford, 2024

## Introduction
In a future day of the workshop we will go into more details about RNA-seq libraries. This homework will go over the tasks that we did in class and provide more practice. We will be using these files on subsequent days of the workshop.

## Homework
<ol>
  <h3><li>Open the homework script</h3>
    Copy the <code>example_process_rnaseq.sbatch</code> script from the GitHub repo day04 <code>scripts</code> folder to your scratch day4 working directory and make the necessary edits to do the following tasks. 
  </li>
  <h3><li>FastQC</h3>
    Evaluate the fastq files in
    <code>/scratch/Shares/public/sread2024/homework_data_files/day4/</code>.
    <br />
    <br />
    How is the quality of these sequence libraries? Things we want to look at are:
    <ul>
      <li>GC Content (Is the library contaminated?)</li>
      <li>Adaptor Content (Did the sequencer read into our adaptors?)</li>
      <li>Read Duplication (Is our sample overamplified? Depends on library type…)</li>
      <li>Sequence Quality/N content (How confident was the sequencer in calling each base?)</li>
      <li>Sequence Quality based on flow cell location (Was there a sequencing failure?)</li>
      <li>Base Identity at each location (Was there any bias in amplification/ligation?)</li>
    </ul>
  </li>
  <h3><li>Trimmomatic</h3>
    Trim the sequence library to remove adapters. Save the output to have the suffix <code>_trim.fastq</code> to track the trimming.
    <ul>
      <li>
        When trimming, consider if you are trimming single-ended (SE) or paired-ended (PE) reads. Consider other parameters listed for trimming:
        <ul>
          <li><code>ILLUMINACLIP</code>: Cut adapter and other illumina-specific sequences from the read</li>
          <li><code>SLIDINGWINDOW</code>: Performs a sliding window trimming approach. It starts scanning at the 5’ end and clips the read once the average quality within the window falls below a threshold</li>
          <li><code>LEADING</code>: Cut bases off the start of a read, if below a threshold quality</li>
          <li><code>TRAILING</code>: Cut bases off the end of a read, if below a threshold quality</li>
          <li><code>CROP</code>: Cut the read to a specified length by removing bases from the end</li>
          <li><code>HEADCROP</code>: Cut the specified number of bases from the start of the read</li>
          <li><code>MINLEN</code>: Drop the read if it is below a specified length</li>
        </ul>
      </li>
    </ul>
  </li>
  <h3><li>FastQC of trimmed files</h3>
    Run FastQC on the trimmed fastq and reevaluate the fastq.  
    <br /><br />
    Take a look at the fastQC of the trimmed fastq and ask yourself, are these files trimmed well? Can you adjust the parameters to make the trimming more stringent to remove adapter content? 
  </li>
  <h3><li>HISAT2</h3>
    Edit the sbatch script to map the two corresponding paired-end fastq files. For example, <code>chr21Ethan_repA.RNA.end1.fastq</code> and <code>chr21Ethan_repA.RNA.end2.fastq</code>
    <br /><br />
    Mapping has different parameters to change the mapping efficiency. What would happen if you alter the script which currently has <code>--very-fast</code> to <code>--very-sensitive</code>?  
  </li>
  <h3><li>SAMTOOLS</h3>
    Now we have a huge SAM file.  
    <ul>
      <li>
        Convert the SAM file to a BAM file to produce a compressed binary file that takes up less space.
        <ul><li>  
          What is the size difference between SAM versus BAM? What is the difference between the two filetypes that contributes to the size difference?
        </li></ul>
      </li>
      <li>
        Sort and index your BAM files.
      </li>
      <li>
        Remove the unnecessary large files.
      </li>
    </ul>
  </li>
  <h3><li>IGV</h3>
    <ul>
      <li>Transfer your BAM files to your local computer.</li>
      <li>
        Open IGV and visualize the sorted BAM files.
        <ul>
          <li>Do you know which reference genome your sequence reads were aligned and mapped to?</li> <li>Since both hg19 and hg38 are different versions of the human genome, are they interchangeable in IGV?</li>
        </ul>
      </li>
    </ul>
  </li>
</ol>