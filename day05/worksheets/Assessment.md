# Part 3: Assessment!
This is NOT an actual test, but a way for you to *assess* what you've learned and if practicing over the weekend would be wise given we'll be moving forward next week. Keep track of how long it takes you to complete the full pipeline for processing FASTQ files. 

You can choose either ONE of the fastq files in /scratch/Shares/public/sread2023/day5/assessment_fastq OR do ALL of them with a loop. Either way, adapt the scripts you made in your first week and run the following steps on your dataset:
1.	Run a QC check on the raw FASTQ files
2.	Trim the FASTQ files and run a QC check again
3.	Map the FASTQ files to hg38
4.	Compress the mapped files into BAM files
5.	Generate bedGraph files from BAM files
6.	Generate TDF files from bedGraphs 
7.	Transfer the TDF to your local computer and visualize using the IGV Web App

**NOTE**: You are ALLOWED and encouraged to use your copy of your previous scripts with minor edits of names, etc. When doing pipelines, you will often be given code that you need to adapt rather than have to reinvent the wheel.

**BONUS**: Combine the scripts into a single pipeline that takes the FASTQ files all the way to TDFs. 
