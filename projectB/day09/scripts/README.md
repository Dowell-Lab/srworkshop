# Short Read Day 9: Peak-centric scripts
- Author: Jessica Huynh-Westfall

The scripts for today's module will be focusing on peak-centric sequencing data. Briefly, the scripts in the directory will help you evaluate your sequencing data, quality control using input. Using MACS to call peaks in your data, and exploring downstream analysis such as bedtools intersect to consider removing blacklist regions, motif discovery and visualizing your data.

1. d9_preprocess_chipseq.sbatch - Run fastQC on your sequencing data and alignment of data to reference genome
2. d9_qc_chipseq.sbatch - QC including Preseq, RSeQC, and multiQC 
3. d9_macs.sbatch - Peak calling and some downstream applications
4. d9_bedtools.sbatch - Bedtools
