# Project B Day 8 | ChIP-seq: Analyzing TF ChIP-seq data

Today we process ChIP-seq data. We will start with preprocessing the raw data and proceed to peak calling, motif identification and annotation.
We will be analyzing TP53 ChIP-seq data done in HCT116 cells from Andrysik et al. 2017 (doi:10.1101/gr.220533.117).

## In-class worksheets

1. Check read quality and preprocess raw fastq files

   - Check that quality of reads with `fastqc` 
   - Trim reads using `trimmomatic`

2. Map reads to reference genome with `HISAT2`

3. Assess the library complexity based on the read mapping with `preseq` and generate a QC report summary with `multiqc`

4. Call ChIP-seq peaks with `MACS2`

5. TF motif discovery with `MEME`

6. Compared TF motif against TF databases with `Tomtom`

## Homework

Repeat the inclass worksheets but in MCF7 samples.
