# Project B Day 8 | ChIP-seq: Analyzing TF ChIP-seq data

Today we process ChIP-seq data. We will start with preprocessing the raw data and proceed to peak calling, motif identification and annotation.
We will be analyzing TP53 ChIP-seq data done in HCT116 cells from Andrysik et al. 2017 (doi:10.1101/gr.220533.117).

Metadata for the samples can be found on NCBI here [SRP083188](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP083188&o=acc_s%3Aa).

| Run (SRR)         | Cell line  | Sample Type     | 
| :---------------- | :-------:  | :-------------: |
| SRR4090089        |  HCT116    | Input           |
| SRR4090090        |  HCT116    | DMSO treated    |
| SRR4090091        |  HCT116    | Nutlin treated  |
| SRR4090092        |  MCF7      | Input           |
| SRR4090093        |  MCF7      | DMSO treated    |
| SRR4090094        |  MCF7      | Nutlin treated  |
| SRR4090095        |  SJSA      | Input           |
| SRR4090096        |  SJSA      | DMSO treated    |
| SRR4090097        |  SJSA      | Nutlin treated  |

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
