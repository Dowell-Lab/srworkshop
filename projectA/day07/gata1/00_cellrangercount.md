# Day 7: Cell Ranger Count

**Author:** Chris Ozeroff · Short Read Workshop

## Overview

In this tutorial you will take a single-cell RNA-sequencing dataset and map it with the **Cell Ranger** pipeline. Cell Ranger performs alignment, filtering, and unique molecular identifier (UMI) and barcode counting. It outputs several files and directories, including a **count matrix** directory that we can analyze later with the R package **Seurat**.

> **Why a subsampled file?** Mapping a full scRNA-seq dataset takes a long time on a large dataset. To keep things moving, you will run Cell Ranger on a FASTQ that is ~10% the size of the real file from the paper. In the next worksheet you will be pointed to full-sized count matrices to analyze in Seurat.

> **What this worksheet is really for.** The count matrix you produce here is **not** used in the following worksheets — those read full-size matrices downloaded from GEO. What you are actually walking away with is a **working template**: a script you can point at your own FASTQs, on your own cluster, when you have your own data. Build it deliberately, and keep it.

## What you'll do

- Set up a working directory for the run
- Build a SLURM batch script from scratch, block by block
- Submit the job and confirm it is running

Rather than copying a finished script, you will assemble it yourself so you understand what each piece does. Paste each block below into your file as you go.

---

## 1. Log on and set up your directory

Log onto the AWS, pull the latest from GitHub:

```bash
cd ~/srworkshop
git pull
```

Make a directory for this run, and an `eando` (error-and-output) directory inside it for your log files:

```bash
mkdir -p ~/workshop-day7/cellranger_count/eando
cd ~/workshop-day7/cellranger_count
```

> **Make `eando` now, not later.** SLURM writes the `.out` and `.err` files the instant the job starts. If the directory named in your `--output`/`--error` paths does not already exist, the job fails immediately. That is why we create it before submitting rather than inside the script.

(If your clone of the repo lives somewhere other than `~/workshop-day7`, adjust these paths to match.)

---

## 2. Build the batch script

Open a new file in vim:

```bash
vim cellrangerCount.sbatch
```

Add each block below in order. Read the short explanation before each one.

### 2a. Shebang and SLURM directives

```bash
#!/bin/bash
#SBATCH --job-name=cellrangerCount_male19            # Job name
#SBATCH --nodes=1                                    # Number of nodes
#SBATCH --ntasks=1                                   # Number of CPUs (tasks)
#SBATCH --time=18:00:00                              # Time limit hrs:min:sec
#SBATCH --partition=short                            # Partition/queue on the server
#SBATCH --mem=10gb                                   # Memory limit
#SBATCH --output=/Users/<github_username>/workshop-day7/cellranger_count/eando/cellrangercount.%j.out
#SBATCH --error=/Users/<github_username>/workshop-day7/cellranger_count/eando/cellrangercount.%j.err
```

Two things to customize here:

- **Replace `<github_username>`** in both the `--output` and `--error` paths with your own username. Confirm your home path with `echo $HOME`.
- **`--ntasks=8`.** Eight is the number of cores we have available on the AWS. On your home university supercomputer you can request more. This has to match `--localcores` in the run command below.

> **Heads up:** the `~` shortcut does *not* reliably expand inside `#SBATCH` lines, which is why these paths are written out in full. Use `~` in the commands you type at the shell, but full absolute paths in the directives.

### 2b. Load Cell Ranger

```bash
module load cellranger/7.2.0
```

### 2c. Print job context

These lines write useful information to your `.out` log, which makes debugging much easier later.

```bash
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes
```

### 2d. The Cell Ranger count command

```bash
cellranger count --id=T21BM_male19 \
  --fastqs=/scratch/Shares/public/sread/data_files/day7a/fastq/sampled_fastq \
  --transcriptome=/scratch/Shares/public/sread/cookingShow/day7a/genomes/refdata-gex-GRCh38-2020-A \
  --sample=DSOX19_1 \
  --include-introns=true \
  --localcores=1      # must match --ntasks above
```

### 2e. Closing line

```bash
echo Job finished at `date +"%T %a %d %b %Y"`
```

Save and quit vim with `:wq`.

### What do all these flags in cellranger mean?

| Flag | What it does |
| --- | --- |
| `--id` | Names the output directory. Here Cell Ranger writes everything into `T21BM_male19/`. |
| `--fastqs` | Path to the directory holding the FASTQ files. |
| `--transcriptome` | Path to the reference transcriptome (GRCh38). |
| `--sample` | The sample prefix that matches the FASTQ filenames. Here that prefix is `DSOX19_1`. |
| `--include-introns` | Counts reads mapping to introns as well as exons. |
| `--localcores` | Number of CPUs Cell Ranger may use. Must equal `--ntasks`. |

---

## 3. Check your script before submitting

Run through this quick checklist:

- [ ] `<github_username>` is replaced in **both** the `--output` and `--error` paths
- [ ] `--ntasks` and `--localcores` are both `8` and match each other
- [ ] The `eando` directory exists (`ls ~/workshop-day7/cellranger_count/eando`)

<details>
<summary>Stuck? Click to see the complete script.</summary>

```bash
#!/bin/bash
#SBATCH --job-name=cellrangerCount_male19            # Job name
#SBATCH --nodes=1                                    # Number of nodes
#SBATCH --ntasks=8                                   # Number of CPUs (tasks)
#SBATCH --time=18:00:00                              # Time limit hrs:min:sec
#SBATCH --partition=short                            # Partition/queue on the server
#SBATCH --mem=10gb                                   # Memory limit
#SBATCH --output=/Users/<github_username>/workshop-day7/cellranger_count/eando/cellrangercount.%j.out
#SBATCH --error=/Users/<github_username>/workshop-day7/cellranger_count/eando/cellrangercount.%j.err

module load cellranger/7.2.0

echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

cellranger count --id=T21BM_male19 \
  --fastqs=/scratch/Shares/public/sread/data_files/day7a/fastq/sampled_fastq \
  --transcriptome=/scratch/Shares/public/sread/cookingShow/day7a/genomes/refdata-gex-GRCh38-2020-A \
  --sample=DSOX19_1 \
  --include-introns=true \
  --localcores=1      # must match --ntasks above

echo Job finished at `date +"%T %a %d %b %Y"`
```

</details>

---

## 4. Submit the job

```bash
sbatch cellrangerCount.sbatch
```

---

## 5. Confirm it is running

Check the queue for your job:

```bash
squeue -u $USER
```

Give it a minute, then peek at the error log to confirm nothing has gone wrong:

```bash
cat eando/cellrangercount.*.err
```

If your job is still listed in the queue and the `.err` file has no error messages, you are good to move on.

---

## 6. Move on to the next worksheet

You do not need to wait for the job to finish. Continue to the next worksheet while Cell Ranger runs.

> **What you'll have when it finishes (~2 hrs):** an output directory named `T21BM_male19/`. Inside `T21BM_male19/outs/` is a `filtered_feature_bc_matrix/` directory containing the `matrix`, `barcodes`, and `features` `.tsv` files that Seurat reads.
>
> Overview of all Cell Ranger outputs:
> https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-overview
