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
mkdir -p /scratch/Users/$USER/workshop-day7/cellranger_count/eando
cd /scratch/Users/$USER/workshop-day7/cellranger_count
```

> **Make `eando` now, not later.** SLURM writes the `.out` and `.err` files the instant the job starts. If the directory named in your `--output`/`--error` paths does not already exist, the job fails immediately.

---

## 2. Build the batch script

Open a new file in vim:

```bash
vim cellrangerCount.sbatch
```

Add each block below in order. Read the short explanation before each one.

### 2a. Shebang and SLURM directives

```
#!/bin/bash
#SBATCH --job-name=cellrangerCount_male19
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=25gb
#SBATCH --time=18:00:00
#SBATCH --partition=short
#SBATCH --output=/scratch/Users/<your_username>/workshop-day7/cellranger_count/eando/cellrangercount.%j.out
#SBATCH --error=/scratch/Users/<your_username>/workshop-day7/cellranger_count/eando/cellrangercount.%j.err
```

Two things to customize here:

- **Replace `<github_username>`** in both the `--output` and `--error` paths with your own username.
- **Where did 8 CPUs and 25 GB come from?** They are not arbitrary — they were measured, and the obvious guess turns out to be very wrong. Use them for now; the reasoning is in the **appendix**, once your job is submitted and you have moved on.

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
  --localmem=24 \
  --localcores=<how_many_cores> # hint: number must match --ntasks above
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
#SBATCH --mem=25gb                                   # Memory limit
#SBATCH --output=/scratch/Users/<github_username>/workshop-day7/cellranger_count/eando/cellrangercount.%j.out
#SBATCH --error=/scratch/Users/<github_username>/workshop-day7/cellranger_count/eando/cellrangercount.%j.err

module load cellranger/7.2.0

echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

cellranger count --id=T21BM_male19 \
  --fastqs=/scratch/Shares/public/sread/data_files/day7a/fastq/sampled_fastq \
  --transcriptome=/scratch/Shares/public/sread/cookingShow/day7a/genomes/refdata-gex-GRCh38-2020-A \
  --sample=DSOX19_1 \
  --include-introns=true \
  --localmem=24 \
  --localcores=8      # must match --ntasks above

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

Later, when the run completes, `T21BM_male19/outs/` contains:

- **`web_summary.html`** — open this first, always. Estimated cells, reads per cell, **sequencing saturation**, fraction of reads in cells. It is the first thing a reviewer would look at.
- **`filtered_feature_bc_matrix/`** — `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`. This is what Seurat reads. *Filtered* means Cell Ranger has already decided which barcodes are cells and which are empty droplets — a decision made **for** you, by an algorithm, with a threshold. Which is exactly why the next worksheet does its own QC on top of it.
- **`raw_feature_bc_matrix/`** — every barcode, before that decision.

Full output reference: <https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-overview>

---
# Appendix — What did I actually use? And how do I choose `--ntasks` and `--mem` for Cell Ranger

When the job finishes, read what `/usr/bin/time -v` recorded:

```bash
grep -E "Maximum resident|Elapsed|Percent of CPU" eando/cellrangercount.*.err
```

- **`Maximum resident set size (kbytes)`** ÷ 1,048,576 = peak GB. Compare that to the 25 GB you requested.
- **`Percent of CPU this job got`** — 800% means you used all eight cores. 105% means you used one and wasted seven.

This ten-second habit is the difference between requesting resources by superstition and requesting them by measurement — and it is the only honest way to size your *next* job.

> On many clusters you would use `seff <jobid>` for this. **It does not work here**: this cluster has no SLURM accounting database (`sacct` reports *"Slurm accounting storage is disabled"*). Worth knowing — the tool everyone recommends is unavailable on a large fraction of real clusters, whereas `/usr/bin/time -v` works everywhere.

---

## The obvious guess was wrong

A previous version of this worksheet requested **`--mem=10gb`**. It ran fine. It ran fine every single time it was tested.

It would have brought the cluster down on the first day 25 people ran it at once.

## Step 1 — What does the cluster have?

Never trust `lscpu`: you are on the login node, and your job runs on a different machine.

```bash
sinfo -o "%P %.10l %.6D %C %m"
```

Here: **10 nodes on `short`, each with 8 CPUs and ~30 GB of RAM.**

## Step 2 — What does the job actually use?

This is the step everyone skips, and it is the only one that matters. Measuring the real run gives peak memory per stage:

```
14.8 GB   ALIGN_AND_COUNT          ← the peak
13.3 GB   DETECT_COUNT_CHEMISTRY
 2.6 GB   WRITE_GENE_INDEX
 1.3 GB   FILTER_BARCODES
 0.5 GB   BARCODE_CORRECTION
 ...
```

**Alignment peaks at nearly 15 GB** — half again more than the 10 GB that was being requested. That memory is mostly STAR loading the GRCh38 index into RAM, so it is a floor set by the *reference genome*, not something a smaller FASTQ shrinks. Subsampling the reads did not help at all.

## Step 3 — The arithmetic, for 25 people

Here is the thing to understand: **SLURM packs jobs onto nodes according to what you *request*, not what you *use*.** An under-request is not a harmless understatement. It is an instruction to the scheduler to overload the node.

| `--mem` | Jobs SLURM packs per node | RAM actually needed | Result |
| --- | --- | --- | --- |
| 10 GB | 3 | 3 × 14.8 = **44 GB** | Node has 30. **Kernel OOM.** |
| 15 GB | 2 | 2 × 14.8 = **29.6 GB** | Node has 30.4. **No margin.** |
| **25 GB** | **1** | **14.8 GB** | **Safe.** |

Hence **25 GB, and the whole node's 8 cores**. One job per node, ten nodes, so about a third of the class runs at a time.

## Why this is the real lesson

The 10 GB request passed every test it was ever given — because it was always tested by one person on an idle cluster. Nothing competed for RAM, so nothing broke.

It would have failed the first time 25 of us hit the alignment stage simultaneously. And it would not have failed cleanly: three students' jobs land on one node, together demand 44 GB of a 30 GB machine, and the **kernel's OOM killer** starts terminating processes semi-randomly. No tidy SLURM error. No message saying "memory." Just jobs dying for reasons nobody can reconstruct.

> **A resource request that works in testing can still be catastrophically wrong in production, and the only way to know is to measure.**
>
> And the corollary: **under-requesting memory is not modest. It is a lie to the scheduler, and other people pay for it.**

## How to measure your own

That is what `/usr/bin/time -v` was doing at the beggining . For a finished job:

```bash
grep -E "Maximum resident|Elapsed|Percent of CPU" eando/cellrangercount.*.err
```

Peak memory in KB ÷ 1,048,576 = GB. **Request that, plus ~30% headroom.** Then check `Percent of CPU` — if you asked for 8 cores and got 110%, you used one core and wasted seven.

Cell Ranger also records this per stage, which is where the table above came from:

```bash
for f in $(find T21BM_male19 -name "_jobinfo"); do
  python3 -c "
import json
d = json.load(open('$f'))
kb = d['rusage']['children']['ru_maxrss']
print(f\"{kb/1048576:6.2f} GB  {d['name'].split('.')[-3]}\")
" 2>/dev/null
done | sort -rn | head
```

## For your own data

10x recommend **at least 8 CPUs and 64 GB, preferably 16 and 128**, with negligible return past 32 cores or 256 GB. A real run (~10,000 cells at ~30,000 reads/cell) takes hours, not the ~30 minutes you saw here.

Do not copy our numbers. Copy the method: **find out what the cluster has → measure what your job needs → divide by how many people are competing → then verify against a real run.**

