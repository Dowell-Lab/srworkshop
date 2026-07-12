# Day 7: Cell Ranger Count

**Author:** Chris Ozeroff · Short Read Workshop

## Overview

Cell Ranger takes raw single-cell FASTQ files and turns them into a **count matrix** — a table of genes × cells that Seurat can read. Along the way it demultiplexes cell barcodes, aligns reads, decides which barcodes are real cells, and counts UMIs.

In this worksheet you will build a SLURM batch script for `cellranger count` from scratch, block by block, and submit it.

> **What this worksheet is really for.** The count matrix you produce here is **not** used in the following worksheets — those read full-size matrices downloaded from GEO. What you are actually walking away with is a **working template**: a script you can point at your own FASTQs, on your own cluster, when you have your own data. Build it deliberately, and keep it.

> **Why a subsampled FASTQ?** A real Cell Ranger takes a long time to run. Yours runs on a FASTQ subsampled to ~10% so that it finishes in a reasonable time and 25 of us can run at once.

## What you'll do

- Pull today's materials from GitHub (in **home**)
- Set up a working directory on **scratch**, where the job actually runs
- Build the batch script block by block
- Submit it and confirm it started

---

## 1. Log on and pull today's materials

The course repo lives in your **home** directory. This is how you get each day's worksheets and scripts:

```bash
cd ~/srworkshop
git pull
```

> **The split, for every day of this course:**
>
> | | **Home** (`~/srworkshop`) | **Scratch** (`/scratch/Users/$USER/`) |
> |---|---|---|
> | Holds | The git repo: worksheets, template scripts | Everything you run, and everything a job writes |
> | Why | Version-controlled, backed up, small | Big, fast, disposable — built for heavy job I/O |
>
> **You `git pull` in home. You run jobs on scratch.** Cell Ranger writes a BAM file and a deep tree of intermediates — easily several GB even from our subsampled FASTQ. Home directories are quota-limited and slow, and heavy job I/O is exactly what scratch exists for.

---

## 2. Set up your working directory on scratch

Make a directory for this run, and an `eando` (**e**rror **and** **o**utput) directory inside it for your SLURM logs:

```bash
mkdir -p /scratch/Users/$USER/workshop-day7/cellranger_count/eando
cd /scratch/Users/$USER/workshop-day7/cellranger_count
pwd
```

Note what `pwd` prints — that is the path you will paste into the script in a moment.

> **Make `eando` now, not later.** SLURM opens the `.out` and `.err` files the *instant* the job starts. If the directory named in your `--output`/`--error` paths does not already exist, the job dies immediately — and because it died before your script ran, there is no log telling you why. Silent failures are the expensive ones.

---

## 3. Build the batch script

You are now inside `/scratch/Users/<your_username>/workshop-day7/cellranger_count`. Open a new file:

```bash
vim cellrangerCount.sbatch
```

Add each block below in order. Read the short explanation before each one.

### 3a. Shebang and SLURM directives

```bash
#!/bin/bash
#SBATCH --job-name=cellrangerCount_male19
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=10gb
#SBATCH --time=18:00:00
#SBATCH --partition=short
#SBATCH --chdir=/scratch/Users/<your_username>/workshop-day7/cellranger_count
#SBATCH --output=/scratch/Users/<your_username>/workshop-day7/cellranger_count/eando/cellrangercount.%j.out
#SBATCH --error=/scratch/Users/<your_username>/workshop-day7/cellranger_count/eando/cellrangercount.%j.err
```

**Replace `<your_username>`** in all three paths. It is whatever `pwd` printed above; `echo $USER` will also tell you.

| Directive | What it does |
| --- | --- |
| `--ntasks=2` | CPUs. Must match `--localcores` below. |
| `--mem=10gb` | Memory. Must match `--localmem` below. |
| `--time=18:00:00` | Walltime ceiling. Ask for too little and SLURM kills you mid-run. Overestimate while learning. |
| `--chdir` | The directory the job runs in. |
| `--output` / `--error` | Log files. `%j` expands to the job ID, so each submission gets its own log instead of overwriting the last. |

> **`--chdir` is the one people skip, and it causes real confusion.** Cell Ranger writes its output directory wherever the job happens to be running — which by default is *wherever you were standing when you typed `sbatch`*. Your `#SBATCH` log paths are absolute and go where you told them; the job's working directory is not, and nothing warns you they differ. Setting `--chdir` explicitly means your output lands on scratch no matter where you submit from.

> **`~` and `$USER` do not expand in `#SBATCH` lines.** SLURM's parser reads those lines, not bash, so it has no idea what they mean. Use `~` and `$USER` freely at the command prompt, but write **full, literal, absolute paths** in the directives.

### 3b. Load Cell Ranger

```bash
module load cellranger/7.2.0
```

> Pin the **version**, not just `cellranger`. Six months from now, reproducing your counts means knowing exactly which version produced them.

### 3c. Print job context

```bash
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes
echo Working directory: `pwd`
```

These look like filler. They are not — they land in your `.out` log and tell future-you which node, which version, and which job produced a given result. Do this in every script you write.

### 3d. The Cell Ranger command

```bash
/usr/bin/time -v cellranger count \
  --id=T21BM_male19 \
  --fastqs=/scratch/Shares/public/sread/data_files/day7a/fastq/sampled_fastq \
  --transcriptome=/scratch/Shares/public/sread/cookingShow/day7a/genomes/refdata-gex-GRCh38-2020-A \
  --sample=DSOX19_1 \
  --include-introns=true \
  --localcores=2 \
  --localmem=10
```

### 3e. Closing line

```bash
echo Job finished at `date +"%T %a %d %b %Y"`
```

Save and quit vim with `:wq`.

---

## 4. What every flag does

| Flag | What it does |
| --- | --- |
| `--id` | Names the output directory. Here everything lands in `T21BM_male19/`. |
| `--fastqs` | Path to the **directory** holding the FASTQs — not to a file. Cell Ranger searches it. |
| `--transcriptome` | Path to a prebuilt 10x reference (GRCh38). 10x ships these; you don't build one unless you work in a non-model organism (`cellranger mkref`, and that's its own afternoon). |
| `--sample` | The filename **prefix**. Here `DSOX19_1`, matching files like `DSOX19_1_S1_L001_R1_001.fastq.gz`. Get this wrong and Cell Ranger reports "no FASTQs found" while staring straight at them — the single most common Cell Ranger error. |
| `--include-introns` | Counts intronic reads too, capturing pre-mRNA. Matters for nuclei, and for progenitors with lots of nascent transcription. It changes your numbers, so report it in your methods. |
| `--localcores` | CPUs Cell Ranger may use. **Must equal `--ntasks`.** |
| `--localmem` | Memory in GB Cell Ranger may use. **Must equal `--mem`.** |

> **Why `--localcores` and `--localmem` are not optional.** Cell Ranger never asks SLURM what it was allocated. By default it inspects the whole machine and helps itself to **all available cores and 90% of detected memory**. SLURM then kills the job for exceeding its allocation, and the error message will not obviously say "memory." 10x say the same in their own docs: always set both.

> **Two bash gotchas in that command block:**
>
> - The `\` escapes the newline, letting one long command span several lines. **A space *after* the backslash breaks it** — bash treats it as the end of the command, and everything below becomes garbage. If you get a bizarre error, check for a trailing space.
> - **Do not put `# comments` on continuation lines.** It happens to work on the *last* line and silently breaks the command anywhere else. Put comments above the block.

> **What is `/usr/bin/time -v` doing there?** It records what the job actually used and prints it at the bottom of your `.err` log — most usefully `Maximum resident set size (kbytes)`, the peak memory. It costs nothing, and it is how you find out whether your `--mem` request was sensible. See §8.

---

## 5. How many cores and how much memory?

Short answer for today: **2 CPUs and 10 GB.** Here is where that comes from, because one day you will have to make this call yourself.

Ask SLURM what exists. Do **not** trust `lscpu` — you are on the login node, and your job runs somewhere else entirely:

```bash
sinfo -o "%P %.10l %.6D %C %m"
```

On this cluster that reports roughly **8 CPUs and 30 GB per node**, across 10 nodes on `short`.

Now the arithmetic. There are ~25 of us. At 10 GB each, a 30 GB node fits **3 jobs** — capping us at 30 concurrent slots, which is enough for everyone. And since *memory* has already limited us to 3 jobs per node, each of those jobs can take 2 of the node's 8 CPUs for free. Hence **2**.

Notice what that means: **memory is the binding constraint here, not CPUs.** Asking for more RAM than you need directly costs your classmates a job slot.

> **These are not realistic numbers.** For real data, 10x recommend **at least 8 CPUs and 64 GB, preferably 16 and 128** — and even then a normal run (~10,000 cells at ~30,000 reads/cell) takes *hours*. We use 2 cores and 10 GB on a 10%-subsampled FASTQ so that 25 jobs fit on a small teaching cluster.
>
> **Do not copy these numbers into your own work. Copy the reasoning:** find out what the cluster has, find out what the software needs, divide by how many people are competing for it — then *measure what your job actually used* and size the next one from that.

---

## 6. Check your script, then submit

- [ ] `<your_username>` replaced in **all three** paths (`--chdir`, `--output`, `--error`)
- [ ] Those paths point to **scratch**, not home
- [ ] `--ntasks=2` and `--localcores=2` agree
- [ ] `--mem=10gb` and `--localmem=10` agree
- [ ] `eando/` exists
- [ ] No trailing spaces after any `\`

<details>
<summary>Stuck? Click to see the complete script.</summary>

```bash
#!/bin/bash
#SBATCH --job-name=cellrangerCount_male19
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=10gb
#SBATCH --time=18:00:00
#SBATCH --partition=short
#SBATCH --chdir=/scratch/Users/<your_username>/workshop-day7/cellranger_count
#SBATCH --output=/scratch/Users/<your_username>/workshop-day7/cellranger_count/eando/cellrangercount.%j.out
#SBATCH --error=/scratch/Users/<your_username>/workshop-day7/cellranger_count/eando/cellrangercount.%j.err

module load cellranger/7.2.0

echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes
echo Working directory: `pwd`

/usr/bin/time -v cellranger count \
  --id=T21BM_male19 \
  --fastqs=/scratch/Shares/public/sread/data_files/day7a/fastq/sampled_fastq \
  --transcriptome=/scratch/Shares/public/sread/cookingShow/day7a/genomes/refdata-gex-GRCh38-2020-A \
  --sample=DSOX19_1 \
  --include-introns=true \
  --localcores=2 \
  --localmem=10

echo Job finished at `date +"%T %a %d %b %Y"`
```

</details>

Submit:

```bash
sbatch cellrangerCount.sbatch
```

---

## 7. Confirm it is running

```bash
squeue -u $USER
```

Read the **`ST`** column:

| State | Means |
| --- | --- |
| `CF` | **CONFIGURING** — normal here. Wait. |
| `R` | Running |
| `PD` | Pending — waiting for resources to free up |
| `CG` | Completing |

> **`CF` is not an error, and it is not your fault.** This cluster's compute nodes are *dynamic* — note the node name, `short-dy-short-1`. They do not exist until someone asks for one, so SLURM has to boot a machine before your job can start. Expect a couple of minutes in `CF`, especially if you are among the first to submit. **Wait. Do not resubmit.**

Once it flips to `R`, check the log:

```bash
tail -5 eando/cellrangercount.*.err
```

> **A non-empty `.err` is not a failure.** Cell Ranger writes its normal progress to stderr, so text there is expected — you should see stages scrolling past. *Read* it rather than panicking. Real failures announce themselves with `error` or `Pipestance failed`.

If it is in the queue and the log shows progress, **you are done here. Move on to the next worksheet — do not sit and watch it.**

---

## 8. Later: what did it actually use?

When the job finishes, read what `/usr/bin/time -v` recorded:

```bash
grep -E "Maximum resident|Elapsed|Percent of CPU" eando/cellrangercount.*.err
```

- **`Maximum resident set size (kbytes)`** ÷ 1,048,576 = peak GB. Compare that to the 10 GB you requested.
- **`Percent of CPU this job got`** — 200% means you used both cores. 105% means you used one and wasted the other.

This ten-second habit is the difference between requesting resources by superstition and requesting them by measurement — and it is the only honest way to size your *next* job.

> On many clusters you would use `seff <jobid>` for this. **It does not work here**: this cluster has no SLURM accounting database (`sacct` reports *"Slurm accounting storage is disabled"*). Worth knowing — the tool everyone recommends is unavailable on a large fraction of real clusters, whereas `/usr/bin/time -v` works everywhere.

---

## 9. What you get, and how to reuse it

When the run completes, `T21BM_male19/outs/` contains:

- **`web_summary.html`** — open this first, always. Estimated cells, reads per cell, **sequencing saturation**, fraction of reads in cells. It is the first thing a reviewer would look at.
- **`filtered_feature_bc_matrix/`** — `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`. This is what Seurat reads. *Filtered* means Cell Ranger has already decided which barcodes are cells and which are empty droplets — a decision made **for** you, by an algorithm, with a threshold. Which is exactly why the next worksheet does its own QC on top of it.
- **`raw_feature_bc_matrix/`** — every barcode, before that decision.

Full output reference: <https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-overview>

### Adapting this script for your own data

This is the part that outlives the workshop. To run Cell Ranger on your own experiment, you change these — and nothing else:

| Change | To |
| --- | --- |
| `--id` | A name for your sample's output directory |
| `--fastqs` | The directory holding your FASTQs |
| `--sample` | Your FASTQ filename prefix |
| `--transcriptome` | The reference for your species |
| `--ntasks` / `--localcores` | What your cluster allows — likely 8–16 |
| `--mem` / `--localmem` | 64 GB or more for a full-size run |
| `--chdir`, `--output`, `--error` | Your own paths |

Everything else stays exactly as it is. That is the whole point of building it yourself.
