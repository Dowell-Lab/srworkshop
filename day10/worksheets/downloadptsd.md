# PTSD FASTQ Download Pipeline (Fiji, SLURM)

This repository contains scripts to download PTSD FASTQ files on the Fiji cluster using SLURM.  
You will:

1. Copy the scripts to your scratch space.
2. Edit paths and email addresses.
3. Run a driver script that submits one SLURM job per FASTQ URL.

---

## 1. Set up your scratch directory

Log in to Fiji:

```bash
ssh your_identikey@fiji.colorado.edu
```

Go to your scratch directory and create a project folder:

```bash
cd /scratch/Users/your_identikey
mkdir -p disease_datasets/ptsd
cd disease_datasets/ptsd
```

Copy the four files provided by your PI (ask for the exact source path if different):

```bash
cp /scratch/Users/ispe1418/disease_datasets/ptsd/run_wholeline.sh .
cp /scratch/Users/ispe1418/disease_datasets/ptsd/header.txt .
cp /scratch/Users/ispe1418/disease_datasets/ptsd/acommandsbatch.sbatch .
cp /scratch/Users/ispe1418/disease_datasets/ptsd/ptsd_download_links .
```

Create subdirectories and move scripts:

```bash
mkdir -p scripts fastq eando
mv run_wholeline.sh header.txt acommandsbatch.sbatch ptsd_download_links scripts/
cd scripts
```

Now all scripts live in:

```text
/scratch/Users/your_identikey/disease_datasets/ptsd/scripts
```

---

## 2. Edit `run_wholeline.sh` (paths)

Open the file with `vim`:

```bash
vim run_wholeline.sh
```

At the top you should see something like:

```bash
indir=/scratch/Users/ispe1418/disease_datasets/ptsd/scripts/
wholefileoflines=${indir}ptsd_download_links
outdir=/scratch/Users/ispe1418/disease_datasets/ptsd/fastq/
script_dir=/scratch/Users/ispe1418/disease_datasets/ptsd/scripts/
```

Change **all** occurrences of `ispe1418` to your own username:

```bash
indir=/scratch/Users/your_identikey/disease_datasets/ptsd/scripts/
wholefileoflines=${indir}ptsd_download_links
outdir=/scratch/Users/your_identikey/disease_datasets/ptsd/fastq/
script_dir=/scratch/Users/your_identikey/disease_datasets/ptsd/scripts/
```

Do **not** change the loop or the `sbatch --export=... acommandsbatch.sbatch` line.

Save and quit:

- Press `Esc`
- Type `:wq`
- Press Enter

---

## 3. Edit `header.txt` (SBATCH options)

Open:

```bash
vim header.txt
```

You should see:

```bash
#!/bin/bash

#SBATCH --job-name=get_data_wget
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --mem=50mb
#SBATCH --output=/scratch/Users/ispe1418/disease_datasets/ptsd/eando/%x_%j.out
#SBATCH --error=/scratch/Users/ispe1418/disease_datasets/ptsd/eando/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ispe1418@colorado.edu
```

You must change:

1. Output and error paths (username):

   ```bash
   #SBATCH --output=/scratch/Users/your_identikey/disease_datasets/ptsd/eando/%x_%j.out
   #SBATCH --error=/scratch/Users/your_identikey/disease_datasets/ptsd/eando/%x_%j.err
   ```

2. Email address:

   ```bash
   #SBATCH --mail-user=your_identikey@colorado.edu
   ```

Leave all other `#SBATCH` lines as‑is unless instructed otherwise by your PI.

Save and quit with `:wq`.

---

## 4. `ptsd_download_links` (usually no edits)

This file contains one `wget` command per FASTQ file, for example:

```bash
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR541/003/SRR5412743/SRR5412743.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR541/005/SRR5412675/SRR5412675.fastq.gz
...
```

Normally, you **do not** change this file.  
Only edit it if you want to download a different set of URLs.

To look at it:

```bash
vim ptsd_download_links
```

Quit without saving:

```text
:q
```

---

## 5. Make scripts executable

From the `scripts` directory:

```bash
chmod +x run_wholeline.sh
chmod +x acommandsbatch.sbatch
```

---

## 6. Run the download driver script

In the `scripts` directory:

```bash
bash run_wholeline.sh
```

This script will:

- Count how many lines are in `ptsd_download_links`.
- Submit one SLURM job per line using `acommandsbatch.sbatch`.

To monitor your jobs:

```bash
squeue -u your_identikey
```

FASTQ files will appear here:

```bash
/scratch/Users/your_identikey/disease_datasets/ptsd/fastq
```

Log files will appear here:

```bash
/scratch/Users/your_identikey/disease_datasets/ptsd/eando
```

---

## 7. Tiny vim cheat sheet

- Open a file: `vim filename`
- Move: arrow keys
- Start editing: press `i`
- Stop editing: press `Esc`
- Save and quit: `:wq` then Enter
- Quit without saving: `:q!` then Enter
