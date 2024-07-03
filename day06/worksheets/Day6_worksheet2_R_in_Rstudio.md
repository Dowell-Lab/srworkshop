# Short Read Day 6 Worksheet | R in RStudio

Authors: Georgia Barone (2023), Rutendo Sigauke (2024)

## Section A : Learning R in RStudio

-  Before beginning this worksheet, make sure RStudio is downloaded on your local computer.

1. Go to the [srworkshop](https://github.com/Dowell-Lab/srworkshop/tree/main) on GitHub to update your local repository.

This section can be done on your personal computer or on the AWS vizualization instance.

2. Click on the folder that is called “scripts”

3. Click on the script titled “Learning_R.R”

4. Download and open Learning_R.R in RStudio

![load in RStudio](ws2_open_file_in_Rstudio.png)

5. Complete Learning_R.R

## Section B : Writing an R script to submit on a supercomputer

1. Log in to AWS

2. Move into the `srworkshop` repository and update from repository by running `git pull`

![Git pull](ws2_git_pull.png)

3. Go into /scratch/Users/<your_username>/ and use the mkdir command to make a folder called `day6`

- Inside your `day6` folder, make the folders `scripts`, `eofiles`, and `results`

![Initialize project](ws2_initialize_folders.png)

4. Use the rsync or scp command to copy Learning_R_submit_aws.R and Submit_Rscript.sbatch
from `/Users/<your_username>/srworkshop/day06/scripts`, into the scripts directory in
`/scratch/Users/<your_username>/day6` you just made.

5. Go to `/scratch/Users/<your_username>/day6/scripts`. Use vim to open and edit the
`Learning_R_submit_aws.R` file. Add your own working directory path (it should be something
along the lines of: `/scratch/Users/<your_username>/day6/results`).

![Edit R Script](ws2_edit_Rscript_vim.png)

6. Look through the rest of the script before saving and exiting vim, to make sure you know what
the code is doing and where your output will be saving to.

7. Now use vim to open/edit the `Submit_Rscript.sbatch` file. This is the sbatch script we will be
using to submit our code to the supercomputer. Edit the script by adding your job name, email,
eofiles path, and path to Learning_R_submit_aws.R.

- Note: the command Rscript runs R scripts or R commands directly from the bash shell

![Sbatch R Script](ws2_sbatch_to_submit_R_script.png)

8. Once you are happy with your `Learning_R_submit_aws.R` and `Submit_Rscript.sbatch` scripts,
submit the `Submit_Rscript.sbatch` script.

9. If the script worked, back up your results to your home directory
(`/Users/<your_username>/srworkshop/day06`).

10. Use rsync to copy mtcars.csv & mtcars_mpg_wt_scatterplot.png from AWS to your local
computer to view.