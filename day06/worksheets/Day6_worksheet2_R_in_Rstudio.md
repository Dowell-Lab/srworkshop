# Short Read Day 6 Worksheet | R in RStudio

Authors: Georgia Barone (2023), Rutendo Sigauke (2024)

## Section A : Learning R in RStudio

-  Before beginning this worksheet, make sure RStudio is downloaded on your local computer.

1. Go to the srworkshop [day 6](https://github.com/Dowell-Lab/srworkshop/tree/main) GitHub page 

2. Click on the folder that is called “scripts”

3. Click on the script titled “Learning_R.R”

4. Download and open Learning_R.R in RStudio

5. Complete Learning_R.R

## Section B : Writing an R script to submit on a supercomputer

1. Log in to AWS

2. cd into the srworkshop repository and run a `git pull`

3. cd into /scratch/Users/<your_username>/

4. Use the mkdir command to make a folder called `day6`

5. Inside your `day6` folder, make the folders `scripts`, `eofiles`, and `results`

6. Use the rsync or scp command to copy Learning_R_submit_aws.R and Submit_Rscript.sbatch
from `/Users/<your_username>/srworkshop/day06/scripts`, into the scripts directory in
`/scratch/Users/<your_username>/day6` you just made.

7. Go to `/scratch/Users/<your_username>/day6/scripts`. Use vim to open and edit the
`Learning_R_submit_aws.R` file. Add your own working directory path (it should be something
along the lines of: `/scratch/Users/<your_username>/day6/results`).

8. Look through the rest of the script before saving and exiting vim, to make sure you know what
the code is doing and where your output will be saving to.

9. Now use vim to open/edit the `Submit_Rscript.sbatch` file. This is the sbatch script we will be
using to submit our code to the supercomputer. Edit the script by adding your job name, email,
eofiles path, and path to Learning_R_submit_aws.R.

- Note: the command Rscript runs R scripts or R commands directly from the bash shell

10. Once you are happy with your `Learning_R_submit_aws.R` and `Submit_Rscript.sbatch` scripts,
submit the `Submit_Rscript.sbatch` script.

11. If the script worked, back up your results to your home directory
(`/Users/<your_username>/srworkshop/day06`).

12. Use rsync to copy mtcars.csv & mtcars_mpg_wt_scatterplot.png from AWS to your local
computer to view.