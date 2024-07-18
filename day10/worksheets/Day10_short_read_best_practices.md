# Short Read Analysis Best Practices
By Mary Allen and Robin Dowell
 
## Goal of any analysis: Discover some cool results while maintaining integrity and reproducibility.

## Download your data from a sequencing facility or a repository
- Lock the permissions of raw files to be read-only so they can never be edited, even accidentally (i.e. `chmod 444`)!
- If you generated the raw data, keep it somewhere that is backed up!
- Note where you got the data from (facility, machine, collaborator, database, etc.)
    - If published data, note the citation of the paper
- Make a meta data table, especially if you generated the data
    - How it was generated
    - Cell types and genetic modifications
    - Perturbations
    - Library kit and sequencing specifics
- If you generated data, consider uploading to GEO/SRA immediately. You can keep the project private for years until you publish.

## Run your analysis
- Set up your storage system optimally
    -`/scratch/` vs `/Users/` - fast vs backed up
    - Make a directory on `/scratch` for each of your projects (`/scratch/Shares/<labname>/` or `/scratch/Users/<username>/`)
    - Make an input and output directory in your project directory
        - Rsync your raw data to your input directory on scratch
            - Scratch is not backed up!
        - Keep a README or NOTES file with the path of the raw data and when you copied it over. If you do any massaging of your data, record it.
    - Make a scripts directory 
        - Use version control, even if just locally
        - Be sure your scripts are backed up, either on the cluster or on GitHub
- All software you run should be in a script (not on command line!)
    - Make a README file that documents everything you would put in your lab notebook
        - Where does the data live (raw and important intermediates)?
        - Which scripts did you run on it (and why)?
        - What files did you create and where are they?
        - What versions of software, genomes, annotations, etc were used?
- Keep the living room clean!
    - In pipelines, there will often be intermediate files (output of program A used as input for program B). Actively manage these intermediates.
        - When you get intermediate files you want to backup, rsync them to somewhere backed up - otherwise, delete as soon as you don’t need them.   
        - Delete stuff on `/scratch` periodically (Data on scratch costs more and clogs up the system).
- Always, always, sanity check your results
    - Make sure files are around the size you expect
        - Did the program work at all?
        - Did the program process ALL of the data?
    - For non-binary files, look at a subset of lines to see if they make sense
        - Are they interpretable?
        - Can you manually verify results with true positives?
    - If possible, ALWAYS look at the data in IGV!
        - Mapping data
        - Coordinate data like bedfiles
    - Don’t believe people when they say their program does a specific thing (check!)

## Publish your data
- Upload the raw data, meta-data and the final processed files to the NIH GEO/SRA databases
    - Uploading to GEO should automatically upload to the SRA
    - Again, it's not a bad idea to upload soon after receiving data from the sequencer
- You must report all manipulations of data (manual or through analysis tools)
- All versions of all programs used must be noted and the paper associated with the program should be cited in the methods section. Example: “We used Tfit (Azofeifa 2017) v 1.1 to identify eRNAs.”
    - For reproducibility, you MUST note all the flags/options you used if not standard (i.e. defaults are assumed)
        - You can do this in your scripts on GitHub and provide the “source” code for version information (must provide GitHub link in methods and intend to maintain the repository).
        - You can also use Jupyter notebooks to document analysis, plotting, manipulations and versions. Provide the Jupyter notebook in the methods.
        - Or you can simply note all the flags in the methods section - this works well if there are few of them, but gets a bit unwieldy if you used a lot of software and unique options.