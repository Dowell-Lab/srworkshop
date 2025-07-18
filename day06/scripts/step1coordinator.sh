CSV_FILE="/scratch/Shares/public/sread2025/cookingShow/day6/populationRNAseq/genesymbolchr.csv"

# Skip the header and loop through each line of the CSV file
tail -n +2 "$CSV_FILE" | while IFS=',' read -r index ensembl chr symbol
do
    # Remove any leading or trailing spaces (if any)
    symbol=$(echo "$symbol" | xargs)

    # Export SYMBOL as a variable to pass to the sbatch script
    export agenename="$symbol"

    # Submit the job via sbatch, passing the variable agenename
    sbatch --export=agenename "/Users/maallen3/srworkshop/day06/scripts/step1_gene_Submit_Rscript.sbatch"

    # Optional: Print the symbol being passed
    echo "Submitting job for gene: $symbol"
break #once the whole thing works for the first gene--- remove this break. It only lets the script loop once until then. 
done
