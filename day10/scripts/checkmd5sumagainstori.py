import pandas as pd
import sys


def addcolcheck(checkfile, datadatframe):
    if datadatframe.find(".tsv")>-1:
        df = pd.read_csv(datadatframe, sep="\t")
    else:
        df = pd.read_csv(datadatframe)
    md5sumsdf = pd.read_csv(checkfile, sep=" ", names=["fastq_md5", "fastq_file"])

    print("dfcolumns", df.columns)
    print("md5columns", md5sumsdf.columns)
    # Make sure fastq_md5 is string
    df["fastq_md5"] = df["fastq_md5"].astype(str)

    # Split on semicolon into lists, then explode to get one md5 per row
    df = df.assign(
        fastq_md5=df["fastq_md5"].str.split(";")
    ).explode("fastq_md5")

    # Strip whitespace just in case
    df["fastq_md5"] = df["fastq_md5"].str.strip()
    md5sumsdf["fastq_md5"] = md5sumsdf["fastq_md5"].astype(str).str.strip()

    # Now do a regular exact match
    df["md5summatch"] = df["fastq_md5"].isin(md5sumsdf["fastq_md5"])

    print("number of fastq files whose md5sums match", df["md5summatch"].value_counts())
    outname = checkfile + ".in_original_check.exploded.csv"
    df.to_csv(outname, index=False, sep="\t")
    print(outname)




#checkfile = "/Users/allenma/CNOT_and_brain/EDC3/qc/afterdownload.md5sum.txt"
#datadatframe= "/Users/allenma/CNOT_and_brain/EDC3/metadata/filereport_read_run_PRJNA436576.tsv"

if __name__ == "__main__":
    checkfile = sys.argv[1]
    datadatframe = sys.argv[2]
    addcolcheck(checkfile, datadatframe)
