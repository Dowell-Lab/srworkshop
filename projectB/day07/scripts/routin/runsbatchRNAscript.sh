indir=/Shares/down/interferonbeta/analysis_Jan2023_MAA/outfiles/RNA-seq/BedgraphsandBigwigs_human/mapped/bams/
outdir=/Shares/down/interferonbeta/analysis_May2023_MAA/outfiles/RNA-seq/annotation_gene_counts/


for pathandfilename in `ls ${indir}*sorted.bam`; do
rootname=`basename $pathandfilename .sorted.bam`
sbatch --export=rootname=$rootname,indir=$indir,outdir=$outdir sbatchRNAscript.sh 
echo $rootname
done


