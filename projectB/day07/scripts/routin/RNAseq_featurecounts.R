#Feature counts for RNA-seq

library("Rsubread")

args = commandArgs(trailingOnly=TRUE)

refseqfile="/scratch/Shares/dowell/genomes/hg38/hg38_refseq_genenames_included.gtf"

#prefix="RNA-BSA-Eric-1"
prefix=args[1]
indir=args[2]
outdir=args[3]
bamfile=paste(indir, prefix, ".sorted.bam", sep="")

fc <- featureCounts(files=bamfile,
                    annot.ext=refseqfile,
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="exon",
                    GTF.attrType="transcript_id",
                    useMetaFeatures=TRUE,
                    allowMultiOverlap=TRUE,
                    largestOverlap=TRUE,
                    countMultiMappingReads=TRUE,
                    isPairedEnd=TRUE,
                    strandSpecific=2,
                    nthreads=8)


fc$annotation["TranscriptID"] <- fc$annotation$GeneID
write.table(x=data.frame(fc$annotation[,c("GeneID","TranscriptID","Length")],
                         fc$counts,stringsAsFactors=FALSE),
            file=paste0(outdir,prefix,".stranded.transcript_counts.txt"),
            quote=FALSE,sep="\t",
            row.names=FALSE)

fc <- featureCounts(files=bamfile,
                    annot.ext=refseqfile,
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="exon",
                    GTF.attrType="transcript_id",
                    useMetaFeatures=TRUE,
                    allowMultiOverlap=TRUE,
                    largestOverlap=TRUE,
                    countMultiMappingReads=TRUE,
                    isPairedEnd=TRUE,
                    strandSpecific=1,
                    nthreads=8)

fc$annotation["TranscriptID"] <- fc$annotation$GeneID
write.table(x=data.frame(fc$annotation[,c("GeneID","TranscriptID","Length")],
                         fc$counts,stringsAsFactors=FALSE),
            file=paste0(outdir,prefix,".oppstranded.transcript_counts.txt"),
            quote=FALSE,sep="\t",
            row.names=FALSE)


fc <- featureCounts(files=bamfile,
                    annot.ext=refseqfile,
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="exon",
                    GTF.attrType="transcript_id",
                    useMetaFeatures=TRUE,
                    allowMultiOverlap=TRUE,
                    largestOverlap=TRUE,
                    countMultiMappingReads=TRUE,
                    isPairedEnd=TRUE,
                    strandSpecific=0,
                    nthreads=8)

fc$annotation["TranscriptID"] <- fc$annotation$GeneID
write.table(x=data.frame(fc$annotation[,c("GeneID","TranscriptID","Length")],
                         fc$counts,stringsAsFactors=FALSE),
            file=paste0(outdir,prefix,".unstranded.transcript_counts.txt"),
            quote=FALSE,sep="\t",
            row.names=FALSE)
