### Normalization and correctioon####
### Unaltered RNA-seq####
library(DESeq2)
library(tidyverse)
RNAmetadata=read.table("/Users/samuelhunter/OneDrive - UCB-O365/sread2022/extra_talk/RNAinfo.txt", 
                       sep="\t", header=TRUE)
common_names <-"refseq_to_common_id.txt"
masterannotationdf=read.table("/Users/samuelhunter/OneDrive - UCB-O365/sread2022/extra_talk/annotation.csv",
                              sep=",",header = TRUE)
masterannotationdf_only21 <- masterannotationdf %>% filter(chr=="chr21")
aneuploidygenes <-masterannotationdf_only21[["name"]] #Only chr21 genes
common_ids <- read.table(paste0("/Users/samuelhunter/OneDrive - UCB-O365/sread2022/extra_talk/",common_names),
                         sep="\t")
baseploidy <- 2
alt_ploidy <-3
lesschrs <- c("chr20","chr21")
RNAcountdat <- read.csv("/Users/samuelhunter/OneDrive - UCB-O365/sread2022/extra_talk/RNAseqcount.txt",
                        sep=",", skip=0)
rownames(RNAcountdat) <- RNAcountdat$X
RNAcountdat  <- RNAcountdat[,-c(1)]
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, 
                                  design = ~ biological_rep + Person)

person1="Ethan"
person2="Eric"
ddsFull <- DESeq(ddsFull)
RNAddsres<-results(ddsFull,contrast=c("Person",person1,person2))
resdata <- as.data.frame(RNAddsres)

RNAddsresshrunk <- lfcShrink(ddsFull,
                       contrast = c("Person",person1,person2), res=RNAddsres, type = 'normal')
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=1)
head(fullresdata)
common_fullresdata <- merge(fullresdata,common_ids,by.x=1,by.y=1)
head(common_fullresdata)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- plyr::ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  theme_classic(base_size = 30) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, 
              position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, 
              position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, 
               geom="crossbar", size=0.25, color="orange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank())



#Volcano Plot
signif_fullresdata <- fullresdata[fullresdata$padj<.01,]
notsignif_fullresdata <- fullresdata[fullresdata$padj>=.01,]
fullresdata_chr21 <- fullresdata[fullresdata$chr=="chr21",]

fullresdata <- na.omit(fullresdata)
head(fullresdata)
ggplot() + 
  geom_point(data=fullresdata,aes(x=log2FoldChange, y=-log10(pvalue)),size=1,alpha=0.5) +
  geom_point(data=signif_fullresdata,aes(x=log2FoldChange, y=-log10(pvalue)),size=1,alpha=0.5,color="red") +
  geom_point(data=fullresdata_chr21,aes(x=log2FoldChange, y=-log10(pvalue)),size=1,alpha=0.5,color="green") +
  ylab(paste0("-log10(pvalue)")) +
  geom_vline(xintercept = log2(1.5), linetype="dotted", color="blue" ) +
  geom_vline(xintercept = log2(1), linetype="dotted", color="black" ) +
  xlab("log2FC") +
  theme_classic(base_size = 24) +
  theme(legend.title = element_blank())


### Setting up our normalization matrix...####
RNAmetadata
controlgenes <-ifelse(rownames(ddsFull) %in% aneuploidygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes
ploidy_trisomy21 = ifelse(rownames(ddsFull) %in% aneuploidygenes, alt_ploidy, baseploidy)
ploidy_typical <- rep(baseploidy, nrow(ddsFull))

normFactors <- matrix(do.call(cbind, mget(paste0(ddsFull$ploidy))),
                      ncol=ncol(ddsFull),nrow=nrow(ddsFull),
                      dimnames=list(1:nrow(ddsFull),
                                    1:ncol(ddsFull)))
head(RNAcountdat)
normFactors <- normFactors/baseploidy

unique(normFactors)


# Normalizing by Ploidy #

ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, 
                                  colData = RNAmetadata, design = ~ biological_rep + Person)
ddsCollapsed_normfactor <-estimateSizeFactors(ddsFull,normMatrix=normFactors, 
                                              controlGenes=controlgenes) 
DESeq(ddsFull)
ddsCollapsed_normfactor <- estimateDispersionsGeneEst(ddsCollapsed_normfactor) 
ddsCollapsed_normfactor<-estimateDispersionsFit(ddsCollapsed_normfactor) 
ddsCollapsed_normfactor <- estimateDispersionsMAP(ddsCollapsed_normfactor)  
ddsCollapsed_normfactor <- nbinomWaldTest(ddsCollapsed_normfactor) 

person1 = "Ethan"
person2 = "Eric"
RNAddsres<-results(ddsCollapsed_normfactor,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=1)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]

ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  theme_classic(base_size = 30) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, 
              position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, 
              position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, 
               geom="crossbar", size=0.25, color="orange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank())

###RNA-seq Change Hypothesis ####

# Changing alternative hypothesis #
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, 
                                  colData = RNAmetadata, design = ~ biological_rep + Person)
ddsFull <- DESeq(ddsFull)
person1 = "Ethan"
person2 = "Eric"
RNAddsres<-results(ddsFull,  contrast=c("Person", person1, person2),
                   altHypothesis = "lessAbs", lfcThreshold=log2(1.5))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=1)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
nrow(fullresdata)
nrow(fullresdata[fullresdata$padj<0.01 & fullresdata$chr == "chr21",])
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

signif_medresdata <- medresdata[medresdata$padj<.1,]
notsignif_medresdata <- medresdata[medresdata$padj>=.1,]

ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC (RNA-seq (T21 vs D21))")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  theme_classic(base_size = 30) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, 
              position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, 
              position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, 
               geom="crossbar", size=0.25, color="orange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank())

