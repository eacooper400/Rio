source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)

install.packages("RColorBrewer")
library(RColorBrewer)

install.packages("pheatmap")
library(pheatmap)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Rio_Analysis/Rio_RNASeq/")

### Check which samples are outliers to be excluded
directory="HTseq/OutlierFix_samples/"
sampleFiles=list.files(directory)
bl=grep("-BL-", sampleFiles)
sampleFiles=sampleFiles[-c(bl)]

sampleNames=gsub("\\.counts\\.txt", '', sampleFiles)
sampleTable=data.frame(sampleName=sampleNames, filename=sampleFiles, condition=sampleNames)
dds=DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
rld=rlog(dds, blind=TRUE)

sampleDists=dist(t(assay(rld)))
sampleDistMatrix=as.matrix(sampleDists)
rownames(sampleDistMatrix)=rld$condition
colors=colorRampPalette(rev(brewer.pal(9,"Blues")))(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

### Completely remove outlier samples from the analysis
### This includes ALL seed tissue or floral tissue samples
### and several leaf tissue (to make balanced matrix)
###remove=grep("-(S|F)[1-3]$", sampleNames)
remove=c(grep("PR22-RI-L2", sampleNames), grep("Rio-RI-L3", sampleNames), grep("PR22-SD-L3", sampleNames), grep("Rio-SD-L3", sampleNames))
sampleNames2=sampleNames[-c(remove)]
sampleNames2=gsub("PR22-V-L1", "Rio-V-L1", sampleNames2)
sampleNames2=gsub("PR22-ANT-L1", "Rio-ANT-L1", sampleNames2) 
sampleFiles2=sampleFiles[-c(remove)]

l=strsplit(sampleNames2, split="-")
sampleGens=unlist(lapply(l, `[[`, 1))
sampleTimes=unlist(lapply(l, `[[`, 2))
sampleTimes=gsub("V", 0, sampleTimes)
sampleTimes=gsub("RI", 1, sampleTimes)
sampleTimes=gsub("FL", 2, sampleTimes)
sampleTimes=gsub("ANT", 3, sampleTimes)
sampleTimes=gsub("SD", 4, sampleTimes)
temp=unlist(lapply(l, `[[`, 3))
l2=strsplit(temp, split='')
sampleTissues=unlist(lapply(l2, `[[`, 1))
sampleTissues=gsub("(S|F)", "M", sampleTissues)
sampleReps=unlist(lapply(l2, `[[`, 2))

sampleTable=data.frame(sampleName=sampleNames2, filename=sampleFiles2, Genotype=sampleGens, Time=sampleTimes, Rep=sampleReps, Tissue=sampleTissues)

### Separate the data by tissue type
intTable=sampleTable[sampleTable$Tissue=="I",1:5]
leafTable=sampleTable[sampleTable$Tissue=="L",1:5]
merTable=sampleTable[sampleTable$Tissue=="M",1:5]

### Construct the DESeq object for each tissue
ddsI=DESeqDataSetFromHTSeqCount(sampleTable=intTable, directory=directory, design= ~ Genotype + Time + Genotype:Time)
ddsL=DESeqDataSetFromHTSeqCount(sampleTable=leafTable, directory=directory, design= ~ Genotype + Time + Genotype:Time)
ddsM=DESeqDataSetFromHTSeqCount(sampleTable=merTable, directory=directory, design= ~ Genotype + Time + Genotype:Time)

### Pre-filter the data by removing genes with 0 or 1 counts
ddsI=ddsI[rowSums(counts(ddsI))>1,]
ddsL=ddsL[rowSums(counts(ddsL))>1,]
ddsM=ddsM[rowSums(counts(ddsM))>1,]

### Perform likelihood ratio tests,
### First remove just the genotype differences,
### but leave genes with significant differences over time
### and genes with a genotype-specific interaction with time
ddsI.1=DESeq(ddsI, test="LRT", reduced= ~ Genotype)
ddsL.1=DESeq(ddsL, test="LRT", reduced= ~ Genotype)
ddsM.1=DESeq(ddsM, test="LRT", reduced= ~ Genotype)

### Next, remove Genotype and Time differences,
### to get just genes with a genotype-specific time interaction
ddsI.2=DESeq(ddsI.1, test="LRT", reduced= ~ Genotype + Time)
ddsL.2=DESeq(ddsL.1, test="LRT", reduced= ~ Genotype + Time)
ddsM.2=DESeq(ddsM.1, test="LRT", reduced= ~ Genotype + Time)

### Print out all results to tables
write.table(results(ddsI.1), file="DESeq_TimeSeries_2/Internode_TimeDiff.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
write.table(results(ddsI.2), file="DESeq_TimeSeries_2/Internode_GxTixn.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
write.table(results(ddsL.1), file="DESeq_TimeSeries_2/Leaf_TimeDiff.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
write.table(results(ddsL.2), file="DESeq_TimeSeries_2/Leaf_GxTixn.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
write.table(results(ddsM.1), file="DESeq_TimeSeries_2/Meristem_TimeDiff.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
write.table(results(ddsM.2), file="DESeq_TimeSeries_2/Meristem_GxTixn.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

### To test for significance at a single timepoint:
### for example, Rep. Initiation in internode:
resRI <- results(ddsI.2, name="GenotypeRio.Time1", test="Wald")

### Cluster genes by their profiles and plot:
betas <- coef(ddsI.2)
topGenes <- head(order(results(ddsI.2)$padj),30)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),cluster_col=FALSE)

betas <- coef(ddsI.1)
topGenes <- head(order(results(ddsI.1)$padj),30)

betas <- coef(ddsL.2)
topGenes <- head(order(results(ddsL.2)$padj),30)

betas <- coef(ddsM.2)
topGenes <- head(order(results(ddsM.2)$padj),30)

### Save genotype and tissue specific variance stabilized counts for WGCNA
intTableR=intTable[intTable$Genotype=="Rio",]
intTableP=intTable[intTable$Genotype=="PR22",]
leafTableR=leafTable[leafTable$Genotype=="Rio",]
leafTableP=leafTable[leafTable$Genotype=="PR22",]
merTableR=merTable[merTable$Genotype=="Rio",]
merTableP=merTable[merTable$Genotype=="PR22",]

ddsIR=DESeqDataSetFromHTSeqCount(sampleTable=intTableR, directory=directory, design= ~ Time)
ddsIP=DESeqDataSetFromHTSeqCount(sampleTable=intTableP, directory=directory, design= ~ Time)
ddsLR=DESeqDataSetFromHTSeqCount(sampleTable=leafTableR, directory=directory, design= ~ Time)
ddsLP=DESeqDataSetFromHTSeqCount(sampleTable=leafTableP, directory=directory, design= ~ Time)
ddsMR=DESeqDataSetFromHTSeqCount(sampleTable=merTableR, directory=directory, design= ~ Time)
ddsMP=DESeqDataSetFromHTSeqCount(sampleTable=merTableP, directory=directory, design= ~ Time)

vsdIR=varianceStabilizingTransformation(ddsIR, blind=TRUE)
vsdIP=varianceStabilizingTransformation(ddsIP, blind=TRUE)
vsdLR=varianceStabilizingTransformation(ddsLR, blind=TRUE)
vsdLP=varianceStabilizingTransformation(ddsLP, blind=TRUE)
vsdMR=varianceStabilizingTransformation(ddsMR, blind=TRUE)
vsdMP=varianceStabilizingTransformation(ddsMP, blind=TRUE)

write.table(assay(vsdIR), file="DESeq_TimeSeries_2/IntRio_vsd.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
write.table(assay(vsdIP), file="DESeq_TimeSeries_2/IntPR22_vsd.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
write.table(assay(vsdLR), file="DESeq_TimeSeries_2/LeafRio_vsd.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
write.table(assay(vsdLP), file="DESeq_TimeSeries_2/LeafPR22_vsd.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
write.table(assay(vsdMR), file="DESeq_TimeSeries_2/MerRio_vsd.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
write.table(assay(vsdMP), file="DESeq_TimeSeries_2/MerPR22_vsd.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

### Save all of the differential expression objects to access later to compare with WGCNA
save(ddsI.1, ddsI.2, ddsL.1, ddsL.2, ddsM.1, ddsM.2, file="DESeq_TimeSeries_2/all_dds_objects.RData")
