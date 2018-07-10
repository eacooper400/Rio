### Installing Required Packages (only done once)
source("https://bioconductor.org/biocLite.R")
biocLite("EBSeq")
biocLite("EBSeqHMM")

### Load Required Pacakges
library(data.table)
library(DESeq2)
library(RColorBrewer)
library(EBSeq)
library(EBSeqHMM)

### Set WorkDir and load my own R functions
setwd("~/Desktop/Rio_Analysis/Rio_RNASeq/")
source("scripts/RNAseq_fxns.R")

### Set variables that define conditions
my.tissues=c("Internode", "Leaf", "Meristem")
my.genotypes=c("PR22", "Rio")
my.times=c("V", "RI", "FL", "ANT", "SD")
my.colors=c("#D95F02","#1F78B4")

### PART ONE: USE DESeq RESULT TO SCREEN FOR GENES
### DIFFERENTIALLY EXPRESSED BETWEEN THE TWO GENOTYPES
lnames=load("DEseq2_objectsBYTissue.RData")
load("DESeq2_altMeristem.RData")
deseq.list=list(ddsI, ddsL, ddsM)
names(deseq.list)=my.tissues
for (t in 1:3) {

    ## PART 1a: Screen for genes that are differentially expressed between
    ## genotypes at ANY time point
    d=deseq.list[[t]]
    table1=dge.anyTime(dds=d, times=my.times)

    ## Part 1b: Screen for genes with significant GxT interactions 
    table1$GT.pvalue=dge.interaction(dds=d, genes=rownames(table1))
    table2=subset(table1, table1$GT.pvalue<0.05)

    ## Part 1c: Classify each remaining gene based on WHEN it is
    ## differentially expressed between the genotypes
    ## Remove any genes with NA values
    table3=pattern.category(table2[,1:5])
    table3$GT.pvalue=table2$GT.pvalue
    to.remove=grep("NA", table3$category)
    table3=table3[-(to.remove),]
    
    ## Part 1d: Retrieve the raw counts for the filtered genes
    ## Save output to files to access later
    m=match(rownames(table3), rownames(d))
    filt.counts=counts(d)[m,]
    orderSamples=c("PR22-V-*", "PR22-RI*", "PR22-FL*", "PR22-ANT*", "PR22-SD*", "Rio-V-*", "Rio-RI*", "Rio-FL*", "Rio-ANT*", "Rio-SD*")
    orderColumns=unlist(lapply(orderSamples, FUN=function(x) grep(x, colnames(filt.counts))))
    Counts=filt.counts[,orderColumns]

    out1=paste(c("EBSeq/Sig_GxT_Genes", "Pvalues", my.tissues[t], "txt"), collapse=".")
    out2=gsub("Pvalues", "Counts", out1)
    write.table(table3, file=out1, quote=TRUE, row.names=TRUE, col.names=TRUE, sep="\t")
    write.table(Counts, file=out2, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")    
}

### PART TWO: USE EBSeq TO FIND GENES DIFFERENTIALLY EXPRESSED
### OVER TIME WITHIN EACH GENOTYPE AND CLUSTER THEM BY EXPRESSION PATTERN
for (t in 1:3) {
    infile=paste(c("EBSeq/Sig_GxT_Genes", "Counts", my.tissues[t], "txt"), collapse=".")
    my.counts=read.table(infile, header=TRUE, stringsAsFactors=FALSE)

    ## Part 2a: Separate the data by genotype, and format counts for EBSeq
    for (g in 1:2) {
        gcounts=my.counts[,grep(my.genotypes[g], colnames(my.counts))]
        gcounts=as.matrix(gcounts)
        condVector=unlist(lapply(strsplit(colnames(gcounts), split="\\."), `[[`, 2))
        Conditions=factor(condVector, levels=my.times)

        ## Part 2b: Adjust the counts by library size
        Sizes=MedianNorm(gcounts)
        gnorm <- GetNormalizedMat(gcounts, Sizes)
        out1=paste(c("EBSeq/Normalized_Matrix", my.genotypes[g], my.tissues[t], "txt"), collapse=".")
        write.table(gnorm, out1, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

        ## Part 2c: Run EBSeqHMM Test to get Posterior Probs for each gene belonging to each possible pattern cluster
        EBSeqHMMGeneOut <- EBSeqHMMTest(Data=gcounts, sizeFactors=Sizes, Conditions=Conditions, UpdateRd=10)
        out2=gsub("Normalized", "PostProb", out1)
        write.table(EBSeqHMMGeneOut$PPMatZ, out2, quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

        ## Part 2d: Distinguish genes Diff. Expressed over time (FDR=0.05) from those that appear constant over time
        ## Save the "best Cluster" for each gene (including ones not DE)
        GeneDECalls <- GetDECalls(EBSeqHMMGeneOut, FDR=.05)
        gclusters=merge(x=gcounts, y=GeneDECalls, by.x=0, by.y=0, all.x=TRUE)
        out3=gsub("PostProb_Matrix", "Clusters", out2)
        write.table(gclusters, out3, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

        ## Part 2e: For each cluster, plot the DE genes that have been assigned with confidence
        GeneConfCalls <- GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05,cutoff=.5, OnlyDynamic=FALSE)
        cluster.size=unlist(lapply(GeneConfCalls$EachPath, nrow))
        to.plot=cluster.size[which(cluster.size>0)]
        plotfile=gsub("txt", "pdf", out3)
        pdf(plotfile, width=8, height=11)   
        layout(matrix(seq(1,6), byrow=TRUE, ncol=2))
        for (i in 1:length(to.plot)) {
            path=GeneConfCalls$EachPath[[names(to.plot)[i]]]
            m1=match(rownames(path), rownames(gnorm))
            path.data=gnorm[m1,]
            if (length(m1)==1) {
                path.data=t(as.matrix(gnorm[m1,], nrow=1, byrow=T))
            }
            ebseq.clusterPlot(counts=path.data, times=get.times(gnorm), plot.color=my.colors[g], name=(names(to.plot)[i]))
        }
        dev.off()
    }
}

### PART THREE: COMPARE CLUSTERING FOR THE TWO GENOTYPES
### SCORE GENES BY HOW MUCH THEIR PATTERNS CHANGE
for (t in 1:3) {
    in1=paste(c("EBSeq/Clusters", my.genotypes[1], my.tissues[t], "txt"), collapse=".")
    in2=paste(c("EBSeq/Clusters", my.genotypes[2], my.tissues[t], "txt"), collapse=".")
    pr22.clust=read.table(in1, header=TRUE, stringsAsFactors=FALSE, sep="\t")
    rio.clust=read.table(in2, header=TRUE, stringsAsFactors=FALSE, sep="\t")
    rio.clust=rio.clust[, c(1,(ncol(rio.clust)-1):(ncol(rio.clust)))]
    pr22.clust$Most_Likely_Path[is.na(pr22.clust$Most_Likely_Path)]="EE-EE-EE-EE"
    rio.clust$Most_Likely_Path[is.na(rio.clust$Most_Likely_Path)]="EE-EE-EE-EE"

    ## Part 3a: Create a score for every possible pair of clusters
    all.clusters=unique(c(pr22.clust$Most_Likely_Path, rio.clust$Most_Likely_Path))
    all.clusters=all.clusters[!(is.na(all.clusters))]
    cmatrix=score.clusters(all.clusters)
    
    ## Part 3b: For each gene, generate a "difference score"
    ## based on its prob. of belonging to a pair of clusters
    ## and the difference score of those clusters    
    clust.compare=data.frame(pr22.clust$Row.names, pr22.clust$Most_Likely_Path, pr22.clust$Max_PP, rio.clust$Most_Likely_Path, rio.clust$Max_PP)
    score.rows=match(clust.compare[,2], rownames(cmatrix))
    score.cols=match(clust.compare[,4], colnames(cmatrix))
    clust.compare$Cscore=vapply(1:nrow(clust.compare), FUN=function(x,y,z,m) m[y[x],z[x]],
                                FUN.VALUE=double(1), USE.NAMES=FALSE, y=score.rows, z=score.cols, m=cmatrix)
    clust.compare$DiffScore=clust.compare$pr22.clust.Max_PP * clust.compare$rio.clust.Max_PP * clust.compare$Cscore

    ## Part 3c: Generate a final table, with the GxT interaction terms,
    ## the expression patterns for each genotype, and the Difference score
    ## for each gene
    in3=paste(c("EBSeq/Sig_GxT_Genes", "Pvalues", my.tissues[t], "txt"), collapse=".")
    gbyt=read.table(in3, stringsAsFactors=FALSE, header=TRUE)
    final=merge(clust.compare, gbyt, by.x=1, by.y=0)
    rownames(final)=final[,1]
    final[,1]=NULL
    final=final[,c(1,3,6,13)]
    colnames(final)=c("PR22_Path", "Rio_Path", "DiffScore", "GbyT_pValue")
    out1=paste(c("EBSeq/ClusterCompare_by_Gene", my.tissues[t], "txt"), collapse=".")
    write.table(final, out1, quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

    ## Part 3d: Find genes with a different pattern at Rep. Init.
    pr22.pattern=as.vector(final$PR22_Path, mode="character")    
    p.ri=unlist(lapply(strsplit(pr22.pattern, split="-"), `[[`, 1))
    rio.pattern=as.vector(final$Rio_Path, mode="character")
    r.ri=unlist(lapply(strsplit(rio.pattern, split="-"), `[[`, 1))
    ri.change=vapply(1:length(p.ri), FUN=function(x,y,z) length(grep(y[x], z[x])), FUN.VALUE=integer(1), USE.NAMES=FALSE, y=p.ri, z=r.ri)
    ri.genes=final[which(ri.change==0),]
    out2=gsub("by_Gene", "RepInitChange", out1)
    write.table(ri.genes, out2, quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")
}
    
new=apply(final.I, 1, FUN=function(x) paste(x[1:2], collapse="-"))
sort(table(new))

### Create a score for each cluster pattern
in1=paste(c("EBSeq/PostProb_Matrix", my.genotypes[1], my.tissues[1], "txt"), collapse=".")
in2=paste(c("EBSeq/PostProb_Matrix", my.genotypes[2], my.tissues[1], "txt"), collapse=".")
in3=paste(c("EBSeq/PostProb_Matrix", my.genotypes[1], my.tissues[2], "txt"), collapse=".")
in4=paste(c("EBSeq/PostProb_Matrix", my.genotypes[2], my.tissues[2], "txt"), collapse=".")
in5=paste(c("EBSeq/PostProb_Matrix", my.genotypes[1], my.tissues[3], "txt"), collapse=".")
in6=paste(c("EBSeq/PostProb_Matrix", my.genotypes[2], my.tissues[3], "txt"), collapse=".")

prI=read.table(in1, stringsAsFactors=FALSE, header=TRUE)
rioI=read.table(in2, stringsAsFactors=FALSE, header=TRUE)
prL=read.table(in3, stringsAsFactors=FALSE, header=TRUE)
rioL=read.table(in4, stringsAsFactors=FALSE, header=TRUE)
prM=read.table(in5, stringsAsFactors=FALSE, header=TRUE)
rioM=read.table(in6, stringsAsFactors=FALSE, header=TRUE)

combine.probs <- function(x,y) {
    m=merge(x,y, by.x=0, by.y=0)
    rownames(m)=m[,1]
    m[,1]=NULL
    x=x[match(rownames(m), rownames(x)),]
    y=y[match(rownames(m), rownames(y)),]
    xc=as.vector(colnames(x), mode="character")
    yc=as.vector(colnames(y), mode="character")
    z=matrix(0, nrow=nrow(x), ncol=(ncol(x) * ncol(y)))
    rownames(z)=rownames(x)
    patterns=c()
    index=1
    for (i in 1:ncol(x)) {
        for (j in 1:ncol(y)) {
            full=paste(c(xc[i], yc[j]), collapse=".")
            patterns=c(patterns, full)
            z[,index]=x[,i]*y[,j]
            index=index+1
        }
    }
    colnames(z)=patterns
    return(z)
}
comb.I=combine.probs(prI, rioI)
comb.L=combine.probs(prL, rioL)
comb.M=combine.probs(prM, rioM)

gtI=final.I$GbyT_pValue[match(rownames(comb.I), rownames(final.I))]
gtL=final.L$GbyT_pValue[match(rownames(comb.L), rownames(final.L))]
gtM=final.M$GbyT_pValue[match(rownames(comb.M), rownames(final.M))]

prob.I=apply(comb.I, 2, FUN=function(x,y) x*y, y=((-1)*log(gtI)))
prob.L=apply(comb.L, 2, FUN=function(x,y) x*y, y=((-1)*log(gtL)))
prob.M=apply(comb.M, 2, FUN=function(x,y) x*y, y=((-1)*log(gtM)))

clust.scoresI=apply(prob.I,2,sum, na.rm=TRUE)
clust.scoresL=apply(prob.L,2,sum, na.rm=TRUE)
clust.scoresM=apply(prob.M,2,sum, na.rm=TRUE)

top25.I=clust.scoresI[which(clust.scoresI>=15)]
top25.L=clust.scoresL[which(clust.scoresL>=0.3)]
top25.M=clust.scoresM[which(clust.scoresM>=0.03)]
            
