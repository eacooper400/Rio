setwd("~/Desktop/Rio_Analysis/Rio_RNAseq/")
options(stringsAsFactors=FALSE)
library(msir)

### Read in and get the list of differentially expressed
### genes in each tissue
int = read.table("EBSeq/Sig_GxT_Genes.Pvalues.Internode.txt", header=T)
leaf = read.table("EBSeq/Sig_GxT_Genes.Pvalues.Leaf.txt", header=T)
mer = read.table("EBSeq/Sig_GxT_Genes.Pvalues.Meristem.txt", header=T)

int = int[which(int$GT.pvalue<0.001),]
leaf = leaf[which(leaf$GT.pvalue<0.001),]
mer = mer[which(mer$GT.pvalue<0.001),]

### Read in the GO terms for all genes
### Get the list of DEGs within each relevant
### tissue that falls into each category
go = read.table("geneID2GO.txt", header=F)
go.cm = go[grep("GO:0005975", go$V2),1]
go.it = go[grep("GO:0055085", go$V2),1]
go.pp = go[grep("GO:0006468", go$V2),1]
go.mt = go[grep("GO:000701(7|8)", go$V2),1]
go.pa = go[grep("GO:0016887", go$V2),1]
go.ly = go[grep("GO:0016829", go$V2),1]

int.go.cm = int[which(rownames(int) %in% go.cm),]
leaf.go.it = leaf[which(rownames(leaf) %in% go.it),]
int.go.pp = int[which(rownames(int) %in% go.pp),]
int.go.mt = int[which(rownames(int) %in% go.mt),]
int.go.pa = int[which(rownames(int) %in% go.pa),]
leaf.go.pa = leaf[which(rownames(leaf) %in% go.pa),]
mer.go.pa = mer[which(rownames(mer) %in% go.pa),]
leaf.go.ly = leaf[which(rownames(leaf) %in% go.ly),]

### Read in the pan-transcriptome info,
### and further divide DEGs by genomic
### background in relevant subsets
pan = read.table("../Transcriptome_Blat/Blat_Comp_wBtxGenes_andRILgeno.txt", header=T, sep="\t")
g1 = pan$RIL.genotype[match(rownames(int.go.cm), pan$GeneID)]
g2 = pan$RIL.genotype[match(rownames(leaf.go.it), pan$GeneID)]
g3 = pan$RIL.genotype[match(rownames(int.go.pp), pan$GeneID)]
g4 = pan$RIL.genotype[match(rownames(int.go.mt), pan$GeneID)]
g5 = pan$RIL.genotype[match(rownames(int.go.pa), pan$GeneID)]
g6 = pan$RIL.genotype[match(rownames(leaf.go.pa), pan$GeneID)]
g7 = pan$RIL.genotype[match(rownames(mer.go.pa), pan$GeneID)]
g8 = pan$RIL.genotype[match(rownames(leaf.go.ly), pan$GeneID)]

int.go.cm.R = int.go.cm[which(g1=="R"),]
leaf.go.it.B = leaf.go.it[which(g2=="B"),]
leaf.go.it.R = leaf.go.it[which(g2=="R"),]
int.go.pp.R = int.go.pp[which(g3=="R"),]
int.go.mt.B = int.go.mt[which(g4=="B"),]
int.go.mt.R = int.go.mt[which(g4=="R"),]
int.go.pa.R = int.go.pa[which(g5=="R"),]
leaf.go.pa.B = leaf.go.pa[which(g6=="B"),]
mer.go.pa.B = mer.go.pa
leaf.go.ly.R = leaf.go.ly[which(g8=="R"),]

### Read in the variance stabilized counts for all genes,
### sort them into their respective subsets
vsd = read.table("varianceStabilized.counts.txt", header=T)
vsd = vsd[, -(grep("BL", colnames(vsd)))]

vsd1 = vsd[match(rownames(int.go.cm.R), rownames(vsd)),]
vsd1 = vsd1[, grep("I[1-3]$", colnames(vsd1))]
vsd2 = vsd[match(rownames(leaf.go.it.B), rownames(vsd)),]
vsd2 = vsd2[, grep("L[1-3]$", colnames(vsd2))]
vsd3 = vsd[match(rownames(leaf.go.it.R), rownames(vsd)),]
vsd3 = vsd3[, grep("L[1-3]$", colnames(vsd3))]
vsd4 = vsd[match(rownames(int.go.pp.R), rownames(vsd)),]
vsd4 = vsd4[, grep("I[1-3]$", colnames(vsd4))]
vsd5 = vsd[match(rownames(int.go.mt.B), rownames(vsd)),]
vsd5 = vsd5[, grep("I[1-3]$", colnames(vsd5))]
vsd6 = vsd[match(rownames(int.go.mt.R), rownames(vsd)),]
vsd6 = vsd6[, grep("I[1-3]$", colnames(vsd6))]
vsd7 = vsd[match(rownames(int.go.pa.R), rownames(vsd)),]
vsd7 = vsd7[, grep("I[1-3]$", colnames(vsd7))]
vsd8 = vsd[match(rownames(leaf.go.pa.B), rownames(vsd)),]
vsd8 = vsd8[, grep("L[1-3]$", colnames(vsd8))]
vsd9 = vsd[match(rownames(mer.go.pa.B), rownames(vsd)),]
vsd9 = vsd9[, grep("(M|S|F)[1-3]$", colnames(vsd9))]
vsd10 = vsd[match(rownames(leaf.go.ly.R), rownames(vsd)),]
vsd10 = vsd10[, grep("L[1-3]$", colnames(vsd10))]

### Read in the info. on the day when each sample was collected
brix = read.csv("brix.csv", header=T)
brix$sampleID = gsub("-", ".", brix$sampleID)

### Separately get the fitted lines for each subset (vsd table)
my.subsets = list(vsd1, vsd2, vsd3, vsd4, vsd5, vsd6, vsd7, vsd8, vsd9, vsd10)
for (i in 1:10) {
    d = my.subsets[[i]]
    times = brix$DaysSinceVeg[match(colnames(d), brix$sampleID)]

    dR = d[, grep("Rio", colnames(d))]
    dP = d[, grep("PR22", colnames(d))]
    timesR = times[grep("Rio", colnames(d))]
    timesP = times[grep("PR22", colnames(d))]

    
    plot(c(1,1), c(2,2), type="n", xlim=c(0,72), ylim=c(0,21),
         xlab="Days Since Vegetative", ylab="Normalized Count")

    for (j in 1:nrow(dR)) {
        yR = dR[j,][order(timesR)]
        yP = dP[j,][order(timesP)]
        lines(sort(timesR), yR, col="dodgerblue4")
        lines(sort(timesP), yP, col="darkorange")
        }
    ##xR = rep(timesR, nrow(dR))
    ##xP = rep(timesP, nrow(dP))
    ##yR = as.vector(t(dR))
    ##yP = as.vector(t(dP))
    ##lR = loess.sd(x = xR, y = yR, nsigma=1.96)
    ##plot(lR$x, lR$y, type="l", lty=1, ylim=c(0,20))
    ##lines(lR$x, lR$upper, lty=2)
    ##lines(lR$x, lR$lower, lty=2)
    ##lP = loess.sd(x = xP, y=yP, nsigma=1.96)
    ##lines(lP$x, lP$y, col="red")
    ##lines(lP$x, lP$upper, col="red", lty=2)
    ##lines(lP$x, lP$lower, col="red", lty=2)
