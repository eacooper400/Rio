setwd("~/Desktop/Rio_Analysis/New_Mummer_032018/snps/")
options(stringsAsFactors=F)

snpEff = read.table("snpEff_genes.txt", header=T)
temp = snpEff[,c(1,6,7,5)]
colnames(temp) = c("GeneID", "LOW", "MODERATE", "HIGH")

int = read.table("../../Rio_RNAseq/EBSeq/Sig_GxT_Genes.Pvalues.Internode.txt", header=T)
int = int[which(int$GT.pvalue<0.001),]

leaf = read.table("../../Rio_RNAseq/EBSeq/Sig_GxT_Genes.Pvalues.Leaf.txt", header=T)
leaf = leaf[which(leaf$GT.pvalue<0.001),]

mer = read.table("../../Rio_RNAseq/EBSeq/Sig_GxT_Genes.Pvalues.Meristem.txt", header=T)
mer = mer[which(mer$GT.pvalue<0.001),]

pan = read.table("../../Transcriptome_Blat/Blat_Comp_wBtxGenes_andRILgeno.txt", header=T, sep="\t")

m1 = match(rownames(int), pan$GeneID)
m2 = match(rownames(leaf), pan$GeneID)
m3 = match(rownames(mer), pan$GeneID)

g1 = pan$RIL.genotype[m1]
g2 = pan$RIL.genotype[m2]
g3 = pan$RIL.genotype[m3]

b1 = pan$Btx623_GeneID[m1]
b2 = pan$Btx623_GeneID[m2]
b3 = pan$Btx623_GeneID[m3]

int = int[which(g1=="B"),]
b1 = b1[which(g1=="B")]

leaf = leaf[which(g2=="B"),]
b2 = b2[which(g2=="B")]

mer = mer[which(g3=="B"),]
b3 = b3[which(g3=="B")]

int.snp = temp[match(b1, temp$GeneID),]
leaf.snp = temp[match(b2, temp$GeneID),]
mer.snp = temp[match(b3, temp$GeneID),]

t1 = which(is.na(int.snp$GeneID))
t2 = which(is.na(leaf.snp$GeneID))
t3 = which(is.na(mer.snp$GeneID))

int.snp = int.snp[-t1,]
leaf.snp = leaf.snp[-t2,]
mer.snp = mer.snp[-t3,]

int = int[-t1,]
leaf = leaf[-t2,]
mer = mer[-t3,]

int.snp = int.snp[order(int$GT.pvalue),]
leaf.snp = leaf.snp[order(leaf$GT.pvalue),]
mer.snp = mer.snp[order(mer$GT.pvalue),]

x1 = apply(int.snp[,3:4], 1, sum)
x2 = apply(leaf.snp[,3:4], 1, sum)
x3 = apply(mer.snp[,3:4], 1, sum)

int.snp = int.snp[which(x1>0),]
leaf.snp = leaf.snp[which(x2>0),]
mer.snp = mer.snp[which(x3>0),]

int = int[which(x1>0),]
leaf = leaf[which(x1>0),]
mer = mer[which(x1>0),]

library(RColorBrewer)
sequential <- brewer.pal(3, "BuGn")
par(mar=c(6.1,4.1,4.1,2.1))
barplot(t(mer.snp[1:50,2:4]), names.arg=mer.snp$GeneID[1:50], col=sequential, border="grey", las=2, cex.names=0.65, ylab="SNP Counts")
legend("topright", legend=colnames(temp)[2:4], fill=sequential, title="SNP Impact", border="grey", bty="n")

barplot(t(int.snp[1:50,2:4]), names.arg=int.snp$GeneID[1:50], col=sequential, border="grey", las=2, cex.names=0.65, ylab="SNP Counts")
legend("topright", legend=colnames(temp)[2:4], fill=sequential, title="SNP Impact", border="grey", bty="n")

barplot(t(leaf.snp[1:50,2:4]), names.arg=leaf.snp$GeneID[1:50], col=sequential, border="grey", las=2, cex.names=0.65, ylab="SNP Counts")
legend("topright", legend=colnames(temp)[2:4], fill=sequential, title="SNP Impact", border="grey", bty="n")

new = rbind(int.snp, leaf.snp, mer.snp)
