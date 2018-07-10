install.packages("zoo")
library(zoo)

### Find Recombination breaks points in PR22 RIL
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Rio_Analysis/")
source("../Teaching/CompGenWS_Spring17/LC_functions.R")
vcf=my.read.vcf(file="all.vcf", header=TRUE, comment.char='')
###nrow(vcf)
###[1] 1812578

missing=c(grep("^\\.", vcf$BTx3197), grep("^\\.", vcf$PR22), grep("^\\.", vcf$Rio))
length(unique(missing))
###[1] 142735
vcf2=vcf[-unique(missing),]
nrow(vcf2)
###[1] 1669843

vcf=vcf2
rio.poly=c(grep("^0/1", vcf$Rio), grep("^1/1", vcf$Rio))
length(rio.poly)
###[1] 55538
vcf2=vcf[-c(rio.poly),]
nrow(vcf2)
###[1] 1614305

vcf=vcf2
btx.poly=grep("^0/1", vcf$BTx3197)
length(btx.poly)
###[1] 121978
vcf2=vcf[-c(btx.poly),]
nrow(vcf2)
###[1] 1492327
vcf=vcf2

no.info=intersect(grep("^0/0", vcf$BTx3197), grep("^0/0", vcf$Rio))
length(no.info)
###[1] 18066
vcf2=vcf[-c(no.info),]
nrow(vcf2)
###[1] 1474261
vcf=vcf2

new.table=cbind(as.vector(vcf$CHROM, mode="character"),
                as.vector(vcf$POS, mode="character"),
                as.vector(vcf$REF, mode="character"),
                as.vector(vcf$ALT, mode="character"),
                as.vector(vcf$BTx3197, mode="character"),
                as.vector(vcf$PR22, mode="character"),
                as.vector(vcf$Rio, mode="character"))
colnames(new.table)=c("CHROM", "POS", "REF", "ALT", "BTx3197", "PR22", "Rio")

get.genotype=function(x) {
    y=unlist(strsplit(x, split=":"))[1]
    z=unlist(strsplit(y, split="/"))
    return(as.numeric(z[1]) + as.numeric(z[2]))
}

test=vapply(new.table[,5], FUN=get.genotype, FUN.VALUE=double(1), USE.NAMES=FALSE)
BTx3197=test
new.table[,5]=BTx3197

test=vapply(new.table[,6], FUN=get.genotype, FUN.VALUE=double(1), USE.NAMES=FALSE)
PR22=test
new.table[,6]=PR22

test=vapply(new.table[,7], FUN=get.genotype, FUN.VALUE=double(1), USE.NAMES=FALSE)
Rio=test
new.table[,7]=Rio

my.chr=levels(as.factor(new.table[,1]))
my.chr=my.chr[grep("chromosome_", my.chr)]

my.strings=c()
for (c in 1:10) {
    data=new.table[grep(my.chr[c], new.table[,1]),6]
    string=paste(as.vector(data), collapse="")
    my.strings=c(my.strings,string)
}
write.table(new.table, file="Recombinat_SNPs.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

#####################################
### Part 2: go through PR22 in sliding windows,
### calculate proportion of Rio to BTx sites in each 15 SNP window

calcProp <- function(x) {
    A <- length(grep("0", x))
    B <- length(grep("2", x))
    prop=A
    if (B > 0) { prop=(A/B) }
    return(prop)
}

classifyProp <- function(x) {
    if (x > 2) { return("R") }
    else {
        if (x < 0.25) { return("B") }
        else { return("H") }
    }
}

simpleSNPs=read.table("Recombinat_SNPs.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
simpleSNPs=simpleSNPs[-(which(simpleSNPs$PR22==1)),]
my.chr=levels(as.factor(simpleSNPs$CHROM))
clist=c()
plist=c()
codes=c()
for (c in 1:10) {
    data=simpleSNPs[grep(my.chr[c], simpleSNPs$CHROM),]
    window.starts=seq(from=1, to=nrow(data), by=14)
    for (w in 1:(length(window.starts))) {
        v=as.vector(data[window.starts[w]:(window.starts[w]+15),6], mode="numeric")
        prop=calcProp(v)
        code=classifyProp(prop)
        clist=c(clist,my.chr[c])
        plist=c(plist,data[window.starts[w],2])
        codes=c(codes,code)
    }
}
new.table=cbind(clist,plist,codes)
colnames(new.table)=c("CHROM", "POS", "TYPE")

write.table(new.table,"Recomb_SNPs_15win.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

#### Identify break points
### whenever you go from A->B or B->A
my.bps=vector("list", 10)
for (c in 1:10) {
    sub=new.table[grep(my.chr[c], new.table[,1]),]
    y=as.vector(sub[,3], mode="character")
    p=as.vector(sub[,2], mode="numeric")
    ah=p[which(rollapply(y, 2, identical, c("R", "H")))]    
    ab=p[which(rollapply(y, 2, identical, c("R", "B")))]
    ba=p[which(rollapply(y, 2, identical, c("B", "R")))]
    bh=p[which(rollapply(y, 2, identical, c("B", "H")))]
    ha=p[which(rollapply(y, 2, identical, c("H", "R")))]
    hb=p[which(rollapply(y, 2, identical, c("H", "B")))]
    all=c(ah,ab,ba,bh,ha,hb)
    names(all) = c(rep("RH", length(ah)), rep("RB", length(ab)), rep("BR", length(ba)), rep("BH", length(bh)), rep("HR", length(ha)), rep("HB", length(hb)))
    all=sort(all)
    my.bps[[c]] = all
}
    
