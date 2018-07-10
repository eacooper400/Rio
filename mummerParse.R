setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Rio_Analysis/Mummer/")

mummer=read.table("../../../Rio_BTx2_diff.1coords", header=FALSE, stringsAsFactors=FALSE)
colnames(mummer)=c("S1", "E1", "S2", "E2", "L1", "L2", "IDY", "G1", "G2", "Cov1", "Cov2", "Chr1", "Chr2")
mummer$Chr2=gsub("chromosome_", "", mummer$Chr2)
mummer$Chr2=as.numeric(mummer$Chr2)
mummer=mummer[complete.cases(mummer),]

btx.names=unlist(lapply(strsplit(btx$V9, split=";"), `[[`, 2))
btx.names=gsub("Name=", "", btx.names)

btx.chr=gsub("Chr", "", btx$V1)
btx.chr=gsub("super_", "1", btx.chr)
btx.chr=gsub("^0", "", btx.chr)
btx.chr=as.numeric(btx.chr)
btx.df=data.frame(btx.names, btx.chr, btx$V4, btx$V5, stringsAsFactors=FALSE)

get.genes=function(chr,start, end, df) {
    x=subset(df, df[,2]==chr)
    y=subset(x, x[,3]>=start)
    z=subset(y, y[,4]<=end)
    return(z[,1])
}       

btx.genes=vector("list", nrow(mummer))
for (i in 1:nrow(mummer)) {
    glist=get.genes(chr=mummer$Chr1[i], start=mummer$S1[i], end=mummer$E1[i], df=btx.df)
    btx.genes[[i]]=glist
}

rio.genes=vector("list", nrow(mummer))
for (i in 1:nrow(mummer)) {
    glist=get.genes(chr=mummer$Chr1[i], start=mummer$S1[i], end=mummer$E1[i], df=blat[,c(9,3,4,5)])
    rio.genes[[i]]=glist
}


### Find Mummer "Block" for every gene:
blocks=rep(NA, nrow(chr3))

for (i in 1:nrow(chr3)) {
    start = chr3$Rio.Start[i]
    if (is.na(start)) {
        next
    }
    p=mumm3[which(mumm3$V4>=start),]
    p2=which(p$V3 <= start)
    if (length(p2) > 0) {
        mbl = p[p2,c(12,1,2)]
        temp.blocks=c()
        for (j in 1:nrow(mbl)) {
            temp.blocks=c(temp.blocks, paste0(mbl[j,], collapse=":"))
        }
        blocks[i]=paste0(temp.blocks, collapse=",")
    }
}
