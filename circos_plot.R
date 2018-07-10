setwd("~/Desktop/Rio_Analysis/Circos_Data/")
library(circlize)

tiff("~/Desktop/Manuscripts/Rio/Figures/Figure_Components/circos_highRes.png", width=7, height=7, units="in", res=300)

### Read in the chromosome lengths
### and the mummmer 1-to-1 alignments
chr=read.csv("chromosomes.csv", header=TRUE)
coords=read.csv("colinear_rect.csv")

### Initialize the plot and the plotting regions
circos.par("track.height" = 0.15)
circos.initialize(factors=chr$chr, xlim=as.matrix(chr[,2:3]))
circos.trackPlotRegion(ylim=c(0,1), bg.col=rep("white", 10), bg.border=rep("white",10))


### Label the Rio chromosomes
chr.center = as.integer(chr$end/2)
y = rep(1, 10) + uy(5, "mm")
for (i in 1:10) {
    circos.text(x=chr.center[i], y=y[i], paste("Rio", i, sep=" "), sector.index=chr$chr[i])
}

### Add rectangles for the alignment blocks
library(RColorBrewer)
my.colors=brewer.pal(9, "Set1")
my.colors = c(my.colors, "hotpink4")

for (i in 1:nrow(coords)) {
    circos.rect(xleft=coords$start[i], ybottom=0,
                xright=coords$end[i], ytop=1,
                sector.index=coords$chr[i], col=my.colors[coords$group[i]], border=NA)
    }

### On the next track, add color coded rectangles for recombination blocks in the RIL
block.colors = c("dodgerblue", "grey", "darkorange")
bpts = read.csv("ril_bpts.csv", header=TRUE)
bpts = bpts[,1:4]
circos.trackPlotRegion(ylim=c(0,1), bg.col=rep("white", 10), track.height=0.08, bg.border=rep("white", 10))
for (i in 1:nrow(bpts)) {
    circos.rect(xleft=bpts$start[i], ybottom=0,
                xright=bpts$end[i], ytop=1,
                sector.index=bpts$chr[i], col=block.colors[bpts$group[i]], border=NA)
    }


### On the next track, add color coded rectangles for inversions
###x = which(coords$start>coords$end)
###inv = coords[x,]
###circos.trackPlotRegion(ylim=c(0,1), bg.col=rep("white", 10), track.height=0.15, bg.border=rep("white", 10))
###for (i in 1:nrow(inv)) {
###    circos.rect(xleft=inv$end[i], ybottom=0,
###                xright=inv$start[i], ytop=1,
###                sector.index=inv$chr[i], col=my.colors[inv$group[i]], border=NA)
###    }

### Read in Assemblytics bed file,
### Find gene regions that have been tandemly duplicated or deleted
###bed = read.table("../New_Mummer_032018/assemblytics/sv_wGeneInfo.bed", header=T, stringsAsFactors=F, fill=T, sep="\t")
###count.genes = function(x) {
###    y = unlist(strsplit(x, ","))
###    return(length(y))
###}
###bed$BTx_Count = unlist(lapply(bed$BTx_Genes, count.genes))
###bed$Rio_Count = unlist(lapply(bed$Rio_Genes, count.genes))
###bed = bed[order(bed$query, bed$query_start),]
###x = which(bed$BTx_Count != bed$Rio_Count)
###bed2 = bed[x,]

### Add a new track in circos
###circos.trackPlotRegion(ylim=c(-200,200), bg.col=rep("white", 10), track.height=0.25, bg.border=rep("white", 10))

win.starts = c()
win.chr = c()

for (i in 1:10) {
    win.starts = c(win.starts, seq(1, chr$end[i], by=1000000))    
    win.chr = c(win.chr, rep(i, length(seq(1, chr$end[i], by=1000000))))
}
win.ends = win.starts+1000000
new = data.frame(Chr=win.chr, Start=win.starts, End=win.ends, Exp=rep(0, length(win.starts)), Con=rep(0, length(win.ends)))

###for (i in 1:nrow(new)) {
###    temp = bed2[which(bed2$query==new$Chr[i]),]
###    temp2 = subset(temp, temp$query_start>=win.starts[i] & temp$query_stop<win.ends[i])
###    exp = temp2[which(temp2$Rio_Count>temp2$BTx_Count),]
###    new$Exp[i] = sum(exp$size)
###    con = temp2[which(temp2$Rio_Count<temp2$BTx_Count),]
 ###   new$Con[i] = sum(con$size)
###}
###new$Exp = new$Exp/1000
###new$Con = (new$Con/1000) * (-1)
###for (i in 1:nrow(new)) {
###    circos.rect(xleft=new$Start[i], ybottom=0, xright=new$End[i], ytop=new$Exp[i], sector.index=chr$chr[new$Chr[i]], col="seagreen3", border=NA)
###    circos.rect(xleft=new$Start[i], ybottom=new$Con[i], xright=new$End[i], ytop=0, sector.index=chr$chr[new$Chr[i]], col="deeppink2", border=NA)
###}

### Read in the vcf file to get SNP density
vcf = read.table("../Mummer/Rio_BTx_Mummer.vcf", header=F, stringsAsFactors=F)
new$SNPs = rep(0, nrow(new))
for (i in 1:nrow(new)) {
    temp = vcf[which(vcf$V1==new$Chr[i]),]
    temp2 = subset(temp, temp$V2>=win.starts[i] & temp$V2<win.ends[i])
    new$SNPs[i]=nrow(temp2)
}
new$SNPs = new$SNPs/1000000

### Add another track for SNP Density
circos.trackPlotRegion(ylim=c(0,.011), bg.col=rep("white", 10), bg.border=rep("white", 10))
for (i in 1:10) {
    temp = new[which(new$Chr==i),]
    circos.lines(temp$Start, temp$SNPs, sector.index=chr$chr[i], area=TRUE, col="grey")
    }

###for (i in 1:10) {
   ### win.starts = seq(1, chr$end[i], by=1000000)
   ### win.ends = win.starts+1000000
   ### temp = bed[which(bed$query==i),]
   ### for (j in 1:length(win.starts)) {
      ###  temp2 = subset(temp, temp$query_start>=win.starts[j] & temp$query_stop<win.ends[j])
        ###exp = (sum(temp2$size[which(temp2$type=="Repeat_expansion")]))/200000
        ##exp = exp/2
        ##con = (sum(temp2$size[which(temp2$type=="Repeat_contraction")]))/200000
        ##con = con/2
        ##circos.rect(xleft=win.starts[j], ybottom=0.5, xright=win.ends[j], ytop=(exp+0.5), sector.index=chr$chr[i], col="mediumorchid3", border=NA)
        ##circos.rect(xleft=win.starts[j], ybottom=(0.5-con), xright=win.ends[j], ytop=0.5, sector.index=chr$chr[i], col="goldenrod2", border=NA)
    ##}
##}

### Plot the ratio of nonsynonymous to synonymous mutations
snpEff = read.table("../New_Mummer_032018/snps/snpEff_genes.txt", header=T)
temp = data.frame(GeneID = snpEff$GeneId, Syn = snpEff$variants_effect_synonymous_variant, Non=snpEff$variants_effect_missense_variant)
temp$Syn[temp$Syn==0] = 1
temp$Non[temp$Non==0] = 1

temp$GeneID = gsub("\\.v*", "", temp$GeneID)
temp$GeneID = gsub("32$", "", temp$GeneID)
temp$Ratio = temp$Non/temp$Syn

geneInfo = read.table("../Annotation_Files/Sbicolor_313_v3.1.gene.gff3", header=F, stringsAsFactors=F)
geneInfo = geneInfo[grep("gene", geneInfo$V3),]
id = unlist(lapply(strsplit(geneInfo$V9, split= ";"), `[[`, 2))
id = gsub("Name=", "", id)
id = gsub("\\.", "", id)
m = match(temp$GeneID, id)
geneInfo$V1 = as.numeric(gsub("Chr", "", geneInfo$V1))
temp$Chr = geneInfo$V1[m]
temp$Pos = geneInfo$V4[m]
temp = temp[which(!(is.na(temp$Chr))),]

circos.trackPlotRegion(ylim=c(0,21), bg.col=rep("white", 10), bg.border=rep("white", 10))
for (i in 1:10) {
    temp2 = temp[which(temp$Chr==i),]
    circos.lines(temp2$Pos, temp2$Ratio, sector.index=chr$chr[i], area=FALSE)
}















### Read in rdiff coordinates to find inversions
qdiff = read.table("../New_Mummer_032018/assemblytics/rioVbtx_assemblytics.qdiff", header=F, stringsAsFactors=FALSE, fill=TRUE)
q = qdiff[which(qdiff$V2=="INV"),]
q$V1 = as.numeric(gsub("chromosome_", "", q$V1))

circos.trackPlotRegion(ylim=c(0,1), bg.col=rep("white", 10), track.height=0.1, bg.border=rep("white", 10))
for (i in 1:nrow(q)) {
    circos.rect(xleft=q$V3[i], ybottom=0, xright=q$V4[i], ytop=1, sector.index=chr$chr[q$V1[i]], col="darkorange4", border=NA)
}
