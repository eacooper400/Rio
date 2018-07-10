### RNA-seq Analysis functions
library(DESeq2)
library(RColorBrewer)

calc.poly <- function(x, a) {
    y=0
    for (i in 1:length(a)) {
        tmp=a[i] * (x**(i-1))
        y=y+tmp
    }
    return(y)
}
calculateRPK <- function(counts,size) {
    rpk=counts/size
    return(rpk)
}
cluster.plot <- function(co1, plot.color, raw, center) {
    times=get.times(raw)
    plot(times, raw[1,], type="n", pch=20, col=plot.color, ylim=c(0,10), xlab="", ylab="", main="")
    for (j in 1:nrow(raw)) {
        y=time.medians(raw[j,],times)        
        points(unique(times), y, type="l", col=plot.color)
    }
    int=mean(co1[,1])
    y3=vapply(1:ncol(co1), FUN=calc.poly, FUN.VALUE=double(1), USE.NAMES=FALSE, a=c(int,center))
    lines(x=1:ncol(co1), y3, lwd=3)
}
dge.anyTime <- function(dds, times) {
    dds=DESeq(dds)
    timeTable=get.timeContrasts(dds=dds, times=times)
    best.p=apply(timeTable, 1, min, na.rm=TRUE)
    sigTable=timeTable[which(best.p<0.05),]
    return(sigTable)
}
dge.interaction <- function(dds, genes) {
    dds=DESeq(dds, test="LRT", reduced= ~Genotype + Time)
    m=match(genes, rownames(dds))
    interactions=results(dds)$padj
    sub.int=interactions[m]
    return(sub.int)
}
ebseq.clusterPlot <- function(counts, times, plot.color, name) {
     ymin=min(counts)
     ymax=max(counts)
     plot(times, times, xlab="Time", ylab="Expression", main=name, type="n", ylim=c(ymin,ymax))
     for (j in 1:nrow(counts)) {
         points(times, counts[j,], pch=20, col=plot.color)
         mids=time.medians(counts=counts[j,], times=times)
         lines(unique(times), mids, col=plot.color)
     }
}
get.betas <- function(counts) {
    my.fit=poly.fit(counts)
    my.betas=lapply(my.fit, FUN=coef)
    return(my.betas)
}
get.contrastList <- function(names) {
    max=length(names)
    starts=3:max
    y=(max/2)-1
    ends=vapply(starts, FUN=function(x,y) x+y, FUN.VALUE=double(1), y=y, USE.NAMES=FALSE)
    s1=starts[which(ends<=max)]
    e1=ends[which(ends<=max)]
    xy.list=vector("list", length(s1))
    for (i in 1:length(s1)) {
        xy.list[[i]]=c(s1[i], e1[i])
    }
    return(xy.list)
}
get.timeContrasts <- function(dds, times=my.times) {    
    res=results(dds, name="Genotype_Rio_vs_PR22")
    veg=res$padj
    all.names=resultsNames(dds)
    con=get.contrastList(names=all.names)
    con.names=lapply(con, function(x) list(all.names[x[1]], all.names[x[2]]))
    temp=lapply(con.names, FUN=run.contrast, dds=dds)
    tab=data.frame(matrix(c(veg,unlist(temp)), nrow=nrow(res), byrow=F))
    rownames(tab)=rownames(dds)
    colnames(tab)=times[1:ncol(tab)]
    return(tab)
}
get.times <- function(counts) {
    info=lapply(colnames(counts), FUN=function(x) unlist(strsplit(x, split="\\.")))
    times=unlist(lapply(info, `[[`, 2))
    times=gsub("V", 1, gsub("RI", 2, gsub("FL", 3, gsub("ANT", 4, gsub("SD", 5, times)))))
    times=as.numeric(times)
    return(times)
}
pattern.category <- function(ptable) {
    temp=ptable
    ptable[ptable>0.05]=0
    ptable[ptable>0]=1
    temp$category=apply(ptable,1,function(x) paste0(x,collapse=''))
    return(temp)
}
pattern.filter <- function(ptable) {
    Cat=pattern.category(ptable)
    to.remove=unique(c(grep("NA", Cat$category), grep("^(0|1)00*", Cat$category)))
    new=Cat[-(to.remove),]
    return(new)
}
poly.fit <- function(counts) {
    times=get.times(counts)
    n=(length(unique(times)) - 1)
    my.fit=t(apply(counts, 1, FUN=function(x, times, n) lm(x~poly(times,n, raw=TRUE)), times=times, n=n))
    return(my.fit)
}
run.contrast <- function(contrast.list, dds) {
    result=results(dds, contrast=contrast.list)
    return(result$padj)
}
score.clusters <- function(cluster.list) {
    cmatrix=matrix(0, nrow=length(cluster.list), ncol=length(cluster.list))
    rownames(cmatrix)=cluster.list
    colnames(cmatrix)=cluster.list
    for (i in 1:length(cluster.list)) {
        c1=unlist(strsplit(cluster.list[i], split="-"))
        c1=as.numeric(gsub("Down", 2, gsub("EE", 1, gsub("Up", 0, c1))))
        for (j in 1:length(cluster.list)) {
             c2=unlist(strsplit(cluster.list[j], split="-"))
             c2=as.numeric(gsub("Down", 2, gsub("EE", 1, gsub("Up", 0, c2))))
             cmatrix[i,j]=sum(abs(c2-c1))
        }
    }
    return(cmatrix)
}
setOrigin <- function(counts) {
    times=get.times(counts)
    veg.median=apply(counts[,which(times==1)], 1, median)
    new=apply(counts, 2, function(x,y) x-y, y=veg.median)
    return(new)
}
standardize01 <- function(counts) {
    xmax=apply(counts, 1, max)
    xmin=apply(counts, 1, min)
    scaleby=(xmax-xmin)/10
    newx=apply(counts, 2, function(x, xmin) x-xmin, xmin=xmin)
    newx=apply(newx, 2, function(x, scale) x/scale, scale=scaleby)
    return(newx)
}
time.medians <- function(counts,times) {
    raw=unname(unlist(counts))
    z=unique(times)
    y=c()
    for (i in 1:length(z)) {
        y=c(y, median(raw[which(times==z[i])]))
    }
    return(y)
}
variance.stabilize <- function(counts) {
    gmeans=apply(counts, 1, mean)
    gsd=apply(counts, 1, sd)
    gvar=gsd**2
    fit=lm(log(gvar)~log(gmeans))
    q=fit$coef[2]
    new=counts**(1-(q/2))
}
wssplot <- function(data, nc=15, seed=1234) {
    wss <- (nrow(data)-1)*sum(apply(data,2,var))
    for (i in 2:nc) {
        set.seed(seed)
        wss[i] <- sum(kmeans(data, centers=i, iter.max=100)$withinss)
    }
    plot(1:nc, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
}











