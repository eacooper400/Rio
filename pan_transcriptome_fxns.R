### Assign a matrix of Blat scores to every Mummer block
blatScoresInMumms = function(Q, T, blat) {
    n = length(Q)
    scoreBlocks = vector("list", n)

    for (i in 1:n) {
        Qblock=Q[[i]]
        Tblock=T[[i]]
        if (length(Qblock) == 0) {
            scoreBlocks[[i]] = as.vector(Tblock, mode="character")
        } else {
            if (length(Tblock) == 0) {
                scoreBlocks[[i]] = as.vector(Qblock, mode="character")
            } else {
                scoreBlocks[[i]] = scoreMatrix(Qblock, Tblock, blat)
            }
        }
    }
    return(scoreBlocks)
}
### Collapse duplicate entries in adjacent blocks (after alignment)
collapse.dups = function(alignTable) {
    temp=alignTable[,10:12]
    x=alignTable[!(duplicated(temp)),]
    altRem=findLooseDups(x[,10:12])
    x=x[-(altRem),]
    return(x)
}
### Find approximate duplicates with NA or missing values
findLooseDups = function(df) {
    to.remove=c()
    dupIDs1=names(table(df[,1]))[which(table(df[,1])>1)]
    dupIDs1=dupIDs1[-1]

    for (i in 1:length(dupIDs1)) {
        my.rows=grep(dupIDs1[i], df[,1])
        y=df[my.rows,]
        check1=grep("-", y[,2])
        if (length(check1)>0) { to.remove=append(to.remove, my.rows[check1]) }
        if ((nrow(y) - length(check1)) > 1) {
            check2 = which(y[,3]==0)
            if (length(check2)>0) { to.remove=append(to.remove, my.rows[check2]) }
        }
    }

    dupIDs2=names(table(df[,2]))[which(table(df[,2])>1)]
    dupIDs2=dupIDs2[-1]

    for (i in 1:length(dupIDs2)) {
        my.rows=grep(dupIDs2[i], df[,2])
        y=df[my.rows,]
        check1=grep("-", y[,1])
        if (length(check1)>0) { to.remove=append(to.remove, my.rows[check1]) }
        if ((nrow(y) - length(check1)) > 1) {
            check2 = which(y[,3]==0)
            if (length(check2)>0) { to.remove=append(to.remove, my.rows[check2]) }
        }
    }
    return(to.remove)
}
### Find which genes have changed position in one genome with respect to another
findPositionShifts = function(alignment, name1, name2) {
    mutations=rep(NA, nrow(alignment))
    start1=grep(paste0(c("GeneStart", name1), collapse="_"), colnames(alignment))
    start2=grep(paste0(c("GeneStart", name2), collapse="_"), colnames(alignment))
    end1=start1+1
    end2=start2+1
    prev1=max(alignment[1,start1], alignment[1,end1])
    prev2=max(alignment[1,start2], alignment[1,end2])
              
    for (i in 1:nrow(alignment)) {
        if (is.na(alignment[i,start1])) {
            mutations[i]="Deletion"
        } else {
            if (is.na(alignment[i,start2])) {
                mutations[i]="Insertion"
            } else {
                if (alignment[i,2] != alignment[i,3]) {
                    mutations[i]="Relocation"
                } else {
                    d1=abs(alignment[i,start1]-prev1)
                    d2=abs(alignment[i,start2]-prev2)
                    prev1=alignment[i,end1]
                    prev2=alignment[i,end2]
                    if (abs(d2-d1)>10000) {
                        mutations[i]="Relocation"
                    } else {
                        mutations[i]="Colinear"
                    }
                }
            }
        }
    }
    return(mutations)
}               
### Convert any non-numeric chromosome characters
fix.chromosomes = function(x) {
    x = gsub("chromosome_", "", x)
    x = gsub("scaffold_", "1", x)
    x = gsub("Chr", "", x)
    x = gsub("super_", "1", x)
    x = as.numeric(x)
    return(x)
}

### Create an F matrix of min possible alignment scores with gap penalty of 1
Fmatrix = function(S) {
    F=matrix(0, nrow=(nrow(S)+1), ncol=(ncol(S)+1))
    F[1,-1] = (-2) * seq(1:ncol(S))
    F[-1,1] = (-2) * seq(1:nrow(S))
    for (i in 2:nrow(F)) {
        for (j in 2:ncol(F)) {
            Match = F[i-1, j-1] + S[i-1, j-1]
            Delete = F[i-1, j] - 2
            Insert = F[i, j-1] - 2
            F[i,j] = max(Match, Insert, Delete)
        }
    }
    return(F)
}

### Find all Genes in Each Block from a Mummer File
geneBlocks=function(mummer, gff, sampleName) {
    blocks=vector("list", nrow(mummer))
    mumm.sample = mummer[, grep(sampleName, colnames(mummer))]
    
    for (i in 1:nrow(mumm.sample)) {
        p1=gff[which(gff$Chr==mumm.sample[i,4]),2:3]
        genes=gff$Name[which(gff$Chr==mumm.sample[i,4])]
        y=unname(unlist(mumm.sample[i,1:2]))
        o=apply(p1, 1, function(x,y) max(x[1],y[1]) <= min(x[2],y[2]), y=sort(y))
        g=which(o)
        if (y[1] < y[2]) {
            blocks[[i]]=as.vector(genes[g], mode="character")
        } else {
            blocks[[i]]=as.vector(genes[sort(g, decreasing=TRUE)], mode="character")
        }
    }
    return(blocks)
}
### Generate a Pairwise Alignment Table for  Mummer Block Score matrix
mummerAlignmentTable = function(Score.matrix, coords) {
    pairwiseAln=data.frame()
    ids=c("Rio", "Sobic")
    if (is.matrix(Score.matrix)) {
        a=traceback(Score.matrix)
        x=matrix(unlist(rep(unname(coords), nrow(a))), nrow=nrow(a), byrow=TRUE)
        colnames(x)=names(coords)
        pairwiseAln = cbind(x,a)
    } else {
        findID = unlist(lapply(ids, function(x,y) length(grep(x, y)), y=Score.matrix))
        tmp = data.frame("AlignmentA"=rep("-",max(findID)), "AlignmentB"=rep("-", max(findID)), "IDY"=rep(NA,max(findID)))
        tmp[,which(findID>0)] = Score.matrix
        x=matrix(unlist(rep(unname(coords), nrow(tmp))), nrow=nrow(tmp), byrow=TRUE)
        colnames(x)=names(coords)
        pairwiseAln = cbind(x,tmp)       
    }
    return(pairwiseAln)
}      

### Parse Blat Results
parse.blat = function(file) {
    blat = read.table(file, header=FALSE, skip=5, sep="\t", stringsAsFactors=FALSE)
    colnames(blat)=c("match", "mismatch", "repMatch", "Ns", "QgapCt", "QgapBases", "TgapCt", "TgapBases", "strand", "Qname", "Qsize", "Qstart", "Qend", "Tname", "Tsize", "Tstart", "Tend", "blockCt", "blockSize", "qStarts", "tStarts")
    blat$Qcoverage=(blat$match + blat$mismatch)/blat$Qsize
    blat$Tcoverage=(blat$match + blat$mismatch)/blat$Tsize
    blat$identity=blat$match/(blat$match + blat$mismatch)
    short = blat[, c(10:11,14:15,22:24)]
    return(short)
}
    

### Parse GFF3 annotation into Gene only table
parse.gff=function(file) {
    temp=read.table(file, header=FALSE, stringsAsFactors=FALSE, sep="\t")
    temp$V1=fix.chromosomes(temp$V1)
    temp=temp[grep("^gene$", temp$V3),]
    x=lapply(temp$V9, function(x) unlist(strsplit(x, split=";")))
    Name=gsub("Name=", "", unlist(lapply(x, `[[`, 2)))
    Name=gsub("\\.g$", "", Name)
    final=data.frame(Chr=temp$V1, Start=temp$V4, End=temp$V5, Name=Name)
    return(final)
}

### Parse Mummer File into relevant table
parse.mummer=function(file) {
    df = read.table(file, header=FALSE, stringsAsFactors=FALSE)
    df[,12] = fix.chromosomes(df[,12])
    df[,13] = fix.chromosomes(df[,13])
    colnames(df) = c("S1", "E1", "S2", "E2", "L1", "L2", "IDY", "G1", "G2", "Cov1", "Cov2", "Chr1", "Chr2")
    small = df[, c(1:6, 12:13)]
    return(small)    
}

### Get a matrix of Blat scores for 2 lists of genes
scoreMatrix = function(q, t, blat) {
    m = matrix(0, nrow=length(q), ncol=length(t))
    rownames(m) = q
    colnames(m) = t

    for (i in 1:length(q)) {
        x = blat[grep(q[i], blat$Qname),]
        for (k in 1:length(t)) {
            y=x[grep(t[k], x$Tname),]
            m[i,k] = max(y$identity)
        }
    }
    m[is.infinite(m)] = 0
    return(m)
}
### Perform traceback to get the best path in F matrix
traceback = function(S) {
    AlignmentA = c()
    AlignmentB = c()
    IDY = c()
    A = rownames(S)
    B = colnames(S)
    F = Fmatrix(S)
    i = length(A) + 1
    j = length(B) + 1

    while ((i > 1) & (j > 1)) {
        Score = F[i,j]
        ScoreDiag = F[i-1, j-1]
        ScoreUp = F[i, j-1]
        ScoreLeft = F[i-1, j]
        if (Score == (ScoreDiag + S[i-1, j-1])) {
            AlignmentA = c(A[i-1], AlignmentA)
            AlignmentB = c(B[j-1], AlignmentB)
            IDY = c(S[i-1, j-1], IDY)
            i = i-1
            j = j-1
            } else {
                if (Score == (ScoreLeft - 2)) {
                    AlignmentA = c(A[i-1], AlignmentA)
                    AlignmentB = c("-", AlignmentB)
                    IDY = c(NA, IDY)
                    i = i-1
                } else {
                    AlignmentA = c("-", AlignmentA)
                    AlignmentB = c(B[j-1], AlignmentB)
                    IDY=c(NA, IDY)
                    j = j-1
                }
            }
    }
    while (i > 1) {
        AlignmentA = c(A[i-1], AlignmentA)
        AlignmentB = c("-", AlignmentB)
        IDY = c(NA, IDY)
        i = i-1
    }
    while (j > 1) {
        AlignmentB = c(B[j-1], AlignmentB)
        AlignmentA = c("-", AlignmentA)
        IDY = c(NA, IDY)
        j = j-1
    }
    result=data.frame(AlignmentA, AlignmentB, IDY)
    return(result)
}

