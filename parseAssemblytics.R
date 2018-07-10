setwd("~/Desktop/Rio_Analysis/New_Mummer_032018/assemblytics/")
options(stringsAsFactors=FALSE)
source("../../scripts/pan_transcriptome_fxns.R")

### Read in Gene Blocks for unique Mums and SV events
block.table = read.table("GeneBlocks_wAssemblyticsCoords.txt", header=T)
block.table2 = read.table("StrucVar_GeneBlocks.txt", header=T)

### Count relocated duplications as genes where the Rio ID
### appears only once, but the BTx ID appears more than once
### Be sure to exclude "-" matches with no transcript in BTx
reloc.dup = intersect(which(table(block.table$AlignmentA)[block.table$AlignmentA]==1), which(table(block.table$AlignmentB)[block.table$AlignmentB]>1 & table(block.table$AlignmentB)[block.table$AlignmentB]<4))
reloc.dup.rio = unique(block.table$AlignmentA[reloc.dup])
reloc.dup.btx = unique(block.table$AlignmentB[reloc.dup])

### Count relocated deletions as genes where the Rio ID
### appears more than once, but the BTx ID only appears once
reloc.del = intersect(which(table(block.table$AlignmentB)[block.table$AlignmentB]==1), which(table(block.table$AlignmentA)[block.table$AlignmentA]>1 & table(block.table$AlignmentA)[block.table$AlignmentA]<4))
reloc.del.rio = unique(block.table$AlignmentA[reloc.del])
reloc.del.btx = unique(block.table$AlignmentB[reloc.del])

                     
### Read in the transcript file for each
rio.genes=parse.gff("../../Annotation_Files/SbicolorRiov2.1.gene.gff3")
btx.genes = parse.gff("../../Annotation_Files/Sbicolor_454_v3.1.1.gene.gff3")

### Determine which genes are not in any Mummer blocks
rio1 = rio.genes$Name[which(is.na(match(rio.genes$Name, block.table$AlignmentA)))]
btx1 = btx.genes$Name[which(is.na(match(btx.genes$Name, block.table$AlignmentB)))]

### Determine which missing Rio genes can be found in assemblytics blocks
m1 = match(rio1, block.table2$AlignmentA)
r = which(!(is.na(m1)))

### Of these, how many match a missing Btx gene?
b = block.table2$AlignmentB[m1[r]]
m2 = match(btx1, b)

### Add Feature Types to the Mummer Table
block.table$Type = rep("Collinear", nrow(block.table))
reloc = which(block.table$BtxChr != block.table$RioChr)
block.table$Type[reloc] = "Relocation"
block.table$Type[reloc.dup] = "Relocation_wDuplication_Rio"
block.table$Type[reloc.del] = "Relocation_wDeletion_Rio"
block.table$Type[which(block.table$AlignmentA=="-")] = "NoTranscript_Rio"
block.table$Type[which(block.table$AlignmentB=="-")] = "NoTranscript_Btx"
block.table$Type[which(block.table$RioS > block.table$RioE)] = "Inversion"

### Set aside genes in expansions, then merge
### the blocks with missing Rio genes into the block.table
to.add = block.table2[m1[r],]
rio.exp = to.add$AlignmentA[which(b=="-")]
to.add$Type = rep("Collinear_Gap", nrow(to.add))
to.add$Type[which(b=="-")] = "Expansion_Rio"
bdup = seq(1, length(b))[-unique(c(m2[which(!(is.na(m2)))], which(b=="-")))]
to.add$Type[bdup] = "Other_Duplication_Rio"
block.table = rbind(block.table, to.add)

### Subtract the found rio genes and found BTx genes
### from the lists of missing genes
rio1 = rio1[-r]
btx1 = btx1[-(which(!(is.na(m2))))]

### Determine which additional btx genes can be found in Assemblytics blocks
block.table3 = block.table2[-(m1[r]),]
m3 = match(btx1, block.table3$AlignmentB)
b = which(!(is.na(m3)))
block.table4 = block.table3[m3[which(!(is.na(m3)))],]
block.table4$Type = rep("Gap", nrow(block.table4))
block.table4$Type[grep("-", block.table4$AlignmentA)] = "Contraction_Rio"
rio.del = block.table4$AlignmentB[grep("-", block.table4$AlignmentA)]
block.table4$Type[-(grep("-", block.table4$AlignmentA))] = "Other_Duplication_BTx"
btx1 = btx1[-(which(!(is.na(m3))))]
block.table = rbind(block.table, block.table4)

### For the remaining missing genes,
### filter out those that map to a chromosome other than 1-10
rio1 = rio1[-(which(rio.genes$Chr[match(rio1, rio.genes$Name)]>10))]
btx1 = btx1[-(which(btx.genes$Chr[match(btx1, btx.genes$Name)]>10))]

### Run Blat on missing genes against whole genome sequence,
### then read in the results
btx.blat = parse.blat("BtxMissing_v_Rio.psl")
rio.blat = parse.blat("RioMissing_v_Btx.psl")
btx.blat = btx.blat[which(btx.blat$Qcoverage>=0.8),]

rawblat = read.table("RioMissing_v_Btx.psl", header=F, skip=5)
rawblat = rawblat[which(rio.blat$Qcoverage>=0.8),]
rio.blat = rio.blat[which(rio.blat$Qcoverage>=0.8),]


### Sort the alignment table
block.table3 = block.table[order(block.table$RioChr, block.table$RioS),]

rblock.m = data.frame(block=rep(NA, length(rio1)),
                      BtxS=rep(NA, length(rio1)),
                      BtxE=rep(NA, length(rio1)),
                      RioS=rio.genes$Start[m],
                      RioE=rio.genes$End[m],
                      BtxL=rep(NA, length(rio1)),
                      RioL=(rio.genes$End[m]-rio.genes$Start[m]),
                      BtxChr=rep(NA, length(rio1)),
                      RioChr=rio.genes$Chr[m],
                      AlignmentA=rio1,
                      AlignmentB=rep("-", length(rio1)),
                      IDY=rep(NA, length(rio1)),
                      Type=rep("Unique", length(rio1)))

