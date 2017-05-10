# TODO: Add comment
# 
# Author: mg14
###############################################################################

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)

source("PCAWG-functions.R")

vcfFileIn <- args[1]
vcfFileOut <- args[2]

print(vcfFileIn)
print(vcfFileOut)

library(VariantAnnotation)
library(Matrix)

s <- strsplit(vcfFileIn,"/")[[1]]
ID <- sub("\\..+", "", s[length(s)])


print(ID)

#' ## CLUSTERS
# Load clusters
clusters <- try(consensusClustersToOld(loadConsensusClusters((ID))))
NO_CLUSTER <- FALSE
if(class(clusters)=='try-error'){
	clusters <- data.frame(cluster=1, n_ssms=1000, proportion=1)
	NO_CLUSTER <- TRUE
}

#' ### 1. Remove spurious superclonal clusters with less than 10% mutations and (max - 10% VAF)
#clusters <- removeSuperclones(clusters)
#clusters <- mergeClusters(clusters, deltaFreq=0.05)
#
#if(all(is.na(purityPloidy[ID,]))) # Missing purity
#	purityPloidy[ID,] <- c(max(clusters$proportion),NA)
purity <- purityPloidy[ID,'purity']

#' ## COPYNUMBER
bb <- loadConsensusCNA(ID, purity=purityPloidy[ID, 'purity'])

if(NO_CLUSTER)
	clusters <- clustersFromBB(bb)

#' ### Mismatch in CN purity and clusters, use DP purity
#if(purity == 1 | abs(max(clusters$proportion) - purity) > 0.05){
#	bb$clonal_frequency <- bb$clonal_frequency * max(clusters$proportion) / purity
#	purity <- max(clusters$proportion)
#	purityPloidy[ID,'purity'] <- purity
#}


#' ### Missing ploidy - recalculate from BB
if(is.na(purityPloidy[ID,"ploidy"]))
	purityPloidy[ID,"ploidy"] <- averagePloidy(bb = bb)

IS_WGD <- classWgd(bb)

#' ## VCF 
#' Load vcf
vcf <- readVcf(vcfFileIn, genome="GRCh37") #, param=ScanVcfParam(which=pos))

# Add ID & gender
meta(header(vcf))$META <- rbind(meta(header(vcf))$META, DataFrame(Value=c(ID, as.character(allGender[ID, "pred_gender"])), row.names=c("ID", "gender")))

# Add driver genes
vcf <- addFinalDriver(vcf, finalDrivers)

# Add TNC
if(!"TNC" %in% rownames(header(vcf)@header$INFO)){
    tnc=scanFa(file=refFile, resize(granges(vcf), 3,fix="center"))
    i = header(vcf)@header$INFO
    exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="String",Description="Trinucleotide context", row.names="TNC"))
    info(vcf)$TNC <- as.character(tnc)
}

#' Add mutation copy numbers
# vcf <-  addMutCn(vcf, bb, clusters)
i = header(vcf)@header$INFO
exptData(vcf)$header@header$INFO <- rbind(i,mcnHeader())
L <- computeMutCn(vcf, bb, clusters=clusters, purity=purity, xmin=3, gender=as.character(allGender[ID, "pred_gender"]), isWgd=IS_WGD, n.boot=500)
info(vcf) <- cbind(info(vcf), L$D)
bb$timing_param <- L$P 


#' Classify mutations
cls <- classifyMutations(vcf, reclassify='all')

info(vcf)$CLS <- cls
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number="1",Type="String",Description="Mutation classification: {clonal [early/late/NA], subclonal}", row.names="CLS"))

#' Timing
bb$n.snv_mnv <- countOverlaps(bb, vcf)
t <- bbToTime(bb)	
mcols(bb) <- DataFrame(mcols(bb), t)

#' Save output
writeVcf(vcf, file=vcfFileOut)
bgzip(vcfFileOut, overwrite=TRUE)
unlink(vcfFileOut)
save(vcf, file=paste0(vcfFileOut,".RData"))


#' ## INDEL
vcfIndelFileIn <- sub("20160830","20161006",gsub("snv_mnv","indel", vcfFileIn))
vcfIndelFileOut <-  sub("20160830","20161006",gsub("snv_mnv","indel", vcfFileOut))
#' Load vcf
vcfIndel <- readVcf(vcfIndelFileIn, genome="GRCh37") #, param=ScanVcfParam(which=pos))
#' Add ID & gender
meta(header(vcfIndel))$META <- rbind(meta(header(vcfIndel))$META, DataFrame(Value=c(ID, as.character(allGender[ID, "pred_gender"])), row.names=c("ID", "gender")))
#' Add driver genes
vcfIndel <- addFinalDriver(vcfIndel, finalDrivers)
#' Add mutation copy numbers
i = header(vcfIndel)@header$INFO
exptData(vcfIndel)$header@header$INFO <- rbind(i,mcnHeader())
L <- computeMutCn(vcfIndel, bb, clusters=clusters, purity=purity, xmin=3, gender=as.character(allGender[ID, "pred_gender"]), isWgd=IS_WGD)
info(vcfIndel) <- cbind(info(vcfIndel), L$D)
#' Classify mutation
clsIndel <- classifyMutations(vcfIndel, reclassify='all')
info(vcfIndel)$CLS <- clsIndel
info(header(vcfIndel)) <- rbind(info(header(vcfIndel)), DataFrame(Number="1",Type="String",Description="Mutation classification: {clonal [early/late/NA], subclonal}", row.names="CLS"))
#' Save
writeVcf(vcfIndel, file=vcfIndelFileOut)
bgzip(vcfIndelFileOut, overwrite=TRUE)
unlink(vcfIndelFileOut)
save(vcfIndel, file=paste0(vcfIndelFileOut,".RData"))

#' Save BB
bb$n.indel <- countOverlaps(bb, vcfIndel)
seqlevels(bb) <- c(1:22, "X","Y")
bb <- sort(bb)
save(bb, file=sub("snv_mnv","cn",sub(".vcf$",".bb_granges.RData",vcfFileOut)))

#' Save clusters & purity
save(clusters, purity, file=sub("snv_mnv","clusters",sub(".vcf$",".clusters_purity.RData",vcfFileOut)))

#' ## PLOT
pdf(file=sub("snv_mnv","other",sub(".vcf$",".pdf",vcfFileOut)), 12,18)
par(mar=c(1,3,3,1), bty="L", mgp=c(2,.5,0), mfrow=c(4,1),cex=1, las=2)
j <- 1
for(v in c('vcf','vcfIndel')){
	vv <- get(v)
	col <- RColorBrewer::brewer.pal(9, "Set1")[c(3,4,2,1,9)]
	plotVcf(vcf = vv, bb = bb, clusters = clusters, col = col, ID = ID,  IS_WGD = IS_WGD, NO_CLUSTER = NO_CLUSTER, title=j==1)
	j <- j+1
}
plotBB(bb, ylim=c(0,10))
try(plotTiming(bb))
dev.off()
#plot(start(vcf) + w[as.character(seqnames(vcf))], qnorm(info(vcf)$pMutCNTail), col=col[cls], xlab='Position', ylab="pTail", pch=16)


