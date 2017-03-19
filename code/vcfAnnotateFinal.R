# TODO: Add comment
# 
# Author: mg14
###############################################################################

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)

source("functions.R")

vcfFileIn <- args[1]
vcfFileOut <- args[2]

print(vcfFileIn)
print(vcfFileOut)

library(VariantAnnotation)
library(Matrix)
library(CoxHD)
library(igraph)

refFile = "/lustre/scratch112/sanger/cgppipe/PanCancerReference/genome.fa.gz" #meta(header(v))["reference",]
refLengths <- scanFaIndex(file=refFile)

dpPath <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/icgc_pancan_full/consensus_clustering/201703_consensus_clustering/consensus_clusters_wm"
dpFiles <- dir(dpPath, pattern="_subclonal_structure.txt", recursive=TRUE)

sampleIds <- sub("_mutation_assignments.txt.gz","",dir(dpPath, pattern="_mutation_assignments.txt.gz"))
sampleIds <- intersect(sampleIds, sub("\\..+","",dir(vcfPath, pattern=".bgz$")))

s <- strsplit(vcfFileIn,"/")[[1]]
ID <- sub("\\..+", "", s[length(s)])


print(ID)

#' ## CLUSTERS
# Load clusters
clusters <- try(loadClusters(ID))
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
#purity <- purityPloidy[ID,'purity']

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
vcf <- addFinalDriver(vcf, driVers)

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
L <- computeMutCn(vcf, bb, clusters=clusters, purity=purity, xmin=3, gender=as.character(allGender[ID, "pred_gender"]), isWgd=IS_WGD, n.boot=10)
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
vcfIndel <- addFinalDriver(vcfIndel, driVers)
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
	cls <- factor(paste(as.character(info(vv)$CLS)), levels = c(levels(info(vv)$CLS), "NA"))
	if(j>1) par(mar=c(3,3,1,1))
	plot(start(vv) + chrOffset[as.character(seqnames(vv))], getAltCount(vv)/getTumorDepth(vv),col=col[cls], xlab='', ylab="VAF", pch=ifelse(info(vv)$pMutCNTail < 0.025 | info(vv)$pMutCNTail > 0.975, 4 , 16), ylim=c(0,1), xlim=c(0,chrOffset["MT"]), xaxt="n", cex=.66)
	if(j==1){
		title(main=paste0(ID,", ", donor2type[sample2donor[ID]], "\nploidy=",round(averagePloidy(bb),2), ", hom=",round(averageHom(bb),2), if(IS_WGD) ", WGD" else "", if(NO_CLUSTER) ", (No clusters available)" else(paste0(", clusters=(",paste(round(clusters$proportion, 2), collapse="; "),")"))), font.main=1, line=1, cex.main=1)
	} 
	abline(v = chrOffset[1:25], lty=3)
	mtext(side=1, line=1, at=chrOffset[1:24] + diff(chrOffset[1:25])/2, text=names(chrOffset[1:24]))
	for(i in seq_along(bb)) try({
					s <- start(bb)[i]
					e <- end(bb)[i]
					x <- chrOffset[as.character(seqnames(bb)[i])]
					y <- bb$timing_param[[i]][,"f"]
					l <- bb$timing_param[[i]][,"pi.s"] * bb$timing_param[[i]][,"P.m.sX"]
					segments(s+x,y,e+x,y, lwd=l*4+.1)
					#text(x=(s+e)/2 +x, y=y, paste(signif(bb$timing_param[[i]][,"m"],2),signif(bb$timing_param[[i]][,"cfi"]/purityPloidy[meta(header(vv))["ID",1],"purity"],2), sep=":"), pos=3, cex=0.5)
				}, silent=TRUE)
	legend("topleft", pch=19, col=col, legend=paste(as.numeric(table(cls)), levels(cls)), bg='white')
	j <- j+1
}
plotBB(bb, ylim=c(0,10))
plotTiming(bb)
dev.off()
#plot(start(vcf) + w[as.character(seqnames(vcf))], qnorm(info(vcf)$pMutCNTail), col=col[cls], xlab='Position', ylab="pTail", pch=16)


