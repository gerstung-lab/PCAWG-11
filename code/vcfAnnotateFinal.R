# TODO: Add comment
# 
# Author: mg14
###############################################################################


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

dpFiles <- dir(dpPath, pattern="_subclonal_structure.txt", recursive=TRUE)

sampleIds <- sub("_mutation_assignments.txt.gz","",dir(dpPath, pattern="_mutation_assignments.txt.gz"))
sampleIds <- intersect(sampleIds, sub("\\..+","",dir(vcfPath, pattern=".bgz$")))

s <- strsplit(vcfFileIn,"/")[[1]]
ID <- sub("\\..+", "", s[length(s)])


print(ID)

#' ## CLUSTERS
# Load clusters
clusters <- loadClusters(ID)

if(all(is.na(purityPloidy[ID,]))) # Missing purity
	purityPloidy[ID,] <- c(max(clusters$proportion),NA)

purity <- purityPloidy[ID,'purity']

# Fix clusters with proportion > purity
w <- clusters$proportion >= purity - 0.075 #ie ~ 1.5 reads 
cl <- as.data.frame(rbind(if(any(!w)) clusters[!w,,drop=FALSE], if(any(w)) colSums(clusters[w,,drop=FALSE])))
cl[nrow(cl),"proportion"] <- purity
clusters <- cl
#clusters <- mergeClusters(clusters, deltaFreq=0.05)


#' ## COPYNUMBER
bb <- loadBB(ID)

if(purity == 1){
	purity <- max(clusters$proportion)
	bb$clonal_frequency <- bb$clonal_frequency * purity
	purityPloidy[ID,'purity'] <- purity
}
	

# Missing ploidy - recalculate from BB
if(is.na(purityPloidy[ID,"ploidy"]))
	purityPloidy[ID,"ploidy"] <- sum(width(bb) * bb$copy_number * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)

#' ## VCF 
#' Load vcf
vcf <- readVcf(vcfFileIn, genome="GRCh37") #, param=ScanVcfParam(which=pos))

# Add ID & gender
meta(header(vcf)) <- rbind(meta(header(vcf)), DataFrame(Value=c(ID, as.character(allGender[ID, "pred_gender"])), row.names=c("ID", "gender")))

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
L <- computeMutCn(vcf, bb, clusters=clusters, purity=purity, xmin=3, gender=as.character(allGender[ID, "pred_gender"]))
info(vcf) <- cbind(info(vcf), L$D)
bb$timing_param <- L$P 


#' Classify mutations
cls <- classifyMutations(vcf, reclassify='all')

info(vcf)$CLS <- cls
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number="1",Type="String",Description="Mutation classification: {clonal [early/late/NA], subclonal}", row.names="CLS"))

#' Save output
writeVcf(vcf, file=vcfFileOut)
bgzip(vcfFileOut, overwrite=TRUE)
unlink(vcfFileOut)
save(bb, file=sub(".vcf$",".bb_granges.RData",vcfFileOut))
save(vcf, file=paste0(vcfFileOut,".RData"))


#' ## INDEL
vcfIndelFileIn <- sub("20160830","20161006",gsub("snv_mnv","indel", vcfFileIn))
vcfIndelFileOut <-  sub("20160830","20161006",gsub("snv_mnv","indel", vcfFileOut))
#' Load vcf
vcfIndel <- readVcf(vcfIndelFileIn, genome="GRCh37") #, param=ScanVcfParam(which=pos))
#' Add ID & gender
meta(header(vcfIndel)) <- rbind(meta(header(vcfIndel)), DataFrame(Value=c(ID, as.character(allGender[ID, "pred_gender"])), row.names=c("ID", "gender")))
#' Add driver genes
vcfIndel <- addFinalDriver(vcfIndel, driVers)
#' Add mutation copy numbers
i = header(vcfIndel)@header$INFO
exptData(vcfIndel)$header@header$INFO <- rbind(i,mcnHeader())
L <- computeMutCn(vcfIndel, bb, clusters=clusters, purity=purity, xmin=3, gender=as.character(allGender[ID, "pred_gender"]))
info(vcfIndel) <- cbind(info(vcfIndel), L$D)
#' Classify mutation
clsIndel <- classifyMutations(vcfIndel, reclassify='all')
info(vcfIndel)$CLS <- clsIndel
info(header(vcfIndel)) <- rbind(info(header(vcfIndel)), DataFrame(Number="1",Type="String",Description="Mutation classification: {clonal [early/late/NA], subclonal}", row.names="CLS"))
#' Save
writeVcf(vcfIndel, file=vcfIndelFileOut)
bgzip(vcfIndelFileOut, overwrite=TRUE)
unlink(vcfIndelFileOut)
save(bb, file=sub(".vcf$",".bb_granges.RData",vcfIndelFileOut))
save(vcfIndel, file=paste0(vcfIndelFileOut,".RData"))


#' ## PLOT

chrOffset <- cumsum(c(0,as.numeric(width(refLengths))))
names(chrOffset) <- c(seqlevels(refLengths), "NA")

for(v in c('vcf','vcfIndel')){
	vv <- get(v)
	pdf(file=sub(".vcf$",".pdf",get(paste0(v,"FileOut"))), 16,8)
	par(mar=c(3,3,1,1), bty="L", mgp=c(2,.5,0))
	col <- RColorBrewer::brewer.pal(9, "Set1")[c(3,4,2,1,9)]
	cls <- factor(paste(as.character(info(vv)$CLS)), levels = c(levels(info(vv)$CLS), "NA"))
	plot(start(vv) + chrOffset[as.character(seqnames(vv))], getAltCount(vv)/getTumorDepth(vv),col=col[cls], xlab='Position', ylab="VAF", pch=ifelse(info(vv)$pMutCNTail < 0.025 | info(vv)$pMutCNTail > 0.975, 4 , 16), ylim=c(0,1), xlim=c(0,chrOffset["MT"]))
	abline(v = chrOffset[1:24], lty=3)
	for(i in seq_along(bb)) try({
					s <- start(bb)[i]
					e <- end(bb)[i]
					x <- chrOffset[as.character(seqnames(bb)[i])]
					y <- bb$timing_param[[i]][,"f"]
					l <- bb$timing_param[[i]][,"pi.s"] * bb$timing_param[[i]][,"P.m.sX"]
					segments(s+x,y,e+x,y, lwd=l*4+.1)
					text(x=(s+e)/2 +x, y=y, paste(signif(bb$timing_param[[i]][,"m"],2),signif(bb$timing_param[[i]][,"cfi"]/purityPloidy[meta(header(vv))["ID",1],"purity"],2), sep=":"), pos=3, cex=0.5)
				})
	legend("topleft", pch=19, col=col, legend=paste(as.numeric(table(cls)), levels(cls)), bg='white')
	dev.off()
}

#plot(start(vcf) + w[as.character(seqnames(vcf))], qnorm(info(vcf)$pMutCNTail), col=col[cls], xlab='Position', ylab="pTail", pch=16)


