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
clusters = wccClusters[[ID]]
purity <- purityPloidy[ID,'purity']

#' ## COPYNUMBER
bb <- loadConsensusCNA(ID, purity=purityPloidy[ID, 'purity'])
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

#' Drivers
info(vcf)$DG <- matchDrivers(vcf, finalDrivers)

#' Save output
writeVcf(vcf, file=vcfFileOut)
bgzip(vcfFileOut, overwrite=TRUE)
unlink(vcfFileOut)
save(vcf, file=paste0(vcfFileOut,".RData"))



#' Save BB
bb$n.indel <- countOverlaps(bb, vcfIndel)
seqlevels(bb) <- c(1:22, "X","Y")
bb <- sort(bb)
save(bb, file=sub("snv_mnv","cn",sub(".vcf$",".bb_granges.RData",vcfFileOut)))

#' Save clusters & purity
save(clusters, purity, file=sub("snv_mnv","clusters",sub(".vcf$",".clusters_purity.RData",vcfFileOut)))

#' ## PLOT
pdf(file=sub("snv_mnv","other",sub(".vcf$",".pdf",vcfFileOut)), 6,6)
plotSample(ID, vcf=vcf, bb=bb)
dev.off()
#plot(start(vcf) + w[as.character(seqnames(vcf))], qnorm(info(vcf)$pMutCNTail), col=col[cls], xlab='Position', ylab="pTail", pch=16)


