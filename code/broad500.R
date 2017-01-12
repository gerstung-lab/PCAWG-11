# TODO: Add comment
# 
# Author: mg14
###############################################################################


args <- commandArgs(trailingOnly = TRUE)

source("functions.R")

vcfFileIn <- args[1]
vcfFileOut <- args[2]
xmin <- as.numeric(args[3])
rho <- as.numeric(args[4])

print(vcfFileIn)
print(vcfFileOut)

library(VariantAnnotation)
library(Matrix)
library(CoxHD)
library(igraph)

r = "/lustre/scratch112/sanger/cgppipe/PanCancerReference/genome.fa.gz" #meta(header(v))["reference",]
refLengths <- scanFaIndex(file=r)

dpPath <- '/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/broad500/Subclonal_Structure'
dpFiles <- dir(dpPath, pattern="subclonal_structure.txt", recursive=TRUE)

s <- strsplit(vcfFileIn,"/")[[1]]
ID <- sub("\\..+", "", s[length(s)])


print(ID)
clusters <- loadClusters(ID)

purityPloidy <- read.table("/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/broad500/pp_table.txt", header=TRUE, sep='\t')
rownames(purityPloidy) <- purityPloidy$sample
purityPloidy <- purityPloidy[,2:3]

if(all(is.na(purityPloidy[ID,]))) # Missing purity
	purityPloidy[ID,] <- c(max(clusters$proportion),NA)

purity <- purityPloidy[ID,"purity"]

clusters$proportion <- clusters$ccf * purity


# Load BB
bbPath <- '/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/broad500/Segments'
bb <- loadBB(ID)
if(length(bb)==0){ # Missing BB, use consensus CN
	cnPath <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/cn/consensus_001/all"
	file <- paste0(cnPath, "/",ID,"_segments.txt")
	tab <- read.table(file, header=TRUE, sep='\t')
	bb <- GRanges(tab$chromosome, IRanges(tab$start, tab$end), strand="*", tab[-3:-1])
	#m <- read.table(gzfile(paste0(basePath, "/0_multiplicity/",ID,"_multiplicity.txt.gz")), header=TRUE)
	#bb <- GRanges(m$chr, IRanges(m$pos, width=1), copy_number=m$tumour_copynumber, major_cn=m$nMaj1, minor_cn=m$nMin1, clonal_frequency=purityPloidy[ID,'purity'])
	#meta(header(vcf)) <- rbind(meta(header(vcf)), DataFrame(Value="False", row.names="Battenberg"))
}
bb$clonal_frequency <- bb$ccf * purity	

# Load vcf
vcf <- readVcf(vcfFileIn, genome="GRCh37") #, param=ScanVcfParam(which=pos))

# Add ID & gender
meta(header(vcf)) <- rbind(meta(header(vcf)), DataFrame(Value=c(ID, as.character(gender[ID, "pred_gender"])), row.names=c("ID", "gender")))


#' Add mutation copy numbers
# vcf <-  addMutCn(vcf, bb, clusters)
i = header(vcf)@header$INFO
exptData(vcf)$header@header$INFO <- rbind(i,mcnHeader())
L <- computeMutCn(vcf, bb, clusters, xmin=xmin, rho=rho)
info(vcf) <- cbind(info(vcf), L$D)
bb$timing_param <- L$P 


#' Classify mutations
cls <- classifyMutations(vcf, reclassify='all')

info(vcf)$CLS <- cls
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number="1",Type="String",Description="Mutation classification: {clonal [early/late/NA], subclonal}", row.names="CLS"))

#' Save output
#fnew <- sub(".vcf",".complete_annotation.vcf",vcfFileOut)
writeVcf(vcf, file=vcfFileOut)
save(bb, file=sub(".vcf$",".bb_granges.RData",vcfFileOut))
bgzip(vcfFileOut, overwrite=TRUE)
save(vcf, file=paste0(vcfFileOut,".RData"))

w <- cumsum(c(0,as.numeric(width(refLengths))))
names(w) <- c(seqlevels(refLengths), "NA")

#' Some quantitative assessment
truth <- read.table(paste0("../broad500/Mut_Assign/",ID,".mutation_assignments.txt"), header=TRUE, sep="\t")
truth <- GRanges(truth$chr, IRanges(truth$pos, width=1), cluster=truth$cluster)
seqlevels(truth) <-paste( 1:22)
truth <- sort(truth)

pdf(file=sub(".vcf$",".pdf",vcfFileOut), 16,8)
par(mar=c(3,3,3,1), bty="L", mgp=c(2,.5,0))
col <- RColorBrewer::brewer.pal(4, "Set1")[c(3,4,2,1)]
plot(start(vcf) + w[as.character(seqnames(vcf))], getAltCount(vcf)/getTumorDepth(vcf),col=col[cls], xlab='Position', ylab="VAF", pch=ifelse(info(vcf)$pMutCNTail < 0.025 | info(vcf)$pMutCNTail > 0.975, 4 , 16))
title(paste0(ID,", " ,round(100*mean(info(vcf)$CNF == clusters$proportion[truth$cluster+1]),1), "% correct, ", round(100*mean((info(vcf)$CNF==purity) == (clusters$proportion[truth$cluster+1]==purity)),1), "% clonal v subclonal"),font=1, line=1)
abline(v = w, lty=3)
for(i in seq_along(bb)) try({
	s <- start(bb)[i]
	e <- end(bb)[i]
	x <- w[as.character(seqnames(bb)[i])]
	y <- bb$timing_param[[i]][,"f"]
	l <- bb$timing_param[[i]][,"pi.s"] * bb$timing_param[[i]][,"P.m.sX"]
	segments(s+x,y,e+x,y, lwd=l*4+.1)
	text(x=(s+e)/2 +x, y=y, paste(signif(bb$timing_param[[i]][,"m"],2),signif(bb$timing_param[[i]][,"cfi"]/purityPloidy[meta(header(vcf))["ID",1],"purity"],2), sep=":"), pos=3, cex=0.5)
})
legend("topleft", pch=19, col=col, legend=levels(cls))
#dev.copy(device=png, file=sub(".vcf$",".png",vcfFileOut), width = 16, height = 8, units = "in", pointsize = 12, res=72)
#dev.off()
dev.off()

#plot(start(vcf) + w[as.character(seqnames(vcf))], qnorm(info(vcf)$pMutCNTail), col=col[cls], xlab='Position', ylab="pTail", pch=16)


