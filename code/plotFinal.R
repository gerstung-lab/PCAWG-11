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


s <- strsplit(vcfFileIn,"/")[[1]]
ID <- sub("\\..+", "", s[length(s)])


library(VariantAnnotation)
library(Matrix)
library(CoxHD)

#' ## CLUSTERS
# Load clusters
clusters <- try(loadClusters(ID))
NO_CLUSTER <- FALSE
if(class(clusters)=='try-error'){
	clusters <- data.frame(cluster=1, n_ssms=1000, proportion=1)
	NO_CLUSTER <- TRUE
}

#' ## VCF 
#' Load vcf
load(file=paste0(vcfFileOut,".RData"))

#' ## INDEL
vcfIndelFileIn <- sub("20160830","20161006",gsub("snv_mnv","indel", vcfFileIn))
vcfIndelFileOut <-  sub("20160830","20161006",gsub("snv_mnv","indel", vcfFileOut))
load(file=paste0(vcfIndelFileOut,".RData"))

#' Load BB
load(file=sub("snv_mnv","cn",sub(".vcf$",".bb_granges.RData",vcfFileOut)))
IS_WGD <- classWgd(bb)
colnames(mcols(bb)) <- sub("star.1","time.star",colnames(mcols(bb)) ) # Fix naming problem
save(bb, file=sub("snv_mnv","cn",sub(".vcf$",".bb_granges.RData",vcfFileOut)))


#' Load clusters & purity
load(file=sub("snv_mnv","clusters",sub(".vcf$",".clusters_purity.RData",vcfFileOut)))

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
plotTiming(bb)#, mcols(bb)[,c("type","time","time.lo","time.up")])
dev.off()
#plot(start(vcf) + w[as.character(seqnames(vcf))], qnorm(info(vcf)$pMutCNTail), col=col[cls], xlab='Position', ylab="pTail", pch=16)


