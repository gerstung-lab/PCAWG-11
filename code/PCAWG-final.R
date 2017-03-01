#' ---
#' title: PCAWG-11 Timing analyses
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 4
#'     number_sections: true
#'     auto_identifiers: true
#'     table_captions: true
#' author: Moritz Gerstung
#' ---

#' # PCAWG-11 Timing analyses

#+ Preliminaries, echo=FALSE
options(width=120)
pdf.options(pointsize=8)
knit_hooks$set(smallMar = function(before, options, envir) {
			if (before) par(mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0)) 
		})
opts_chunk$set(dev=c('my_png','pdf'), fig.ext=c('png','pdf'), fig.width=3, fig.height=3, smallMar=TRUE)
my_png <-  function(file, width, height, pointsize=12, ...) {
	png(file, width = 1.5*width, height = 1.5*height, units="in", res=72*1.5, pointsize=pointsize, ...)
}


#' ## Prelim
#' ### Libraries
library(VariantAnnotation)
setwd("/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/code")
source("functions.R")

#' ## Load data
#' ### SNV and MNV
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_008/snv_mnv"
#finalVcfSnv <- mclapply(dir(p, pattern="*.vcf.RData", full.names=TRUE), function(f){
#			e <- new.env()
#			load(f, envir=e)
#			return(e$vcf)
#		}, mc.cores=6)
#names(finalVcfSnv) <- sub(".conse.+","",dir(p, pattern="*.vcf.RData", full.names=FALSE))
#which(sapply(finalVcfSnv, class)=="try-error")

finalSnv <- list()
j <- 1
for(f in dir(p, pattern="*.vcf.RData", full.names=TRUE)){
	if(j %% 10 ==0) print(j); j <- j+1
	load(f)
	finalSnv[[f]] <- vcf
}
names(finalSnv) <- sub(".conse.+","",dir(p, pattern="*.vcf.RData", full.names=FALSE))

#' ### Copy number profiles
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_008/cn"
finalBB <- list()
for( f in dir(p, pattern="*.bb_granges.RData", full.names=TRUE)){
	load(f)
	finalBB[[f]] <- bb
}
names(finalBB) <- sub(".conse.+","",dir(p, pattern="*.bb_granges.RData", full.names=FALSE))

#' ### Indel
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_008/indel"
finalIndel <- list()
for( f in dir(p, pattern="*.vcf.RData", full.names=TRUE)){
	load(f)
	finalIndel[[f]] <- vcfIndel
}
names(finalIndel) <- sub(".conse.+","",dir(p, pattern="*.vcf.RData", full.names=FALSE))

#' ### Clusters and purity
finalClusters <- list()
finalPurity <- numeric()
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_008/clusters"
for( f in dir(p, pattern="*.RData", full.names=TRUE)){
	load(f)
	finalClusters[[f]] <- clusters
	finalPurity[f] <- purity
}
names(finalClusters) <- names(finalPurity) <- sub(".conse.+","",dir(p, pattern="*.RData", full.names=FALSE))


#' ## Distribution of Mutations
#' ### MAP genotypes
finalGenotypesSnv <- simplify2array(mclapply(finalSnv, getGenotype, mc.cores=6, useNA="always"))
finalGenotypesIndel <- simplify2array(mclapply(finalIndel, getGenotype, mc.cores=6, useNA="always"))
finalGenotypes <- aperm(abind::abind(subs=finalGenotypesSnv,indels=finalGenotypesIndel, along=5), c(1,5,2,3,4))
rm(finalGenotypesSnv,finalGenotypesIndel)

#' ### Probabilistic genotypes
finalGenotypesSnvP <- simplify2array(mclapply(finalSnv, probGenotype, mc.cores=2))
finalGenotypesIndelP <- simplify2array(mclapply(finalIndel, probGenotype, mc.cores=2))
finalGenotypesP <- aperm(abind::abind(subs=finalGenotypesSnvP,indels=finalGenotypesIndelP, along=4), c(1,4,2,3))
rm(finalGenotypesSnvP,finalGenotypesIndelP)

#save.image(file=paste0(Sys.Date(),".RData"))
#save(finalGenotypes, file=paste0(Sys.Date(),"-FinalGenotypes.RData"))

#' ### Overall distribution
f <- function(x) unlist(sapply(seq_along(x), function(i) rep(i, x[i])))
d <- t(asum(finalGenotypesP[,"subs",,], 1))
o <- order(droplevels(donor2type[sample2donor[rownames(d)]]), -d[,1]/rowSums(d))
I <- t(apply(d/rowSums(d), 1, function(x) f(mg14:::roundProp(x * 100,p=100))))
d <- t(asum(finalGenotypesP[,"indels",,], 1))
J <- t(apply(d/rowSums(d), 1, function(x) if(!any(is.nan(x))) f(mg14:::roundProp(x * 100,p=100)) else rep(NA,100)))
s <- cumsum(table(droplevels(donor2type[sample2donor[rownames(d)]][o])))

col <- RColorBrewer::brewer.pal(9, "Set1")[c(3,4,2,1,9)] ## Colors for early-subclonal

par(fig=c(0,1,0,1),mar=c(1,4,1,1)+.1, mgp=c(2,.5,0), mfrow=c(3,1), bty="n", las=2, xpd=FALSE)
image(z=I[o,], x=1:nrow(I), useRaster=TRUE, col=col, xaxt="n", ylab="SNVs")
abline(v=s, col='white')
image(z=J[o,], x=1:nrow(I),useRaster=TRUE, col=col, xaxt="n", ylab="Indels")
abline(v=s, col='white')
par(bty="n", xpd=NA)
plot(NA,NA, xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0, nrow(J)), ylim=c(0,1), xaxs="i", yaxs='i')
d0 <- s - diff(c(0,s))/2
d1 <- mg14:::mindist(d0,30)
segments(d0, 1.12, d0, 1.2)
segments(d0, 1.12, d1, 1.08)
segments(d1, 1.08, d1, 1)
mg14::rotatedLabel(x0 = d1, names(s), y0=1)
legend("bottom", fill=col, legend=paste(dimnames(finalGenotypes)[[3]]), bty="n", horiz=TRUE, title="Mutation timing")

t <- 12/8
dev.copy2pdf(file="finalMutationProp.pdf", width=9*t, height=2.7*t, pointsize=8*t)


#' #### Barplot drivers
p <- asum(finalGenotypesP, c(2,4))
g <- asum(finalGenotypes, c(2,4:5))
g <- g[order(rowSums(g), decreasing=TRUE),]
colnames(g) <- paste(colnames(g))
rownames(g) <- paste(rownames(g)) 
rownames(p) <- paste(rownames(p))
p <- p[rownames(g),]
w <- rowSums(g) > 0
w[1] <- FALSE

par(fig=c(0,1,0,1),mar=c(1,4,1,1)+.1, mgp=c(3,.5,0))
barplot(t(g[w,]), col=col, las=2, legend=TRUE, args.legend=list("topright", bty='n'), ylab="Number of cases", names.arg=rep("",sum(w)), border=NA)
#mg14::rotatedLabel(x=.Last.value,labels=rownames(g)[62:201], cex=.25)

u <- par("usr")
v <- c(
		grconvertX(u[1:2], "user", "ndc"),
		grconvertY(u[3:4], "user", "ndc")
)
v <- c( (v[1]+v[2])/3.33, v[2], (v[3]+v[4])/3, v[4] )
par( fig=v, new=TRUE, mar=c(0,0,0,0) )
barplot(t(g[2:51,]), col=col, las=2,  names.arg=rep("",50))
mg14::rotatedLabel(x=.Last.value,labels=rownames(g)[2:51], cex=.5)

dev.copy2pdf(file="finalDrivers.pdf", width=9, height=5, pointsize=8)

#' #### Barpot drivers - proportions
par(fig=c(0,1,0,1),mar=c(3,4,1,1)+.1, mgp=c(3,.5,0))
w <- rowSums(p) > 0
n <- 50
barplot(t(p /rowSums(p))[,w], width=c(rep(2,n+1), rep(0.2,sum(w)-n-1)), space=c(0,2,rep(0.1,n), rep(0,sum(w)-n-2)), col=col, border=NA, ylab="Fraction of mutations", names.arg=rep("",sum(w)))
mg14::rotatedLabel(x=.Last.value[1:(n+1)],labels=c("Genome-wide", rownames(p)[2:(n+1)]), cex=.5)

s <- 12/8
dev.copy2pdf(file="finalDriversProp.pdf", width=9*s, height=3*s, pointsize=8*s)

#' #### Cumulative effects
#+ genesCumulative
tt <- abind::abind(p[-1,], p[-1,] + 0.5, along=3)

par(mar=c(3,3,1,1), mgp=c(2,.5,0), bty="L")
r <- array(apply(tt/rep(colSums(tt), each=nrow(tt)), 3, function(x) apply(x, 2, function(y) cumsum(sort(y, decreasing=TRUE)))), dim=dim(tt))
plot(cumsum(sort(r[,1,1], decreasing=TRUE)), col=NA, type='s', xlab="Number of different driver genes", ylab="Fraction of driver mutations", log='', ylim=c(0,1), xlim=c(0,664), xaxs='i', yaxs='i')
for(j in 1:4) {
	polygon(c(1:nrow(r), nrow(r):1), c(r[,j,1], rev(r[,j,2])), col=paste0(col[j],"33"), border=NA)
	lines((r[,j,1]+r[,j,2])/2, col=col[j], lwd=2)
	points(y=0,x=which.min((r[,j,1]+r[,j,2])/2 < 0.5), col=col[j], pch="|")
}
legend("bottomright", col=col[1:4], lty=1, paste0(c("clonal [early]", "clonal [late]", "clonal [other]", "subclonal"), ", n=", round(colSums(p[-1,]))[-5]), bty="n")
abline(h=0.5, lty=3)

dev.copy2pdf(file="finalGenesCumulative.pdf", width=4,height=4)

par(mar=c(4,3,2.5,1), mgp=c(2,.5,0), bty="L")
d50 <- apply((r[,,1]+r[,,2])/2 < 0.5, 2, which.min)[c(1,3,2,4)]
b <- barplot(d50,las=2, col=col[c(1,3,2,4)], border=NA, ylab="Genes contributing 50% of driver mutations")
segments(b,apply(r[,,1] < 0.5, 2, which.min)[c(1,3,2,4)],b,apply(r[,,2] < 0.5, 2, which.min)[c(1,3,2,4)])
mg14::rotatedLabel(x=b,labels=c("clonal [early]", "clonal [late]", "clonal [other]", "subclonal")[c(1,3,2,4)])

dev.copy2pdf(file="finalGenes50.pdf", width=3,height=4)



#' ## WGD analyses

#' fix BB
for(ID in names(finalBB))
	finalBB[[ID]]$total_cn <- finalBB[[ID]]$major_cn+ finalBB[[ID]]$minor_cn

finalPloidy <- sapply(finalBB, averagePloidy)
names(finalPloidy) <- names(finalBB)

finalHom <- sapply(finalBB, averageHom)
names(finalHom) <- names(finalBB)

isWgd <- .classWgd(finalPloidy, finalHom)

plot(finalHom, finalPloidy, col=.classWgd( finalPloidy, finalHom)+1, xlim=c(0,1))

#' Adjust CI's
pseudo.count <- 5
for(ID in names(finalBB)){
	finalBB[[ID]]$time.up <- (pseudo.count + finalBB[[ID]]$n.snv_mnv * finalBB[[ID]]$time.up)/(pseudo.count + finalBB[[ID]]$n.snv_mnv)
	finalBB[[ID]]$time.lo <- (0 + finalBB[[ID]]$n.snv_mnv * finalBB[[ID]]$time.lo)/(pseudo.count + finalBB[[ID]]$n.snv_mnv)
}

fracGenomeWgdComp <- t(sapply(finalBB, function(bb) {
					fgw <- try(fractionGenomeWgdCompatible(bb)); 
					if(class(fgw)!='try-error') fgw
					else rep(NA,6)}))
rownames(fracGenomeWgdComp) <- names(finalBB)

wgdStar <- factor(rep(1,nrow(fracGenomeWgdComp)), levels=0:3, labels=c("unlikely","uninformative","likely","very likely"))
wgdStar[fracGenomeWgdComp[,"avg.ci"]<=0.75 & fracGenomeWgdComp[,"nt.total"]/chrOffset["MT"] >= 0.33 ] <- "likely"
wgdStar[fracGenomeWgdComp[,"nt.wgd"]/fracGenomeWgdComp[,"nt.total"] < 0.66] <- "unlikely"
wgdStar[wgdStar=="likely" & fracGenomeWgdComp[,"nt.wgd"]/fracGenomeWgdComp[,"nt.total"] > 0.8 & fracGenomeWgdComp[,"sd.wgd"] < 0.1 &  fracGenomeWgdComp[,"nt.total"]/chrOffset["MT"] > 0.5] <- "very likely"
names(wgdStar) <-  names(finalBB)
prop.table(table(wgdStar[!isWgd]))

wgdPoss <- !isWgd & 2.5 - 1.5 * finalHom <= finalPloidy

wgdStat <- factor(wgdPoss + 2*isWgd - wgdPoss*isWgd, labels=c("absent","possible","present"))
table(wgdStat, wgdStar)

pdf("WGD-timing.pdf", 12,6.5)
j <- 1
for(ID in names(finalBB)[isWgd]){
	if(j%%100 == 0) print(j); j <- j+1
#	for(f in dir("../final/annotated_007/snv_mnv", pattern=paste0(ID,".+bb_granges.RData"), full.names=TRUE)) load(f)
#	t <- bbToTime(bb)
	par(mfrow=c(2,1), mar=c(3,3,2,1), mgp=c(2,.5,0), bty="L", cex=1, las=2)
	plotBB(finalBB[[ID]], ylim=c(0,8))
	title(main=paste0(ID,", ", donor2type[sample2donor[ID]], ", ploidy=",round(finalPloidy[ID],2), ", hom=",round(finalHom[ID],2)), font.main=1, line=0)
	par(mar=c(3,3,2,1))
	plotTiming(finalBB[[ID]])
	abline(h=fracGenomeWgdComp[ID,"time.wgd"], lty=3)
	title(main=paste0("Timeable=", round(fracGenomeWgdComp[ID,2]/chrOffset["MT"]*100), "%, WGD=",round(fracGenomeWgdComp[ID,1]/fracGenomeWgdComp[ID,2]*100), "%, sd.WGD=",round(fracGenomeWgdComp[ID,'sd.wgd'],2), "%, avg.CI=",round(fracGenomeWgdComp[ID,'avg.ci'],2),", verdict=",wgdStar[ID]), font.main=1, line=0)
}
dev.off()

#' BB without VCF
otherBB <- sapply(setdiff(sub("\\..+","",dir(bbPath)), names(finalBB)), function(ID) loadConsensusCNA(ID, purity=purityPloidy[ID, 'purity']))
for(ID in names(otherBB))
	otherBB[[ID]]$total_cn <- otherBB[[ID]]$major_cn+ otherBB[[ID]]$minor_cn

otherPloidy <- sapply(otherBB, averagePloidy)
otherHom <- sapply(otherBB, averageHom)
otherWgd <- .classWgd(otherPloidy, otherHom)
otherPoss <- !otherWgd & 2.5 - 1.5 * otherHom <= otherPloidy
otherStat <- factor(otherPoss + 2*otherWgd - otherPoss*otherWgd, levels=0:2,labels=c("absent","possible","present"))


tab <- data.frame(WGD_call = wgdStat, WGD_timing=wgdStar, ploidy=finalPloidy, hom=finalHom, fracGenomeWgdComp)
tab <- rbind(tab, data.frame(WGD_call=otherStat, WGD_timing=NA, ploidy=otherPloidy, hom=otherHom, nt.wgd=NA, nt.total=NA, time.wgd=NA, sd.wgd=NA,avg.ci=NA, sd.all=NA))
write.table(file="WGD-info.txt", tab, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")