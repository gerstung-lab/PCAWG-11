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

#+ Preliminaries, echo=FALSE
options(width=120)
pdf.options(pointsize=8)
knit_hooks$set(smallMar = function(before, options, envir) {
			if (before) par(mar=c(3,3,1,1), bty="L", mgp=c(2,.5,0), tcl=-0.25, las=1) 
		})
opts_chunk$set(dev=c('my_png','pdf'), fig.ext=c('png','pdf'), fig.width=3, fig.height=3, smallMar=TRUE)
my_png <-  function(file, width, height, pointsize=12, ...) {
	png(file, width = 1.5*width, height = 1.5*height, units="in", res=72*1.5, pointsize=pointsize, ...)
}


#' # Prelim
#' ## Libraries
library(VariantAnnotation)
setwd("~/pcawg/code")
source("PCAWG-functions.R")
MC_CORES=2

#+ evalOff, echo=FALSE
opts_chunk$set(eval=FALSE)
load("2018-02-22-PCAWG-final.RData")
source("PCAWG-functions.R")


#' # Load data
#' ## Whitelist
#' ### SNV and MNV
#+ loadSNV
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_012/snv_mnv"
finalSnv <- list()
j <- 1
for(f in dir(p, pattern="*.vcf.RData", full.names=TRUE)){
	if(j %% 10 ==0) print(j); j <- j+1
	load(f)
	finalSnv[[f]] <- vcf
}
names(finalSnv) <- sub(".conse.+","",dir(p, pattern="*.vcf.RData", full.names=FALSE))

#' ### Copy number profiles
#+ loadBB
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_012/cn"
finalBB <- list()
for( f in dir(p, pattern="*.bb_granges.RData", full.names=TRUE)){
	load(f)
	colnames(mcols(bb)) <- sub("star.1","time.star",colnames(mcols(bb)) ) # Fix naming problem
	finalBB[[f]] <- bb
}
names(finalBB) <- sub(".conse.+","",dir(p, pattern="*.bb_granges.RData", full.names=FALSE))

#' ### Indel
#+ loadIndel
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_012/indel"
finalIndel <- list()
for( f in dir(p, pattern="*.vcf.RData", full.names=TRUE)){
	load(f)
	finalIndel[[f]] <- vcfIndel
}
names(finalIndel) <- sub(".conse.+","",dir(p, pattern="*.vcf.RData", full.names=FALSE))

#' ### Clusters and purity
#+ loadClusters
finalClusters <- list()
finalPurity <- numeric()
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_012/clusters"
for( f in dir(p, pattern="*.RData", full.names=TRUE)){
	load(f)
	finalClusters[[f]] <- clusters
	finalPurity[f] <- purity
}
names(finalClusters) <- names(finalPurity) <- sub(".conse.+","",dir(p, pattern="*.RData", full.names=FALSE))


#' ## Update drivers
#' ### Update finalDrivers
#+ updateDrivers
finalDriversAnnotated <- finalDrivers
d <- info(finalSnv[[3]])[seq_along(finalDriversAnnotated),19:32]
#d[,] <- NA
mcols(finalDriversAnnotated)[colnames(d)] <- d
for(i in seq_along(finalDriversAnnotated)){
	if(finalDriversAnnotated$mut_type[i] %in% c("snv","mnv")){
		v <- finalSnv[[as.character(finalDriversAnnotated$sample[i])]]
	}else{
		v <- finalIndel[[as.character(finalDriversAnnotated$sample[i])]]
	}
	j <- findOverlaps(finalDriversAnnotated[i], v, select='first')
	if(!is.na(j))
		mcols(finalDriversAnnotated)[i,colnames(d)] <- info(v)[j, colnames(d)]
	else
		mcols(finalDriversAnnotated)[i,colnames(d)] <- NA
}


#' ## Graylisted data
#' ### SNV and MNV
#+ loadSnvGray
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_012/graylist/snv_mnv"
finalSnvGray <- list()
j <- 1
for(f in dir(p, pattern="*.vcf.RData", full.names=TRUE)){
	if(j %% 10 ==0) print(j); j <- j+1
	load(f)
	finalSnvGray[[f]] <- vcf
}
names(finalSnvGray) <- sub(".conse.+","",dir(p, pattern="*.vcf.RData", full.names=FALSE))
finalSnv[names(finalSnvGray)] <- finalSnvGray

#' ### Copy number profiles
#+ loadBBGray
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_012/graylist/cn"
finalBBGray <- list()
for( f in dir(p, pattern="*.bb_granges.RData", full.names=TRUE)){
	load(f)
	colnames(mcols(bb)) <- sub("star.1","time.star",colnames(mcols(bb)) ) # Fix naming problem
	finalBBGray[[f]] <- bb
}
names(finalBBGray) <- sub(".conse.+","",dir(p, pattern="*.bb_granges.RData", full.names=FALSE))
finalBB[names(finalBBGray)] <- finalBBGray


#' ### Indel
#+ loadIndelGray
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_012/graylist/indel"
finalIndelGray <- list()
for( f in dir(p, pattern="*.vcf.RData", full.names=TRUE)){
	load(f)
	finalIndelGray[[f]] <- vcfIndel
}
names(finalIndelGray) <- sub(".conse.+","",dir(p, pattern="*.vcf.RData", full.names=FALSE))
finalIndel[names(finalIndelGray)] <- finalIndelGray


#' ### Clusters and purity
#+ loadClustersGray
finalClustersGray <- list()
finalPurityGray <- numeric()
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_012/graylist/clusters"
for( f in dir(p, pattern="*.RData", full.names=TRUE)){
	load(f)
	finalClustersGray[[f]] <- clusters
	finalPurityGray[f] <- purity
}
names(finalClustersGray) <- names(finalPurityGray) <- sub(".conse.+","",dir(p, pattern="*.RData", full.names=FALSE))
finalClusters[names(finalClustersGray)] <- finalClustersGray
finalPurity <- c(finalPurity, finalPurityGray)

whiteList <- seq_along(finalSnv) %in% 1:2703
grayList <- !whiteList

#' # QC
#+ QC
q1 <- sapply(finalSnv, function(vcf) mean(abs(0.5- info(vcf)$pMutCNTail) > 0.495 , na.rm=TRUE))
q5 <- sapply(finalSnv, function(vcf) mean(abs(0.5- info(vcf)$pMutCNTail) > 0.475 , na.rm=TRUE))

#+ QCboxplot, eval=TRUE
par(mfrow=c(1,1))
boxplot(1-q5 ~ donor2type[sample2donor[names(finalSnv)]], las=2, ylab="Fraction of data inside theoretical 95% CI")
abline(h=0.95, lty=3)

#+ QQplots, fig.width=12, fig.height=12, eval=TRUE
#pdf("QQplots.pdf", 4,4, pointsize=8)
par(mfrow=c(5,5))
for(i in seq_along(finalSnv)[1:25]){
	n <- nrow(finalSnv[[i]])
	qqnorm(qnorm(info(finalSnv[[i]])$pMutCNTail[sample(1:n, min(1e4,n))]), main=paste(substr(names(finalSnv)[i],1,8), "Q5 =", signif(q5[i],2), ", Q1 =", signif(q1[i],2)), xlim=c(-5,5), ylim=c(-5,5), pch=16)
	abline(0,1, col='red')
}
#dev.off()


#' # Driver genotypes
#' ## MAP genotypes
#+ finalGenotypes
finalGenotypesSnv <- simplify2array(mclapply(finalSnv[whiteList], getGenotype, mc.cores=MC_CORES, useNA="always"))
finalGenotypesIndel <- simplify2array(mclapply(finalIndel[whiteList], getGenotype, mc.cores=MC_CORES, useNA="always"))
finalGenotypes <- aperm(abind::abind(subs=finalGenotypesSnv,indels=finalGenotypesIndel, along=5), c(1,5,2,3,4))
rm(finalGenotypesSnv,finalGenotypesIndel)

#' ## Probabilistic genotypes
#+ finalGenotypesP
finalGenotypesSnvP <- simplify2array(mclapply(finalSnv[whiteList], probGenotype, mc.cores=MC_CORES))
finalGenotypesIndelP <- simplify2array(mclapply(finalIndel[whiteList], probGenotype, mc.cores=MC_CORES))
finalGenotypesP <- aperm(abind::abind(subs=finalGenotypesSnvP,indels=finalGenotypesIndelP, along=4), c(1,4,2,3))
rm(finalGenotypesSnvP,finalGenotypesIndelP)

#' ## Probabilistic genotypes - tail prob (QC)
#+ finalGenotypesQ
finalGenotypesSnvQ <- simplify2array(mclapply(finalSnv[whiteList], probGenotypeTail, mc.cores=MC_CORES))
finalGenotypesIndelQ <- simplify2array(mclapply(finalIndel[whiteList], probGenotypeTail, mc.cores=MC_CORES))
finalGenotypesQ <- aperm(abind::abind(subs=finalGenotypesSnvQ,indels=finalGenotypesIndelQ, along=3), c(1,3,2))
rm(finalGenotypesSnvQ,finalGenotypesIndelQ)

#' # Save output
#+ saveOut, eval=FALSE
save.image(file=paste0(Sys.Date(),"-PCAWG-final.RData"))
save(finalGenotypes, finalGenotypesP, finalGenotypesQ, file=paste0(Sys.Date(),"-FinalGenotypes.RData"))

#+ evalOn, eval=TRUE, echo=FALSE
opts_chunk$set(eval=TRUE)

#' # Distribution of Mutations
#' ## Duplicated samples
w <- names(finalSnv)
n <- names(which(table(sample2donor[w]) > 1)) # donors
s <- w[w %in% names(sample2donor[sample2donor %in% n])] # multisamples
d <- setdiff(sample2donor[w], sample2donor[levels(finalDrivers$sample)]) # donors w/o driver
u <- sample2donor[s[sample2donor[s] %in% intersect(d,n)]]
selectedSamples <- !w %in% setdiff(s[!s %in% finalDrivers$sample ], names(u)[!duplicated(u)])
uniqueSamples <- !duplicated(sample2donor[names(finalSnv)])

#' ## Overall distribution
#' ### Subs or indels
f <- function(x) unlist(sapply(seq_along(x), function(i) rep(i, x[i])))
d <- t(asum(finalGenotypesP[,"subs",,], 1))
o <- order(droplevels(donor2type[sample2donor[rownames(d)]]), -d[,1]/rowSums(d))
I <- t(apply(d/rowSums(d), 1, function(x) f(mg14:::roundProp(x * 100,p=100))))
d <- t(asum(finalGenotypesP[,"indels",,], 1))
J <- t(apply(d/rowSums(d), 1, function(x) if(!any(is.nan(x))) f(mg14:::roundProp(x * 100,p=100)) else rep(NA,100)))
s <- cumsum(table(droplevels(donor2type[sample2donor[rownames(d)]][o])))

#+ finalMutationsProb, fig.width=9, fig.height=2.7
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

#t <- 12/8
#dev.copy2pdf(file="finalMutationProp.pdf", width=9*t, height=2.7*t, pointsize=8*t)

#' ### Subs + indels
f <- function(x) unlist(sapply(seq_along(x), function(i) rep(i, x[i])))
d <- t(asum(finalGenotypesP[,,,], 1:2))
o <- order(droplevels(donor2type[sample2donor[rownames(d)]]), -d[,1]/rowSums(d))
I <- t(apply(d/rowSums(d), 1, function(x) f(mg14:::roundProp(x * 100,p=100))))
s <- cumsum(table(droplevels(donor2type[sample2donor[rownames(d)]][o])))

col <- RColorBrewer::brewer.pal(9, "Set1")[c(3,4,2,1,9)] ## Colors for early-subclonal

#+ finalMutationPropAll, fig.width=9, fig.height=1.8
par(fig=c(0,1,0,1),mar=c(1,4,1,1)+.1, mgp=c(2,.5,0), mfrow=c(2,1), bty="n", las=2, xpd=FALSE)
image(z=I[o,], x=1:nrow(I), useRaster=TRUE, col=col, xaxt="n", ylab="Point mutations")
abline(v=s, col='white')
par(bty="n", xpd=NA)
plot(NA,NA, xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0, nrow(I)), ylim=c(0,1), xaxs="i", yaxs='i')
d0 <- s - diff(c(0,s))/2
d1 <- mg14:::mindist(d0,30)
segments(d0, 1.12, d0, 1.2)
segments(d0, 1.12, d1, 1.08)
segments(d1, 1.08, d1, 1)
mg14::rotatedLabel(x0 = d1, names(s), y0=1)
legend("bottom", fill=col, legend=paste(dimnames(finalGenotypes)[[3]]), bty="n", horiz=TRUE, title="Mutation timing")

#t <- 12/8
#dev.copy2pdf(file="finalMutationPropAll.pdf", width=9*t, height=1.8*t, pointsize=8*t)


#' ### Barplot drivers
p <- asum(finalGenotypesP[,,,selectedSamples[whiteList]], c(2,4))
g <- asum(finalGenotypes[,,,,selectedSamples[whiteList]], c(2,4:5))
g <- g[order(rowSums(g), decreasing=TRUE),]
colnames(g) <- paste(colnames(g))
rownames(g) <- paste(rownames(g)) 
rownames(p) <- paste(rownames(p))
p <- p[rownames(g),]
w <- rowSums(g) > 0
w[1] <- FALSE

#' Unique loci
l <- sapply(strsplit(rownames(p),"::"), function(x) paste(x[2:3], collapse=":"))
pu <- t(sapply(unique(l), function(u) asum(p[l==u,,drop=FALSE], 1)))
gu <- t(sapply(unique(l), function(u) asum(g[l==u,,drop=FALSE], 1)))
gu <- gu[order(rowSums(gu), decreasing=TRUE),]
pu <- pu[rownames(gu),]
wu <- rowSums(gu) > 0
wu[1] <- FALSE

#+ finalDrivers, fig.width=9, fig.height=5
par(fig=c(0,1,0,1),mar=c(1,4,1,1)+.1, mgp=c(3,.5,0))
barplot(t(gu[wu,]), col=col, las=2, legend=TRUE, args.legend=list("topright", bty='n'), ylab="Number of cases", names.arg=rep("",sum(wu)), border=NA)
#mg14::rotatedLabel(x=.Last.value,labels=rownames(g)[62:201], cex=.25)

u <- par("usr")
v <- c(
		grconvertX(u[1:2], "user", "ndc"),
		grconvertY(u[3:4], "user", "ndc")
)
v <- c( (v[1]+v[2])/3.33, v[2], (v[3]+v[4])/3, v[4] )
par( fig=v, new=TRUE, mar=c(0,0,0,0) )
b <- barplot(t(gu[2:51,]), col=col, las=2,  names.arg=rep("",50))
mg14::rotatedLabel(x=b,labels=rownames(gu)[2:51], cex=.5)
#dev.copy2pdf(file="finalDrivers.pdf", width=9, height=5, pointsize=8)

#' ### Barpot drivers - proportions
#+ finalDriversProp, fig.width=9, fig.height=3
par(fig=c(0,1,0,1),mar=c(3,4,1,1)+.1, mgp=c(3,.5,0))
w <- rowSums(pu) > 0
n <- 50
b <- barplot(t(pu /rowSums(pu))[,w], width=c(rep(2,n+1), rep(0.2,sum(w)-n-1)), space=c(0,2,rep(0.1,n), rep(0,sum(w)-n-2)), col=col, border=NA, ylab="Fraction of mutations", names.arg=rep("",sum(w)))
mg14::rotatedLabel(x=b[1:(n+1)],labels=c("Genome-wide", rownames(pu)[2:(n+1)]), cex=.5)

#s <- 12/8
#dev.copy2pdf(file="finalDriversProp.pdf", width=9*s, height=3*s, pointsize=8*s)

#' ### Cumulative effects
#+ genesCumulative, fig.width=4, fig.height=4
tt <- abind::abind(pu[-1,], pu[-1,] + 0.5, along=3)

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
#dev.copy2pdf(file="finalGenesCumulative.pdf", width=4,height=4)

#+ finalGenes50, fig.width=3, fig.height=4
par(mar=c(4,3,2.5,1), mgp=c(2,.5,0), bty="L")
d50 <- apply((r[,,1]+r[,,2])/2 < 0.5, 2, which.min)[c(1,3,2,4)]
b <- barplot(d50,las=2, col=col[c(1,3,2,4)], border=NA, ylab="Genes contributing 50% of driver mutations")
segments(b,apply(r[,,1] < 0.5, 2, which.min)[c(1,3,2,4)],b,apply(r[,,2] < 0.5, 2, which.min)[c(1,3,2,4)])
mg14::rotatedLabel(x=b,labels=c("clonal [early]", "clonal [late]", "clonal [other]", "subclonal")[c(1,3,2,4)])
#dev.copy2pdf(file="finalGenes50.pdf", width=3,height=4)



#' # Whole-genome duplications
finalPloidy <- sapply(finalBB, averagePloidy)
names(finalPloidy) <- names(finalBB)

finalHom <- sapply(finalBB, averageHom)
names(finalHom) <- names(finalBB)

isWgd <- .classWgd(finalPloidy, finalHom)

#+ wdgHomPloidy, fig.widht=4, fig.height=4
plot(finalHom, finalPloidy, col=.classWgd( finalPloidy, finalHom)+1, xlim=c(0,1))

fracGenomeWgdComp <- t(sapply(finalBB, function(bb) {
					fgw <- try(fractionGenomeWgdCompatible(bb)); 
					if(class(fgw)!='try-error') fgw
					else rep(NA,10)}))
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



#' # Coamplification and WGD
d <- fracGenomeWgdComp
i <- d[,"avg.ci"]<=0.5 & d[,"chr.all"] > 2 #&  fracGenomeWgdComp[,"nt.total"]/chrOffset["MT"] >= 0.1
timingClass <- paste(ifelse(isWgd,"WGD","near-diploid"), ifelse(!i, "uninformative",""))
timingClass[i] <- paste0(timingClass[i], ifelse(d[i,"nt.wgd"]/d[i,"nt.total"] > 0.75,"sync","async"))
#timingClass[i] <- paste0(timingClass[i], cut(fracGenomeWgdComp[i,"nt.wgd"]/fracGenomeWgdComp[i,"nt.total"], c(0,0.5,0.8,1), include.lowest=TRUE))
timingClass <- factor(timingClass)

#+ timingClass, fig.width=4, fig.height=4
#pdf("TimingClass.pdf", 4,4)
c <- c(RColorBrewer::brewer.pal(9, "Pastel1"),"#DDDDDD")
t <- table(timingClass)[c(1:3,6,4:5)]
pie(t, init.angle=90, labels=paste0(names(t), ",\nn=", t), col=c[c(1,1,9,10,2,2)], density=c(36,NA,NA,NA,36,NA))
t <- table(isWgd)
par(new=TRUE)
pie(t, labels=c("",""), col=NA, lwd=5, lty=1, init.angle=90)
#dev.off()

colnames(d) <- c("ntCoamp","ntAmp","timeCoamp","segCoamp","segAmp","chrCoamp","chrAmp", "sdTimeCoamp","avgCiSeg","sdAllSeg")
tab <- data.frame(avgPloidy=finalPloidy, avgHom=finalHom, isWgd=isWgd, d, informative=i, timingClass=timingClass)
#tab <- rbind(tab, data.frame(WGD_call=otherStat, WGD_timing=NA, ploidy=otherPloidy, hom=otherHom, nt.wgd=NA, nt.total=NA, time.wgd=NA, sd.wgd=NA,avg.ci=NA, sd.all=NA))
write.table(file=paste0(Sys.Date(),"-Timing-info.txt"), tab, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")


#' ## Timing examples

#+ timingExamples, fig.width=4, fig.height=4
w <- which(wgdStar=="likely" & !isWgd)
#pdf(paste0(names(w[1]), ".pdf"), 4,4, pointsize=8)
plotSample(w[1])
plotSample(w[2])
plotSample(w[3])
#dev.off()

w <- which(wgdStar=="very likely" & isWgd)
#pdf(paste0(names(w[1]), ".pdf"), 4,4, pointsize=8)
plotSample(w[1])
plotSample(w[2])
plotSample(w[9])
#dev.off()

w <- which(wgdStar=="unlikely" & !isWgd & fracGenomeWgdComp[,"nt.total"]/chrOffset["MT"] > 0.25 & fracGenomeWgdComp[,"avg.ci"] < 0.5)
#pdf(paste0(names(w[1]), ".pdf"), 4,4, pointsize=8)
plotSample(w[1])
plotSample(w[2])
plotSample(w[3])
#dev.off()

#' ## GBM examples
#+ timingExamplesGbm, fig.width=4, fig.height=4
w <- which(fracGenomeWgdComp[,"time.wgd"]<0.1 & fracGenomeWgdComp[,"nt.total"]/chrOffset["MT"] > 0.1 &  !isWgd & donor2type[sample2donor[names(finalBB)]]=="CNS-GBM")
#pdf(paste0(names(w[1]), ".pdf"), 4,4, pointsize=8)
plotSample(w[1])
plotSample(w[2])
plotSample(w[3])
#dev.off()

#' ## Relationship with mutation rates
#' Calculate number of substitutions and deciles per tumour type
n <- nSub <- sapply(finalSnv, nrow)
n[tab$timeCoamp==0] <- NA
q <- unlist(sapply(split(n, donor2type[sample2donor[names(finalSnv)]]), function(x) as.numeric(cut(x, {if(sum(!is.na(x))>1) quantile(x, seq(0,1,0.1), na.rm=TRUE) else 1:10}, include.lowest=TRUE))))
m <- match(names(finalSnv),unlist(split(names(finalSnv), donor2type[sample2donor[names(finalSnv)]])))
t <- tab$timeCoamp
table(decSub=q[m], time=cut(t, seq(0,1,0.1)))

#' Also calculate deciles of timing per tumour type
t[t==0] <- NA
r <- unlist(sapply(split(t, donor2type[sample2donor[names(finalSnv)]]), function(x) as.numeric(cut(x, {if(sum(!is.na(x))>1 & length(unique(x)) > 2) quantile(jitter(x), seq(0,1,0.1), na.rm=TRUE) else 1:10}, include.lowest=TRUE))))
table(decSub=q[m], decTime=r[m])

#' Plot 
#+ timeNsub, fig.height=3, fig.width=4
#pdf("timeNsub.pdf", 3, 2.5, pointsize=8)
par(mar=c(3,4,1,1), bty="n", mgp=c(2,.5,0), las=1, tcl=-.25) 
d <- as.character(donor2type[sample2donor[names(finalSnv)]])
lineCol <- tissueColors
lineCol[grep("Lung", names(lineCol))] <- "black"
plot(t, nSub, log='y', bg=tissueColors[d], pch=21, xlab="Typical amplification time", ylab="", cex=.66, lwd=.5, yaxt="n", ylim=c(100,3e6))
mtext(side=2, "Number of SNVs", line=3, las=3)
u <- round(par("usr")[3:4])
a <- axisTicks(par("usr")[3:4], log=TRUE)
axis(side=2, at=a, labels=prettyNum(a))
b <- sapply(a[-length(a)], function(x) (1:10)*x)
axis(side=2, at=b, labels=rep("", length(b)), tcl=-.1)
#dev.off()

#' ## Overview histograms
#' Some functions
aggregatePerChromosome <- function(bb, isWgd=FALSE){
	.aggregateSegments <- function(m){
		#m <- mcols(bb)
		t <- weighted.mean(m$time, m$n.snv_mnv, na.rm=TRUE)
		n <- sum(m$n.snv_mnv, na.rm=TRUE)
		sd <- sd(m$time, na.rm=TRUE)
		ci <- weighted.mean(m$time.up-m$time.lo, m$n.snv_mnv, na.rm=TRUE)
		w <- sum(m$width[!is.na(m$time)], na.rm=TRUE)
		c(time=t, n=n, sd=sd, ci=ci,w=w)
	}
#	if(!isWgd){
		s <- split(as.data.frame(bb)[,c("time","time.up","time.lo","n.snv_mnv","width")], seqnames(bb))
		r <- t(sapply(s, .aggregateSegments))
		r <- r[c(1:22,"X"),]
#	}else{
		w <- .aggregateSegments(as.data.frame(bb))
		r <- rbind(r,WGD=w)
#	}
	return(r)
}

allChrAgg <- simplify2array(mclapply(finalBB, aggregatePerChromosome, mc.cores=2))


t <- allChrAgg[1:23,"time",!isWgd]
t[allChrAgg[1:23,"w",!isWgd] < diff(chrOffset)[1:23]*.33] <- NA

s <- split(as.data.frame(t(t)), droplevels(donor2type[sample2donor[names(finalSnv)]])[!isWgd])
n <- 10


at <- function(x, n){
	if(sum(!is.na(x))<3) return(rep(sum(!is.na(x))/n,n))
	bw=if(sum(!is.na(x))< 6) 0.5 else "nrd0"
	d <- density(x, n=n, from=1/n/2, to=1-1/n/2, bw=bw, na.rm=TRUE)
	d$y/sum(d$y)*d$n
}

allChrCancerHist <- sapply(s, apply, 2, at, n=n, simplify="array")
u <- split(data.frame(WGD=allChrAgg["WGD","time",isWgd]), droplevels(donor2type[sample2donor[names(finalSnv)]])[isWgd])
wgdCancerHist <- sapply(u, function(x) if(nrow(x)>0){at(x$WGD,n=n)}else{rep(0,n)}, simplify="array")
allChrCancerHist <- abind::abind(allChrCancerHist, All=sapply(sapply(s, as.matrix), at, n=n, simplify="array")/23*5, WGD=wgdCancerHist, along=2)

#' ### Per tumour type
#+ histTiming, fig.height=6, fig.width=6
prgn <- RColorBrewer::brewer.pal(11,"PRGn")
set1 <- RColorBrewer::brewer.pal(9,"Set1")
col <- colorRampPalette(set1[c(4,9,3)])(n)

p <- 0
v <- table(droplevels(donor2type[sample2donor[names(finalSnv)]]))
h <- (allChrCancerHist + p)  / rep(v + p, each=prod(dim(allChrCancerHist)[1:2]))
h <- aperm(h, c(2,3,1))

a <- colMeans(h[c("All","WGD"),,] * c(23/5,1)) %*% 1:n / asum(h* c(23/5,1), c(1,3))
o <- order(-a)
h <- h[,o,]
w <- v[o]>=15 & apply(h, 2, max) > 0.05*8/n
h <- h[,w,]

m <- 0.02
layout(matrix(1:prod(dim(h)[1:2]+1), ncol=dim(h)[1]+1, byrow=TRUE), height=c(rev(apply(h, 2, max))+m, 0.15), width=c(5, rep(1,dim(h)[1])))
par(mar=c(0.05,0.1,0,0.1), xpd=NA)
for(j in dim(h)[2]:0+1) for(i in 0:dim(h)[1]+1) {
		#if(all(h[i,j,]==0)) 
		if(i==1 & j !=1) {plot(NA,NA,xlab="",ylab="", xaxt="n",yaxt="n",xlim=c(0,1),ylim=c(0,1), bty="n")
			text(1,0,dimnames(h)[[2]][j-1],pos=2)
			next
		}
		if(j ==1 ){
			plot(NA,NA,xlab="",ylab="", xaxt="n",yaxt="n",xlim=c(0,1),ylim=c(0,1), bty="n")
			if(i==1) next
			text(0.5,1,dimnames(h)[[1]][i-1],pos=1)
			next
		}
		r <- c(0,max(h[,j-1,]+m))
		par(bty=if(i==2)"L" else "n")
		barplot(h[i-1,j-1,], ylim=r, width=1/n,space=0, col=rev(col), xaxt="n", yaxt="n", xlab="",ylab="", border=NA,xpd=TRUE, yaxs="i", xaxs="i", xlim=c(-0.5/n,1+0.5/n))
		axis(side=1, at=c(-0.5/n,1+0.5/n), labels=c("",""), tcl=-.1)
		if(i>1)
			#abline(v=0, col='grey')
		if(i==2){
			#abline(h=0, col='grey', lty=1)
			axis(side=2, at=c(0,0.05*8/n), labels=c("",""), tcl=-.1)
		}
	}
#dev.copy2pdf(file="histTiming.pdf",width=6, height=6, pointsize=8)


vv <- v[dimnames(h)[[2]]]
vv <- vv/sum(vv)

hh <- matrix(matrix(aperm(h, c(1,3,2)), ncol=length(vv)) %*% vv, nrow=nrow(h))
rownames(hh) <- rownames(h)

#' ### Pan-Can histograms
#+ histTimingPanCan, fig.height=2. fig.width=2
par(mar=c(3,3,1,1), mgp=c(2,.5,0), tcl=-0.5, bty="L", xpd=NA)
barplot(hh["WGD",], space=0, col=rev(col), xlab="Time [mutations]", ylab="Relative frequency", width=0.1, ylim=c(0,.065), yaxs='r', border=NA)
axis(side=1)
barplot(hh["All",], space=0, col=rev(col), xlab="Time [mutations]", ylab="Relative frequency", width=0.1, ylim=c(0,.065), yaxs='r')
axis(side=1)
#dev.copy2pdf(file="histTimingPanCan.pdf",width=2, height=2, pointsize=8)


#' # Real-time WGD & subclones
age <- clinicalData$donor_age_at_diagnosis
names(age) <- clinicalData$icgc_donor_id

typeNa <- gsub("\t","",strsplit("Bone-Cart
						Breast-LobularCa
						Breast-DCIS
						Lymph-NOS
						Myeloid-MDS
						Cervix-AdenoCa", "\n")[[1]])

#' ## Subclones
#' ### Prelim
#' Calculate effective genome size, i.e. time-averaged ploidy from mutation copy numbers
#+ effGenome
effGenome <- unlist(mclapply(finalSnv, function(vcf) {
					w <- info(vcf)$CLS!="subclonal"
					if(donor2type[sample2donor[meta(header(vcf))$META["ID",]]]=="Skin-Melanoma")
						w <- w & isDeaminationNoUV(vcf)
					else
						w <- w & isDeamination(vcf)
					2/avgWeights(vcf[na.omit(w)])
				}, mc.cores=MC_CORES))
names(effGenome) <- names(finalSnv)

subcloneDeam <- t(simplify2array(mclapply(finalSnv, function(vcf) {
							if(donor2type[sample2donor[meta(header(vcf))$META["ID",]]]=="Skin-Melanoma")
								w <- isDeaminationNoUV(vcf)
							else
								w <- isDeamination(vcf)
							p <- info(vcf)$pSub[w]; 
							c(sum(p, na.rm=TRUE), sum(1-p, na.rm=TRUE))})))

d <- droplevels(donor2type[sample2donor[names(finalSnv)]])
typesSubclones <- setdiff(levels(d), c(typeNa, names(which(table(d)<5))))

nClones <- sapply(finalClusters, nrow)

#' ### Mutation rates and age
#' Analyse relation to age, exclude hypermutators and samples with tumour in normal 1%. 
#+ timeSubcloneAge, fig.width=10, fig.height=10
deamRate <- cc <- list()
remove <- character()
par(mfrow=c(6,6), mar=c(3,3,2,1),mgp=c(2,.5,0), tcl=-0.25,cex=1, bty="L", xpd=FALSE, las=1, xpd=FALSE)
for(n in typesSubclones){
	i <- d==n
	tt0 <- subcloneDeam[i,]/cbind(finalPloidy[i], effGenome[i]) / cbind(nClones[i]-1, 1)/3 # 3Gb Haploid genome
	tt0[is.infinite(tt0)|is.nan(tt0)] <- 0
	yy <- rowSums(tt0)
	a <- age[sample2donor[names(finalSnv)[i]]]
	xx <- a
	r <- yy/xx 
	m <- median(r[TiN[names(xx)] <= 0.01 & ! is.na(TiN[names(xx)])],na.rm=TRUE)
	deamRate[[n]] <- r
	try({
				w <- (r-m)^2/m^2 <= 2^2 & TiN[names(xx)] <= 0.01 & ! is.na(TiN[names(xx)])
				remove <- c(remove, names(which(!w)))
				plot(xx, yy, bg=tissueColors[n], col=tissueBorder[n], pch=NA, log='', xlab="Age at diagnosis", ylab="SNVs/Gb", main=n, ylim=c(0,pmin(1000,max(yy, na.rm=TRUE))), xlim=c(0,max(age, na.rm=TRUE)),  cex.main=1)
				#par(xpd=NA)
				#segments(x0=0,y0=0, xx, yy, col=tissueLines[n], lty=tissueLty[n])
				points(xx, yy, bg=tissueColors[n], col=ifelse(w,tissueBorder[n], tissueColors[n]), pch=ifelse(w,21,4))
				abline(0, m, lty=3)
				#lines(c(x0,2*x0), c(0,1))
				#print(paste(n,cor(xx[w],yy[w],use='c'), cor(xx[w],tt0[w,] %*% c(0.2,1), use='c'), sep=": "))
				l <- lm(yy~xx, data=data.frame(yy=yy[w],xx=xx[w], row.names=as.character(seq_along(yy[w]))))
				f <- (summary(l))
				p <- predict(l,newdata=data.frame(xx=seq(0,100,1)), se=TRUE)
				polygon(c(seq(0,100,1),rev(seq(0,100,1))), c(p$fit - 2*p$se.fit, rev(p$fit+2*p$se.fit)), border=tissueBorder[n], col=paste0(tissueColors[n],"44"))
				v <- which(!is.na(xx[w]*yy[w]))
				cc[[n]] <- cbind(f$coefficients, f$cov.unscaled * f$sigma^2, coef(nnls::nnls(cbind(1,xx[w][v]), yy[w][v])))
			})
}
n <- names(deamRate)
q <- sapply(deamRate, function(r){
			m <- median(r[TiN[sample2donor[names(r)]] <= 0.01 & ! is.na(TiN[sample2donor[names(r)]])],na.rm=TRUE)
			w <- (r-m)^2/m^2 <= 2^2 & TiN[sample2donor[names(r)]] <= 0.01 & ! is.na(TiN[sample2donor[names(r)]])
			quantile(r[w], na.rm=TRUE)})
plot(sapply(deamRate, median, na.rm=TRUE), pch=NA , ylab="SNVs/Gb/yr", main="CpG>TpG rate", ylim=c(0, max(q)), cex.main=1, xaxt='n', xlab="Tumour type")
segments(seq_along(deamRate),q["0%",],seq_along(deamRate), q["100%",], col=tissueLines[n], lty=1)
points(sapply(deamRate, median, na.rm=TRUE), pch=21, col=tissueBorder[n], bg=tissueColors[n])

length(remove)

#' Pan-can
#+ timeSubcloneAgePancan, fig.width=2, fig.height=2
tt0 <- subcloneDeam/cbind(finalPloidy, effGenome) / cbind(nClones-1, 1)/3 # 3Gb Haploid genome
tt0[is.infinite(tt0)|is.nan(tt0)] <- 0
m <- sapply(deamRate, function(r){m <- median(r[TiN[sample2donor[names(r)]] <= 0.01 & ! is.na(TiN[sample2donor[names(r)]])],na.rm=TRUE)})
s <- rowSums(tt0)#/m[as.character(donor2type[sample2donor[names(finalSnv)]])] 
s[remove] <- NA
t <- donor2type[sample2donor[names(finalSnv)]]
x <- age[sample2donor[names(finalSnv)]]
plot(x,s, bg=tissueColors[t], pch=21, ylim=c(0,1000), col=tissueBorder[t], cex=tissueCex[t]*2/3, lwd=0.25, xlab="Age", ylab="SNVs/Gb")
p <- predict(loess(s~x), newdata=sort(x, na.last=NA), se=TRUE)
r <- function(x) c(x, rev(x))
polygon(r(sort(x, na.last=NA)), c(p$fit+2*p$se, rev(p$fit-2*p$se)), col="#00000044", border=NA)
lines(sort(x, na.last=NA),p$fit)
#s <- 12/8; dev.copy2pdf(file="timeSubcloneAgePancan.pdf", width=2*s, height=2*s, pointsize=8*s)


#for(nn in names(m)) abline(a=0,b=m[nn],col=tissueLines[nn], lty=tissueLty[nn])


#' Positive intercept?
a <- simplify2array(cc[!names(cc) %in% c("Myeloid-AML","Bone-Epith")])
all(a[1,1,] + 2*a[1,2,]>0 )

#' Positive slope?
all(na.omit(a[2,1,] + 2*a[2,2,]>0))

#+rateOffset, fig.width=2, fig.height=2
plot(a[1,1,], a[2,1,], col=tissueColors[dimnames(a)[[3]]], pch=NA, xlab="Offset", ylab="SNVs/Gb/yr")
segments(a[1,1,], a[2,1,] - a[2,2,],a[1,1,], a[2,1,]+a[2,2,], col=tissueLines[dimnames(a)[[3]]], pch=19)
segments(a[1,1,]-a[1,2,], a[2,1,], a[1,1,]+a[1,2,], a[2,1,], col=tissueLines[dimnames(a)[[3]]], pch=19)
points(a[1,1,], a[2,1,], pch=21, bg=tissueColors[dimnames(a)[[3]]], col=tissueLines[dimnames(a)[[3]]])
abline(h=0, lty=3)
abline(v=0, lty=3)

#' Fraction of mutations due to linear accumulation
#+fracLinear, fig.width=4, fig.height=2
par(mar=c(6,3,1,1))
ma <- sapply(split(age, donor2type[names(age)]), median, na.rm=TRUE)
fm <- pmax(a[2,7,],0)*ma[dimnames(a)[[3]]]/(pmax(0,a[2,7,])*ma[dimnames(a)[[3]]] + pmax(0,a[1,7,]))*100
o <- order(fm)
fmq <- sapply(names(fm), function(n){
			aa <- mvtnorm::rmvnorm(10000, mean=a[,1,n], sigma=a[,5:6,n] )
			aa <- aa[aa[,1]>=0 & aa[,2]>=0,]
			quantile(pmax(aa[,2],0)*ma[n]/(pmax(0,aa[,2])*ma[n] + pmax(0,aa[,1])), c(0.025, 0.975))
		}) *100
barplot(fm[o], col=tissueColors[dimnames(a)[[3]]][o], border=tissueLines[dimnames(a)[[3]]][o], las=2,names.arg=rep("",length(fm)) , ylab="Age-attributed mutations [%]") -> b
mg14::rotatedLabel(b, labels=names(fm[o]))
segments(b, fm[o], b, fmq[2,o], col=tissueLines[dimnames(a)[[3]]][o], lwd=2)
segments(b, fmq[1,o], b, fm[o], col=tissueBorder[dimnames(a)[[3]]][o], lwd=2)
abline(h=min(fmq[2,]))
abline(h=max(fmq[1,]))

#' ### Hierarchical Bayesian models of deamination rates
#' Prepare data
library(rstan)
y <- Reduce("c",deamRate)
y <- y[!names(y) %in% remove]
x <- age[sample2donor[names(y)]]
y <- y*x
t <- donor2type[sample2donor[names(y)]]
d <- data.frame(x,y,t)
df <- d
df <- df[rowSums(is.na(df))==0,]
tt <- model.matrix(~droplevels(t)-1, data=df)

data <- list(n = nrow(df),
		p = ncol(tt),
		y = df$y,
		x = tt * df$x,
		t = tt
)

#+ PCAWG-rates.chunk, echo=FALSE
read_chunk('./PCAWG-rates.stan', labels="PCAWG-rates.stan")

#' Model definition for stan
#+ PCAWG-rates.stan, eval=FALSE

#' Fit model
#+ PCAWG-stan
fit <- stan(
		file = "PCAWG-rates.stan",  # Stan program
		data = data,            # named list of data
		chains = 1,             # number of Markov chains
		warmup = 1000,          # number of warmup iterations per chain
		iter = 2000,            # total number of iterations per chain
		cores = 1,              # number of cores 
		refresh = 1000,         # show progress every 'refresh' iterations
		open_progress=FALSE,
		seed=42
)

#' Collect parameters
s <- summary(fit, pars=c("alpha","beta"))$summary
ab <- array(s, dim=c(33,2,10), dimnames=list(levels(droplevels(t)), c("a","b"), colnames(s)))

#' Summary plot
#+ deamAgeBayes, fig.width=2, fig.height=2
plot(x,y, bg=tissueColors[t], pch=21, ylim=c(0,1000), col=tissueBorder[t], cex=tissueCex[t]*2/3, lwd=0.25, xlab="Age", ylab="SNVs/Gb")
for(i in 1:nrow(ab))
	abline(ab[i,1,"50%"], ab[i,2,"50%"], col=tissueLines[levels(droplevels(t))[i]], lty=tissueLty[levels(droplevels(t))[i]])

#' Rate and offset
#+rateOffsetBayes, fig.width=2, fig.height=2
plot(ab[,1,"50%"], ab[,2,"50%"], col=tissueColors[dimnames(ab)[[1]]], pch=NA, xlab="Offset", ylab="SNVs/Gb/yr", xlim=range(ab[,1,c("2.5%","97.5%")]), ylim=range(ab[,2,c("2.5%","97.5%")]))
segments(ab[,1,"50%"], ab[,2,"2.5%"],ab[,1,"50%"], ab[,2,"97.5%"], col=tissueLines[dimnames(ab)[[1]]], pch=19)
segments(ab[,1,"2.5%"], ab[,2,"50%"],ab[,1,"97.5%"], ab[,2,"50%"], col=tissueLines[dimnames(ab)[[1]]], pch=19)
points(ab[,1,"50%"], ab[,2,"50%"], pch=21, bg=tissueColors[dimnames(ab)[[1]]], col=tissueLines[dimnames(ab)[[1]]])
abline(h=0, lty=3)
abline(v=0, lty=3)

a <- extract(fit, pars="alpha")$alpha
b <- extract(fit, pars="beta")$beta
colnames(a) <- colnames(b) <- levels(droplevels(t))

#' Overview
#+ timeSubcloneAgeBayes, fig.width=10, fig.height=10
par(mfrow=c(6,6), mar=c(3,3,2,1),mgp=c(2,.5,0), tcl=-0.25,cex=1, bty="L", xpd=FALSE, las=1, xpd=FALSE)
d <- droplevels(donor2type[sample2donor[names(finalSnv)]])
for(n in typesSubclones){
	i <- d==n
	tt0 <- subcloneDeam[i,]/cbind(finalPloidy[i], effGenome[i]) / cbind(nClones[i]-1, 1)/3 # 3Gb Haploid genome
	tt0[is.infinite(tt0)|is.nan(tt0)] <- 0
	yy <- rowSums(tt0)
	xx <- age[sample2donor[names(finalSnv)[i]]]
	r <- yy/xx 
	m <- median(r[TiN[names(xx)] <= 0.01 & ! is.na(TiN[names(xx)])],na.rm=TRUE)
	try({
				w <- !names(yy) %in% remove
				plot(xx[w], yy[w], bg=tissueColors[n], col=tissueBorder[n], pch=21, log='', xlab="Age at diagnosis", ylab="SNVs/Gb", main=n, ylim=c(0,pmin(1000,max(yy[w], na.rm=TRUE))), xlim=c(0,max(age, na.rm=TRUE)),  cex.main=1)
				#points(xx, yy, bg=tissueColors[n], col=ifelse(w,tissueBorder[n], tissueColors[n]), pch=ifelse(w,21,4))
				abline(0, m, lty=3)
				x0 <- seq(0,100,1)
				p <- apply(sapply(x0, function(x) a[,n] + b[,n]*x), 2, quantile, c(0.025,0.975), na.rm=TRUE)
				polygon(c(x0,rev(x0)), c(p["2.5%",], rev(p["97.5%",])), border=tissueBorder[n], col=paste0(tissueColors[n],"44"))
			})
}

#' Fraction of mutations due to linear accumulation
q <- sapply(colnames(a), function(n){
			w <- which(t==n & !is.na(y))
			f <- sapply(w, function(j) x[j] * b[,n] / (a[,n] + x[j] * b[,n]))
			quantile(rowMeans(f), c(0.025, 0.25, .5,.75,.975))
		})*100

qPanCan=quantile(rowMeans(do.call("cbind",sapply(colnames(a), function(n){
									w <- which(donor2type[sample2donor[names(y)]]==n & !is.na(y))
									f <- sapply(w, function(j) x[j] * b[,n] / (a[,n] + x[j] * b[,n]))
								}))),
		c(0.025, 0.25, .5,.75,.975))*100

#' Plot
#+fracLinearBayes, fig.width=4, fig.height=2
par(mar=c(6,3,1,1))
o <- order(q["50%",])
barplot(q["50%",o], col=tissueColors[colnames(q)][o], border=tissueLines[colnames(q)][o], las=2,names.arg=rep("",length(q["50%",])) , ylab="Age-attributed mutations [%]", ylim=c(0,100)) -> b
mg14::rotatedLabel(b, labels=names(q["50%",o]))
segments(b, q["50%",][o], b, q["97.5%",o], col=tissueLines[colnames(q)][o], lwd=2)
segments(b, q["2.5%",o], b, q["50%",][o], col=tissueBorder[colnames(q)][o], lwd=2)
abline(h=min(q["97.5%",]), lty=3)
abline(h=max(q["2.5%",]), lty=3)
#abline(h=qPanCan["50%"], lty=4)


#' ### Timing
#' Acceleration values to simulate
accel <- c(1,2.5,5,7.5,10)
names(accel) <- paste0(accel, "x")

#' The actual timing
#+ timeSubclones, warning=FALSE
set.seed(42)
timeSubclones <- sapply(typesSubclones, function(l) {
			i <- d==l
			tt0 <- subcloneDeam[i,]/cbind(finalPloidy[i], effGenome[i]) / cbind(nClones[i]-1, 1)
			resB <- sapply(1:1000, function(foo){ ## Assess the impact of Poisson fluctuations on numbers
						tt <- matrix(rpois(length(tt0), lambda=tt0), ncol=ncol(tt0))
						res <- sapply(accel, function(a)  tt[,1]/a/rowSums(tt/rep(c(a,1), each=nrow(tt)))) * age[sample2donor[names(finalSnv)[i]]]
						colnames(res) <- paste0(accel, "x")
						#res[res==0] <- NA
						res}, simplify='array')
			res <- sapply(accel, function(a)  tt0[,1]/a/rowSums(tt0/rep(c(a,1), each=nrow(tt0)))) * age[sample2donor[names(finalSnv)[i]]]
			colnames(res) <- paste0(accel, "x")	
			resCI <- apply(resB,1:2, quantile, c(0.025,0.975), na.rm=TRUE)
			arr <- abind::abind(res, resCI, along=1)
			rownames(arr)[1] <- "hat"
			arr <- aperm(arr, c(2,1,3))
			tt0[is.infinite(tt0)|is.nan(tt0)] <- 0
			r <- which(rowSums(subcloneDeam[i,]) < 50 ) ## Exclude samples with less than 50 subs 
			arr[r,,] <- NA
			return(arr)
		})

#' Plot
#+ realTimeSubclone, fig.width=6, fig.height=2.225
u <- setdiff(names(finalSnv)[uniqueSamples], remove)
par( mar=c(7,3,1,1), mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
qSubclone <- sapply(timeSubclones, function(x) apply(x[,"hat",][rownames(x)%in%u,,drop=FALSE], 2, quantile, c(0.05,0.25,0.5,0.75,0.95), na.rm=TRUE), simplify='array')
a <- "5x"
tSubclonesByType <- sapply(timeSubclones, function(x) x[,,a][setdiff(rownames(x),remove), 1:3, drop=FALSE])
tSubclonesByType[["Ovary-AdenoCa"]] <- timeSubclones[["Ovary-AdenoCa"]][rownames(tSubclonesByType[["Ovary-AdenoCa"]]),,"7.5x"]
m <- qSubclone["50%",a,]#t[1,3,]
m["Ovary-AdenoCa"] <- qSubclone["50%","7.5x","Ovary-AdenoCa"]
m[sapply(tSubclonesByType, function(x) sum(!is.na(x[,1]))) < 5] <- NA
o <- order(m, na.last=NA)
plot(NA,NA, xlim=c(0.5,length(m[o])), ylab="Years before diagnosis", xlab="", xaxt="n", yaxs="i", ylim=c(0,30))
x <- seq_along(m[o])
mg14::rotatedLabel(x, labels=names(sort(m)))
for(i in seq_along(o))try({
				n <- names(m)[o[i]]
				f <- function(x) x/max(abs(x))
				a <- if(n== "Ovary-AdenoCa") "7.5x" else "5x" 
				bwd <- 0.8/2
				j <- if(length(na.omit(tSubclonesByType[[o[i]]][,"hat"]))>1) f(mg14::violinJitter(na.omit(tSubclonesByType[[o[i]]][,"hat"]))$y)/4 + i else i
				tpy <- 2
				segments(j, na.omit(tSubclonesByType[[o[i]]][,"97.5%"]), j, na.omit(tSubclonesByType[[o[i]]][,"2.5%"]), col=mg14::colTrans(tissueLines[n],tpy))
				points(j, na.omit(tSubclonesByType[[o[i]]][,"hat"]), pch=21, col=mg14::colTrans(tissueBorder[n], tpy), bg=mg14::colTrans(tissueColors[n],tpy), 
						cex=tissueCex[n]*2/3, lwd=1)
				rect(i-bwd,qSubclone["25%",a,n],i+bwd,qSubclone["75%",a,n], border=tissueLines[n],  col=paste(tissueColors[n],"44", sep=""))
				segments(i-bwd,qSubclone["50%",a,n],i+bwd,qSubclone["50%",a,n],col=tissueLines[n], lwd=2)
				segments(i,qSubclone["75%",a,n],i,qSubclone["95%",a,n],col=tissueLines[n], lwd=1.5)
				segments(i,qSubclone["5%",a,n],i,qSubclone["25%",a,n],col=tissueLines[n], lwd=1.5)
				f <- function(x) x/max(abs(x))
			})

#par(xpd=TRUE)
#s <- 12/8
#dev.copy2pdf(file="realTimeSubclone.pdf", width=6*s, height=3.5*3/5*s, pointsize=8*s)

sapply(timeSubclones, nrow)

#' Numbers per decade
yy <- do.call("rbind",tSubclonesByType)
yy <- yy[setdiff(rownames(yy), remove),"hat"]
table(cut(yy, seq(0,60,10)))

#' ## WGD
#' ### Functions
wgdTime <- function(vcf, bb, clusters, purity){
	w <- which(info(vcf)$MajCN==2 & info(vcf)$MinCN==2& sapply(info(vcf)$CNID, length)==1 & isDeamination(vcf))
	if(donor2type[sample2donor[meta(header(vcf))$META["ID",]]]=="Skin-Melanoma")
		w <- intersect(w, which(isDeaminationNoUV(vcf)))
	v <- vcf[w]
	if(nrow(v)<=100) return(NULL)
	seqnames(rowRanges(v)) <- factor(rep(1, nrow(v)), levels=seqlevels(v))
	b <- GRanges(1, IRanges(1,max(end(v))), copy_number=4, major_cn=2, minor_cn=2, clonal_frequency=purity)
	computeMutCn(v, b, clusters, purity, isWgd=TRUE, n.boot=10)
}

#+ finalWgdParam
finalWgdParam <- mclapply(names(finalSnv)[isWgd], function(ID){
			wgdTime(finalSnv[[ID]], finalBB[[ID]], clusters=finalClusters[[ID]], purity=finalPurity[ID])
		},  mc.cores=MC_CORES)

finalWgdPi <- sapply(finalWgdParam[!sapply(finalWgdParam, is.null)], function(x) {
			pi <- x$P[[1]][,"P.m.sX"] * x$P[[1]][,"pi.s"]
			pi.up <- (x$P[[1]][,"P.m.sX.up"] * x$P[[1]][,"pi.s"])[1:2]
			pi.lo <- (x$P[[1]][,"P.m.sX.lo"] * x$P[[1]][,"pi.s"])[1:2] 
			pi.clonal <- pi[1:2]
			pi.subclonal <- mean(pi[!x$P[[1]][,"clonalFlag"]]) # or sum for linear progression
			if(length(pi.subclonal)==0 | is.na(pi.subclonal) | is.infinite(pi.subclonal)) pi.subclonal <- 0 
			r <- rbind(cbind(hat=pi.clonal, lo=pi.lo, up=pi.up), pi.subclonal)
			r[1,2:3] <- r[1,3:2]
			r
		}, simplify='array')
dimnames(finalWgdPi)[[1]] <- c("clonal.1","clonal.2","sub")

correctAccel <- function(pi, ta, a){
	t0 <- 2*pi[2]/(2*pi[2]+pi[1]) ##  no acc
	t1 <- t0 + (1-t0) *(a-1)/a*ta #acc before dup
	t2 <- t0 * (ta + a*(1-ta)) ## after
	tWgdClonal <- pmin(t1, t2) # as fraction of clonal
	aEffClonal <- ta + (1-ta)*a # effective rate, avg over clonal
	gEffClonal <- 2*tWgdClonal + 4*(1-tWgdClonal) # effective genome size
	piClonal <- sum(pi[1:2])
	piSub <- sum(pi[-2:-1])
	tSub <- piSub / 4 / a
	tClonal <- piClonal / gEffClonal/ aEffClonal
	tClonal <- tClonal / (tClonal + tSub)
	return(c(tWgd=tWgdClonal * tClonal, tClonal=tClonal))
}

correctAccelRand <- function(pi, ta=seq(0.8,1,0.01), a=seq(1,10,1)){
	x <- sapply(ta, function(taa) sapply(a, function(aa) correctAccel(pi, taa, aa)), simplify='array')
}

#' ### Timing
foo <- apply(finalWgdPi[,,], 2:3,  function(x) correctAccelRand(x, a=accel))
finalWgdPiAdj <- array(foo, dim=c(2, length(accel), length(eval(formals(correctAccelRand)$ta)), dim(foo)[-1]))
dimnames(finalWgdPiAdj)[[1]] <- c('t.WGD','t.subclonal')
dimnames(finalWgdPiAdj)[[2]] <- paste0(accel, "x")
dimnames(finalWgdPiAdj)[[4]] <- dimnames(finalWgdPi)[[2]]
dimnames(finalWgdPiAdj)[[5]] <- names(finalBB)[isWgd][!sapply(finalWgdParam, is.null)]

#finalWgdPiAdj <- sapply(1:ncol(finalWgdPi), function(i) correctAccelRand(finalWgdPi[,i], a=accel), simplify='array')
#dimnames(finalWgdPiAdj)[[2]] <- paste0(accel, "x")

n <- dimnames(finalWgdPiAdj)[[5]]
finalWgdTime <- finalWgdPiAdj[,,,,n] * rep(age[sample2donor[n]], each=2)


d <- droplevels(donor2type[sample2donor[n]])
s <- setdiff(levels(d), c(typeNa, names(which(table(d)<3))))
timeWgd <- sapply(s, function(l) {
			i <- d==l & ! n %in% c(rownames(purityPloidy)[purityPloidy$wgd_uncertain])#, names(which(q5 > 0.1)))
			a <- (1-finalWgdPiAdj[1,,,,n][,,,i]) * rep(age[sample2donor[n]][i], each = prod(dim(finalWgdPiAdj)[2:4]))
			m <- aperm(mg14:::asum(a, 2)/dim(a)[2])#sum(!is.na(age[sample2donor[n]][i]))
			rownames(m) <- n[i]
			m
			#quantile(m, c(.025,.25,.5,.75,.975), na.rm=TRUE)
			#rowMeans(m)
			#apply(m, 1, median)
		}, simplify=FALSE)

#+ realTimeWgd, fig.height=3, fig.width=4
par( mar=c(7,3,1,1), mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
u <- setdiff(names(finalSnv)[uniqueSamples], remove)
qWgd <- sapply(timeWgd, function(x) apply(x[rownames(x) %in% u,"hat",], 2, quantile, c(0.25,0.5,0.75), na.rm=TRUE), simplify='array')
nWgd <- sapply(timeWgd, function(x) sum(x[rownames(x) %in% u,"hat","1x"]!=0, na.rm=TRUE))

a <- "5x"
m <- qWgd['75%',a,]#t[1,3,]
m["Ovary-AdenoCa"] <- qWgd["75%","7.5x","Ovary-AdenoCa"]
o <- order(m, na.last=NA)
x <- seq_along(m[o])
tWgdByType <- sapply(timeWgd, function(yy) yy[setdiff(rownames(yy), remove),1:3,a])
tWgdByType[["Ovary-AdenoCa"]] <- timeWgd[["Ovary-AdenoCa"]][setdiff(rownames(timeWgd[["Ovary-AdenoCa"]]), remove),,"7.5x"]
plot(NA,NA, xlim=c(0.5,length(m[o])), ylim=c(0,max(do.call('rbind',tWgdByType)[,1], na.rm=TRUE)+5), ylab="Years before diagnosis", xlab="", xaxt="n", yaxs="i")
mg14::rotatedLabel(x, labels=names(sort(m)))
for(i in seq_along(o)){
	n <- names(m)[o[i]]
	f <- function(x) x/max(abs(x))
	a <- if(n== "Ovary-AdenoCa") "7.5x" else "5x" 
	j <- f(mg14::violinJitter(na.omit(tWgdByType[[o[i]]][,"hat"]))$y)/4 + i
	tpy <- 2
	segments(j, na.omit(tWgdByType[[o[i]]][,"up"]), j, na.omit(tWgdByType[[o[i]]][,"lo"]), col=mg14::colTrans(tissueLines[n],tpy))
	points(j, na.omit(tWgdByType[[o[i]]][,"hat"]), pch=21, col=mg14::colTrans(tissueBorder[n], tpy), bg=mg14::colTrans(tissueColors[n],tpy), 
			cex=tissueCex[n]*2/3, lwd=1)
	bwd <- 0.8/2
	rect(i-bwd,qWgd[1,a,n],i+bwd,qWgd[3,a,n], border=tissueLines[n],  col=paste0(tissueColors[n],"44"))
	segments(i-bwd,qWgd[2,a,n],i+bwd,qWgd[2,a,n],col=tissueLines[n], lwd=2)
	#d <- density(na.omit(y[[o[i]]][,"hat"]))
	#polygon(c(d$y/max(d$y)/2.5+i, rev(-d$y/max(d$y)/2.5+i)),c(d$x, rev(d$x)), border=tissueLines[names(m)[o[i]]],  col=mg14::colTrans(tissueColors[names(m)[o[i]]],3))
#	segments(j, na.omit(y[[o[i]]][,"up"]), j, na.omit(y[[o[i]]][,"lo"]), col=tissueColors[names(m)[o[i]]], lwd=2)
#	points(j, na.omit(y[[o[i]]][,"hat"]), pch=16, col="white", cex=0.5)
}
par(xpd=TRUE)
#s <- 12/8
#dev.copy2pdf(file="realTimeWgd.pdf", width=4*s, height=3.5*s, pointsize=8*s)

#' Numbers per decade
yy <- do.call("rbind",tWgdByType)
yy <- yy[setdiff(rownames(yy), remove),"hat"]
table(cut(yy, seq(0,60,10)))

#' #### Assessment of absolute mutation counts
nDeam22 <- sapply(finalWgdParam, function(x) if(!is.null(x$D)) nrow(x$D) else NA)
names(nDeam22) <- names(finalSnv)[isWgd]
w22 <- sapply(finalBB[isWgd], function(bb) {
			w <- bb$major_cn==2 & bb$minor_cn==2 & !duplicated(bb)
			sum(as.numeric(width(bb)[w]), na.rm=TRUE)})
nDeam22 <- nDeam22/w22*1e9

t0 <- 2*finalWgdPi["clonal.2","hat",]/( 2*finalWgdPi["clonal.2","hat",] +  finalWgdPi["clonal.1","hat",])
names(t0) <- dimnames(finalWgdPiAdj)[[5]]

#+ nDeam22Time, fig.width=2, fig.height=2
d <- nDeam22 *  t(sapply(finalWgdParam, function(x) if(!is.null(x$P)) x$P[[1]][1:2,"P.m.sX"] else c(NA,NA)))
rownames(d) <- names(finalSnv)[isWgd]
t0 <- 2*finalWgdPi["clonal.2","hat",]/( 2*finalWgdPi["clonal.2","hat",] +  finalWgdPi["clonal.1","hat",])
names(t0) <- dimnames(finalWgdPiAdj)[[5]]
d[rownames(d) %in% remove] <- NA
y <- d[names(t0),]/6
x <- t0
t <- donor2type[sample2donor[names(t0)]]

plot(x,y[,2], bg=tissueColors[t], pch=21,  col=tissueBorder[t], cex=tissueCex[t]*1, lwd=0.5, xlab="Time", ylab="Early SNVs/Gb", log='')
#points(x,y[,1], bg=tissueColors[t], pch=21,  col=tissueBorder[t], cex=tissueCex[t]*1, lwd=0.5)
p <- predict(loess(y[,2]~x, span=1), newdata=sort(x, na.last=NA), se=TRUE)
r <- function(x) c(x, rev(x))
polygon(r(sort(x, na.last=NA)), c(p$fit+2*p$se, rev(p$fit-2*p$se)), col="#00000044", border=NA)
lines(sort(x, na.last=NA),p$fit)

plot(x,y[,1]+y[,2], bg=tissueColors[t], pch=21,  col=tissueBorder[t], cex=tissueCex[t]*1, lwd=0.5, xlab="Time", ylab="Total SNVs/Gb", log='')
p <- predict(loess(rowSums(y)~x, span=1), newdata=sort(x, na.last=NA), se=TRUE)
r <- function(x) c(x, rev(x))
polygon(r(sort(x, na.last=NA)), c(p$fit+2*p$se, rev(p$fit-2*p$se)), col="#00000044", border=NA)
lines(sort(x, na.last=NA),p$fit)

#s <- 12/8; dev.copy2pdf(file="nDeam22Time.pdf", width=2*s, height=2*s, pointsize=8*s)


#' #### Acceleration adjustment relative to lowest quintile.
#+ accelRelWgd, fig.height=2, fig.width=2
accelRelWgd <- sapply(names(timeWgd), function(n) {
			r <- deamRate[[n]]
			x <- timeWgd[[n]]
			r0 <- quantile(r[!names(r) %in% remove], 0.20, na.rm=TRUE)
			a <- cut(pmin(pmax(1,(r/r0-0.9)/0.1),10), c(0,1.5,3.75,6.25,8.75,20), labels=c("1x","2.5x","5x","7.5x","10x"))
			print(table(a))
			names(a) <- names(r)
			ta <- sapply(rownames(x), function(ss) x[ss, "hat",a[ss]])
			ta
		})
par(mar=c(3,3,1,1),mgp=c(2,.5,0), tcl=-0.25,cex=1, bty="L", xpd=FALSE, las=1)
x <- Reduce("c",accelRelWgd)
y <- Reduce("c",sapply(timeWgd, function(x) x[,"hat","5x"]))

t <- donor2type[sample2donor[names(x)]]
plot(y+runif(length(y)),x+runif(length(y)), pch=21, bg=tissueColors[t], col=tissueBorder[t], cex=tissueCex[t]*1, lwd=0.5, xlab="Time (constant acceleration)", ylab="Time (sample-specific accel.)", log='')
#s <- 12/8; dev.copy2pdf(file="accelRelWgd.pdf", width=2*s, height=2*s, pointsize=8*s)

#' Quantiles
#+ qAccelRelWgd, fig.height=2, fig.width=2
par(mar=c(3,3,1,1),mgp=c(2,.5,0), tcl=-0.25,cex=1, bty="L", xpd=FALSE, las=1)
u <- setdiff(names(finalSnv)[uniqueSamples], remove)
qAccelRelWgd <- sapply(accelRelWgd, function(x){
			quantile(x[names(x) %in% u], na.rm=TRUE)
		})
t <- colnames(qAccelRelWgd)
plot(qWgd["50%","5x",], qAccelRelWgd["50%",], pch=NA,  xlab="Time [years], constant acceleration", ylab="Time [years], sample-specific accel.", xlim=c(0,30), ylim=c(0,30))
abline(0,1, lty=3)
segments(qWgd["25%","5x",], qAccelRelWgd["50%",],qWgd["75%","5x",], qAccelRelWgd["50%",], lty=tissueLty[t], col=tissueLines[t])
segments(qWgd["50%","5x",], qAccelRelWgd["25%",],qWgd["50%","5x",], qAccelRelWgd["75%",], lty=tissueLty[t], col=tissueLines[t])
points(qWgd["50%","5x",], qAccelRelWgd["50%",], bg=tissueColors[t], pch=21,  col=tissueBorder[t], cex=tissueCex[t], lwd=0.5)
#s <- 12/8; dev.copy2pdf(file="qAccelRelWgd.pdf", width=2*s, height=2*s, pointsize=8*s)



#' Mutations per year vs time
#+ mutYearTime, fig.height=8, fig.width=8
par(mfrow=c(5,5), mar=c(3,3,2,1),mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
for(n in names(tWgdByType)){
	a <- age[sample2donor[rownames(tWgdByType[[n]])]]
	yy <- nDeam22[rownames(tWgdByType[[n]])]/a
	xx <- 1-t0[rownames(tWgdByType[[n]])]#y[[n]][,"hat"]/a
	try({
				l <- lm(yy ~ xx)
				x0 <- -l$coef[2]/l$coef[1]
				#print(x0)
				plot(xx, yy, bg=tissueColors[n], col=tissueBorder[n],, pch=21, log='', xlab="1-time [mutations]", ylab="SNVs/Gb/yr", main=n, ylim=c(0,max(yy, na.rm=TRUE)), xlim=c(0,max(xx, na.rm=TRUE)), cex.main=1)
				#abline(l, lty=3)
				#abline(l$coef[1], l$coef[1])
				m <- median(yy/(1+xx),na.rm=TRUE)
				abline(m, m, lty=3)
				#lines(c(x0,2*x0), c(0,1))
			})
}

#' Age at diagnosis
#+ mutAgeWgd, fig.height=8, fig.width=8
par(mfrow=c(5,5), mar=c(3,3,2,1),mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1, xpd=FALSE)
deamRateWgd <- list()
for(n in names(tWgdByType)){
	a <- age[sample2donor[rownames(tWgdByType[[n]])]]
	yy <- nDeam22[rownames(tWgdByType[[n]])]/(2-t0[rownames(tWgdByType[[n]])])
	xx <- a
	r <- yy/xx 
	m <- median(r,na.rm=TRUE)
	deamRateWgd[[n]] <- r
	try({
				w <- (r-m)^2/m^2 <= 2^2 
				plot(xx, yy, bg=tissueColors[n], col=tissueBorder[n], pch=NA, log='', xlab="Age at diagnosis", ylab="SNVs/Gb", main=n, ylim=c(0,max(yy, na.rm=TRUE)), xlim=c(0,max(age, na.rm=TRUE)),  cex.main=1)
				par(xpd=NA)
				segments(x0=0,y0=0, xx, yy, col=tissueLines[n], lty=tissueLty[n])
				points(xx, yy, bg=tissueColors[n], col=ifelse(w,tissueBorder[n], tissueColors[n]), pch=ifelse(w,21,4))
				par(xpd=FALSE)
				#abline(l, lty=3)
				#abline(l$coef[1], l$coef[1])
				abline(0, m, lty=3)
				#lines(c(x0,2*x0), c(0,1))
			})
}
n <- names(deamRateWgd)
q <- sapply(deamRateWgd, function(r){
			m <- median(r,na.rm=TRUE)
			w <- (r-m)^2/m^2 <= 2^2 
			range(r[w], na.rm=TRUE)})
plot(sapply(deamRateWgd, median, na.rm=TRUE), pch=NA , ylab="SNVs/Gb/yr", main="CpG>TpG rate", ylim=c(0, max(q)), cex.main=1, xaxt='n', xlab="Tumour type")
segments(seq_along(deamRateWgd),q[1,],seq_along(deamRateWgd), q[2,], col=tissueLines[n], lty=1)
points(sapply(deamRateWgd, median, na.rm=TRUE), pch=21, col=tissueBorder[n], bg=tissueColors[n])

#' Conceptual plot
#+ concept, fig.height=2, fig.width=2
#par(mfrow=c(1,1), mar=c(3,3,1,1), mgp=c(2,0.5,0), bty="L")
x0 <- 70
y0 <- 1
a <- 5
plot(x0,y0, xlab="Time [yr]", ylab="Time [fraction of mutations]", pch=19, xlim=c(0,80), ylim=c(0,y0*1.1))
r <- seq(0.66,1,0.001)
yy <- c(0,y0*r/(a*(1-r)+r))
xx <- c(0,r)*x0
polygon(xx, yy, col="grey", border=NA)
for(i in seq(51, 341,30))
	lines(c(0,xx[i],x0), c(0,yy[i],y0), col='darkgrey')
t <- 3/5*y0
ta <- function(t0, a, ta){
	t1 <- t0 + (1-t0) *(a-1)/a*ta #acc before dup
	t2 <- t0 * (ta + a*(1-ta)) ## after
	tf <- pmin(t1, t2) # as fraction of clonal
	return(tf)
}
tf <- sapply(r, function(rr) ta(t/y0, a, rr))
tmax <- xx[which.min(abs(yy-t))]
tmin <- t/ y0*x0
lines(c(x0,x0,0), c(0,y0,y0), lty=3)
lines(c(tmax,tmax,0), c(0,t,t), lty=2)
lines(c(tmin,tmin,0), c(0,t,t), lty=2)
#d <- density(tf*x, from=tmin, to=tmax)
lines(quantile(tf, c(0.025,0.975))*x0, c(0,0), lwd=2)
points(c(tmin,median(tf)*x0), c(0,0), pch=c(1,19))
mtext(side=1, at=x0, text="Diagnosis", line=2)
text(x=0, y=t, labels="WGD", pos=4, adj=c(0,1))
#lines(d$x,d$y/max(d$y)*50)
#s <- 12/8; dev.copy2pdf(file="concept.pdf", width=2*s, height=2*s, pointsize=8*s)




#+ realTimeWgdAccel, fig.height=2, fig.width=2
par(mar=c(3,3,1,1), mgp=c(2,0.5,0), tcl=-0.25, bty="L")
plot(accel, qWgd["50%",,1], type='l', lty=0, ylim=c(0,30), xlab= "5meC rate acceleration", ylab="Median occurrence [years]", yaxs="i", xaxt="n")
axis(side=1, at=accel)
for(j in 1:dim(qWgd)[3]) lines(accel, qWgd["50%",,j], type='l', col=tissueLines[dimnames(qWgd)[[3]][j]], 
			lty=ifelse(nWgd[dimnames(qWgd)[[3]][j]]<=4, 3, tissueLty[dimnames(qWgd)[[3]][j]]))
#s <- 12/8; dev.copy2pdf(file="realTimeWgdAccel.pdf", width=2*s, height=2*s, pointsize=8*s)



#+ realTimeSubcloneWgdScatter
par( mar=c(4,3,1,1), mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
a <- "5x"
plot(qSubclone["50%",a,dimnames(qWgd)[[3]]], qWgd["50%",a,], col=tissueColors[dimnames(qWgd)[[3]]], pch=16, cex=2, xlab="Median time subclones", ylab="Median time WGD", xlim=c(0,5), ylim=c(0,10))

#+ realTimeSubcloneWgd, fig.width=2.5, fig.height=3.5
par( mar=c(3,3,3,10), mgp=c(2,.5,0), tcl=-0.25,cex=1, bty="n", xpd=FALSE, las=1)
w <- "50%"
x <- qSubclone[w,a,]
x["Ovary-AdenoCa"] <- qSubclone[w,"7.5x","Ovary-AdenoCa"]
y <- qWgd[w,a,]
y["Ovary-AdenoCa"] <- qWgd[w,"7.5x","Ovary-AdenoCa"]
plot(c(rep(1, dim(qSubclone)[3]), rep(2, each=dim(qWgd)[3])), c(x,y), bg=tissueColors[c(dimnames(qSubclone)[[3]], dimnames(qWgd)[[3]])], pch=21, cex=1, xaxt="n", ylab="Years before diagnosis", xlab="", xlim=c(0.5,2.5), ylim=c(0, max(y, na.rm=TRUE)))
segments(rep(1, each=dim(qWgd)[3]), x[dimnames(qWgd)[[3]]], rep(2, each=dim(qWgd)[3]), y,col=tissueLines[dimnames(qWgd)[[3]]], lty= ifelse(nWgd <= 5, 3, tissueLty[dimnames(qWgd)[[3]]]))
o <- order(y, na.last=NA)
y0 <- y[o]
y1 <- mg14:::mindist(y[o], diff(par('usr')[3:4])/30)
par(xpd=NA)
mtext(dimnames(qWgd)[[3]][o], at=y1, side=4 )
segments(2.1,y0,2.2,y0)
segments(2.2,y0,2.3,y1)
segments(2.3,y1,2.4,y1)
mg14::rotatedLabel(1:2, labels=c("Subclones","WGD"))
#segments(2.1,y0,2.3,y1)

#par(xpd=TRUE)
#s <- 12/8
#dev.copy2pdf(file="realTimeSubcloneWgd.pdf", width=2.5*s, height=3.5*s, pointsize=8*s)

#save(qWgd, qSubclone, timeWgd, timeSubclones, file=paste0(Sys.Date(),"-realTimeWgdAndSubclones.RData"))

#' # Session
#' ## Objects
l <- ls()
data.frame(variable=l, Reduce("rbind",lapply(l, function(x) data.frame(class=class(get(x)), size=format(object.size(get(x)), units="auto")))))
#' ## Packages
sessionInfo()
devtools::session_info()

#' # Other code
#' All code is available at github.com/gerstung-lab/PCAWG-11
#+ additionalCode, cache=FALSE, echo=FALSE, eval=TRUE
read_chunk('./MutationTime.R', labels="MutationTimer")
read_chunk('./PCAWG-functions.R', labels="PCAWG-functions")
read_chunk('./VCF-annotate.R', labels="VCF-annotate")

#' ## MutationTime.R
#' See https://gist.github.com/mg14/7a8e1aa28cb9ade7e376acdbd2364790
#+ MutationTimer, eval=FALSE

#' ## PCAWG-functions.R
#' All basic functions
#+ PCAWG-functions, eval=FALSE

#' ## VCF-annotate.R
#+ VCF-annotate, eval=FALSE