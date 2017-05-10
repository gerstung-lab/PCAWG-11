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

#+ evalOff, echo=FALSE
opts_chunk$set(eval=FALSE)

#' ## Load data
#' ### SNV and MNV
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_010/snv_mnv"
finalSnv <- list()
j <- 1
for(f in dir(p, pattern="*.vcf.RData", full.names=TRUE)){
	if(j %% 10 ==0) print(j); j <- j+1
	load(f)
	finalSnv[[f]] <- vcf
}
names(finalSnv) <- sub(".conse.+","",dir(p, pattern="*.vcf.RData", full.names=FALSE))

#' ### Copy number profiles
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_010/cn"
finalBB <- list()
for( f in dir(p, pattern="*.bb_granges.RData", full.names=TRUE)){
	load(f)
	colnames(mcols(bb)) <- sub("star.1","time.star",colnames(mcols(bb)) ) # Fix naming problem
	finalBB[[f]] <- bb
}
names(finalBB) <- sub(".conse.+","",dir(p, pattern="*.bb_granges.RData", full.names=FALSE))

#' ### Indel
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_010/indel"
finalIndel <- list()
for( f in dir(p, pattern="*.vcf.RData", full.names=TRUE)){
	load(f)
	finalIndel[[f]] <- vcfIndel
}
names(finalIndel) <- sub(".conse.+","",dir(p, pattern="*.vcf.RData", full.names=FALSE))

#' ### Clusters and purity
finalClusters <- list()
finalPurity <- numeric()
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_010/clusters"
for( f in dir(p, pattern="*.RData", full.names=TRUE)){
	load(f)
	finalClusters[[f]] <- clusters
	finalPurity[f] <- purity
}
names(finalClusters) <- names(finalPurity) <- sub(".conse.+","",dir(p, pattern="*.RData", full.names=FALSE))


#' ## Update drivers
#' ### Update VCF
for(i in seq_along(finalSnv)){
	info(finalSnv[[i]])$DG <- matchDrivers(finalSnv[[i]], finalDrivers)
	info(finalIndel[[i]])$DG <- matchDrivers(finalIndel[[i]], finalDrivers)
	if(i %% 10 ==0) print(i); i <- i+1
}

#' ### Update finalDrivers
finalDriversAnnotated <- finalDrivers
finalDriversAnnotated$sample <- drivers$sample
finalDriversAnnotated$mut_type <- drivers$mut_type
d <- info(finalSnv[[3]])[seq_along(finalDriversAnnotated),19:33]
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


#' ## Load graylisted data
#' ### SNV and MNV
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_010/graylist/snv_mnv"
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
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_010/graylist/cn"
finalBBGray <- list()
for( f in dir(p, pattern="*.bb_granges.RData", full.names=TRUE)){
	load(f)
	colnames(mcols(bb)) <- sub("star.1","time.star",colnames(mcols(bb)) ) # Fix naming problem
	finalBBGray[[f]] <- bb
}
names(finalBBGray) <- sub(".conse.+","",dir(p, pattern="*.bb_granges.RData", full.names=FALSE))
finalBB[names(finalBBGray)] <- finalBBGray


#' ### Indel
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_010/graylist/indel"
finalIndelGray <- list()
for( f in dir(p, pattern="*.vcf.RData", full.names=TRUE)){
	load(f)
	finalIndelGray[[f]] <- vcfIndel
}
names(finalIndelGray) <- sub(".conse.+","",dir(p, pattern="*.vcf.RData", full.names=FALSE))
finalIndel[names(finalIndelGray)] <- finalIndelGray


#' ### Clusters and purity
finalClustersGray <- list()
finalPurityGray <- numeric()
p <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/annotated_010/graylist/clusters"
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

#' ## QC
#+ QC
q1 <- sapply(finalSnv, function(vcf) mean(abs(0.5- info(vcf)$pMutCNTail) > 0.495 , na.rm=TRUE))
q5 <- sapply(finalSnv, function(vcf) mean(abs(0.5- info(vcf)$pMutCNTail) > 0.475 , na.rm=TRUE))

par(mfrow=c(1,1))
boxplot(1-q5 ~ donor2type[sample2donor[names(finalSnv)]], las=2, ylab="Fraction of data inside theoretical 95% CI")
abline(h=0.95, lty=3)

#+ QQplots, fig.width=4, fig.height=4
#pdf("QQplots.pdf", 4,4, pointsize=8)
for(i in seq_along(finalSnv)){
	n <- nrow(finalSnv[[i]])
	qqnorm(qnorm(info(finalSnv[[i]])$pMutCNTail[sample(1:n, min(1e4,n))]), main=paste(substr(names(finalSnv)[i],1,8), "Q5 =", signif(q5[i],2), ", Q1 =", signif(q1[i],2)), xlim=c(-5,5), ylim=c(-5,5), pch=16)
	abline(0,1, col='red')
}
#dev.off()


#' ## Distribution of Mutations
#' ### MAP genotypes
finalGenotypesSnv <- simplify2array(mclapply(finalSnv[whiteList], getGenotype, mc.cores=2, useNA="always"))
finalGenotypesIndel <- simplify2array(mclapply(finalIndel[whiteList], getGenotype, mc.cores=2, useNA="always"))
finalGenotypes <- aperm(abind::abind(subs=finalGenotypesSnv,indels=finalGenotypesIndel, along=5), c(1,5,2,3,4))
rm(finalGenotypesSnv,finalGenotypesIndel)

#' ### Probabilistic genotypes
finalGenotypesSnvP <- simplify2array(mclapply(finalSnv[whiteList], probGenotype, mc.cores=2))
finalGenotypesIndelP <- simplify2array(mclapply(finalIndel[whiteList], probGenotype, mc.cores=2))
finalGenotypesP <- aperm(abind::abind(subs=finalGenotypesSnvP,indels=finalGenotypesIndelP, along=4), c(1,4,2,3))
rm(finalGenotypesSnvP,finalGenotypesIndelP)

#' ### Probabilistic genotypes - tail prob (QC)
finalGenotypesSnvQ <- simplify2array(mclapply(finalSnv[whiteList], probGenotypeTail, mc.cores=2))
finalGenotypesIndelQ <- simplify2array(mclapply(finalIndel[whiteList], probGenotypeTail, mc.cores=2))
finalGenotypesQ <- aperm(abind::abind(subs=finalGenotypesSnvQ,indels=finalGenotypesIndelQ, along=3), c(1,3,2))
rm(finalGenotypesSnvQ,finalGenotypesIndelQ)

#' ### Save output
save.image(file=paste0(Sys.Date(),"-PCAWG-final.RData"))
#save(finalGenotypes, finalGenotypesP, finalGenotypesQ, file=paste0(Sys.Date(),"-FinalGenotypes.RData"))

#+ evalOn, echo=FALSE
opts_chunk$set(eval=TRUE)
load("2017-05-10-PCAWG-final.RData")

#' ### Duplicated samples
w <- names(finalSnv)
n <- names(which(table(sample2donor[w]) > 1)) # donors
s <- w[w %in% names(sample2donor[sample2donor %in% n])] # multisamples
d <- setdiff(sample2donor[w], sample2donor[levels(finalDrivers$sample)]) # donors w/o driver
u <- sample2donor[s[sample2donor[s] %in% intersect(d,n)]]
selectedSamples <- !w %in% setdiff(s[!s %in% finalDrivers$sample ], names(u)[!duplicated(u)])
uniqueSamples <- !duplicated(sample2donor[names(finalSnv)])

#' ### Overall distribution
#' #### Subs or indels
f <- function(x) unlist(sapply(seq_along(x), function(i) rep(i, x[i])))
d <- t(asum(finalGenotypesP[,"subs",,], 1))
o <- order(droplevels(donor2type[sample2donor[rownames(d)]]), -d[,1]/rowSums(d))
I <- t(apply(d/rowSums(d), 1, function(x) f(mg14:::roundProp(x * 100,p=100))))
d <- t(asum(finalGenotypesP[,"indels",,], 1))
J <- t(apply(d/rowSums(d), 1, function(x) if(!any(is.nan(x))) f(mg14:::roundProp(x * 100,p=100)) else rep(NA,100)))
s <- cumsum(table(droplevels(donor2type[sample2donor[rownames(d)]][o])))

#+ finalMutationsProb, fig.width=9, fig.height=8
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

#' #### Subs + indels
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


#' #### Barplot drivers
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
barplot(t(gu[2:51,]), col=col, las=2,  names.arg=rep("",50))
mg14::rotatedLabel(x=.Last.value,labels=rownames(gu)[2:51], cex=.5)
#dev.copy2pdf(file="finalDrivers.pdf", width=9, height=5, pointsize=8)

#' #### Barpot drivers - proportions
#+ finalDriversProp, fig.width=9, fig.height=3
par(fig=c(0,1,0,1),mar=c(3,4,1,1)+.1, mgp=c(3,.5,0))
w <- rowSums(pu) > 0
n <- 50
barplot(t(pu /rowSums(pu))[,w], width=c(rep(2,n+1), rep(0.2,sum(w)-n-1)), space=c(0,2,rep(0.1,n), rep(0,sum(w)-n-2)), col=col, border=NA, ylab="Fraction of mutations", names.arg=rep("",sum(w)))
mg14::rotatedLabel(x=.Last.value[1:(n+1)],labels=c("Genome-wide", rownames(pu)[2:(n+1)]), cex=.5)

#s <- 12/8
#dev.copy2pdf(file="finalDriversProp.pdf", width=9*s, height=3*s, pointsize=8*s)

#' #### Cumulative effects
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



#' ## WGD analyses
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

#pdf("WGD-timing-009.pdf", 12,6.5)
#j <- 1
#for(ID in names(finalBB)[isWgd]){
#	if(j%%100 == 0) print(j); j <- j+1
##	for(f in dir("../final/annotated_007/snv_mnv", pattern=paste0(ID,".+bb_granges.RData"), full.names=TRUE)) load(f)
##	t <- bbToTime(bb)
#	par(mfrow=c(2,1), mar=c(3,3,2,1), mgp=c(2,.5,0), bty="L", cex=1, las=2)
#	plotBB(finalBB[[ID]], ylim=c(0,8))
#	title(main=paste0(ID,", ", donor2type[sample2donor[ID]], ", ploidy=",round(finalPloidy[ID],2), ", hom=",round(finalHom[ID],2)), font.main=1, line=0)
#	par(mar=c(3,3,2,1))
#	plotTiming(finalBB[[ID]])
#	abline(h=fracGenomeWgdComp[ID,"time.wgd"], lty=3)
#	title(main=paste0("Timeable=", round(fracGenomeWgdComp[ID,2]/chrOffset["MT"]*100), "%, WGD=",round(fracGenomeWgdComp[ID,1]/fracGenomeWgdComp[ID,2]*100), "%, sd.WGD=",round(fracGenomeWgdComp[ID,'sd.wgd'],2), "%, avg.CI=",round(fracGenomeWgdComp[ID,'avg.ci'],2),", verdict=",wgdStar[ID]), font.main=1, line=0)
#}
#dev.off()

##' BB without VCF
#otherBB <- sapply(setdiff(sub("\\..+","",dir(bbPath)), names(finalBB)), function(ID) loadConsensusCNA(ID, purity=purityPloidy[ID, 'purity']))
#for(ID in names(otherBB))
#	otherBB[[ID]]$total_cn <- otherBB[[ID]]$major_cn+ otherBB[[ID]]$minor_cn
#
#otherPloidy <- sapply(otherBB, averagePloidy)
#otherHom <- sapply(otherBB, averageHom)
#otherWgd <- .classWgd(otherPloidy, otherHom)
#otherPoss <- !otherWgd & 2.5 - 1.5 * otherHom <= otherPloidy
#otherStat <- factor(otherPoss + 2*otherWgd - otherPoss*otherWgd, levels=0:2,labels=c("absent","possible","present"))


#' ## Coamplification and WGD
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



#fracGenomeWgdComp <- data.frame(fracGenomeWgdComp)
#attach(fracGenomeWgdComp)
#x <- mg14:::violinJitter()
#i <- !is.na(nt.coamp/nt.amp) & avg.ci < 0.5 & chr.all > 1
#
#par(mfrow=c(1,3), bty="n", mgp=c(2,.5,0), mar=c(3,3,1,1)+.1, cex=1)
#mg14:::violinJitterPlot((chr.wgd/chr.all)[i]*100, factor(WGD[i], labels=c("ND","WGD")), cex=1*sqrt(chr.all[i]/23), pch=16, ylim=c(0,100), ylab="Co-amplified chromosomes",  col.pty=rep("#00000088",2), plot.violins=FALSE)
#mg14:::violinJitterPlot((n.wgd/n.all)[i]*100, factor(WGD[i], labels=c("ND","WGD")), cex=2*sqrt(n.all[i]/1000), pch=16, ylim=c(0,100), ylab="Co-amplified segments",  col.pty=rep("#00000088",2), plot.violins=FALSE)
#mg14:::violinJitterPlot((nt.wgd/nt.total)[i]*100, factor(WGD[i], labels=c("ND","WGD")), cex=1*sqrt(nt.total[i]/3e9), pch=16, ylim=c(0,100), ylab="Co-amplified nucleotides",  col.pty=rep("#00000088",2), plot.violins=FALSE)

#' ## Timing examples
p <- function() {
	stackTime <- function(bb, t=seq(0,1,0.01)){
		u <- unique(bb)
		w <- as.numeric(width(u))
		f <- function(x) pmin(pmax(x,0.01),0.99)
		ut <- f((0.5*5+u$time * u$n.snv_mnv)/(5+u$n.snv_mnv))
		uu <- f(u$time.up)
		ul <- f(u$time.lo)
		diff(car::logit(f(t))) * rowSums(sapply(which(!is.na(ut)), function(i) w[i]*dnorm(car::logit(t[-1] - diff(t)/2), mean=car::logit(ut[i]), sd= (car::logit(uu[i]) - car::logit(ul[i]) + 0.05)/4)))#(t <= u$time.up[i] & t >= u$time.lo[i])))
		#rowSums(sapply(which(!is.na(ut)), function(i) w[i]*(t <= u$time.up[i] & t >= u$time.lo[i])))
	}
	layout(matrix(1:3, ncol=1), height=c(4,2,4))
	par(mar=c(0.5,3,0.5,0.5), mgp=c(2,0.25,0), bty="L", las=2, tcl=-0.25, cex=1)
	plotVcf(finalSnv[[w[1]]], finalBB[[w[1]]], finalClusters[[w[1]]], title=FALSE, legend=FALSE, col.grid='white',  xaxt=FALSE, cex=0.33)
	plotBB(finalBB[[w[1]]], ylim=c(0,5), legend=FALSE, type='bar', col.grid='white', col=c("lightgrey", "darkgrey"), xaxt=FALSE)
	par(mar=c(3,3,0.5,0.5))
	plotTiming(finalBB[[w[1]]], legend=FALSE, col.grid=NA)
	s <- stackTime(finalBB[[w[1]]])
	g <- colorRampPalette(RColorBrewer::brewer.pal(4,"Set1")[c(3,2,4)])(100)
	segments(x0=chrOffset["MT"] ,y0=seq(0,1,l=100),x1=chrOffset["MT"] + s/max(s) * 1e8, col=g, lend=3)
	print(w[1])
}

#+ timingExamples, fig.width=4, fig.height=4
w <- which(wgdStar=="likely" & !isWgd)
#pdf(paste0(names(w[1]), ".pdf"), 4,4, pointsize=8)
p()
#dev.off()

w <- which(wgdStar=="very likely" & isWgd)
#pdf(paste0(names(w[1]), ".pdf"), 4,4, pointsize=8)
p()
#dev.off()

w <- which(wgdStar=="unlikely" & !isWgd & fracGenomeWgdComp[,"nt.total"]/chrOffset["MT"] > 0.25 & fracGenomeWgdComp[,"avg.ci"] < 0.5)
#pdf(paste0(names(w[1]), ".pdf"), 4,4, pointsize=8)
p()
#dev.off()

#' GBM examples
#+ timingExamplesGbm, fig.width=4, fig.height=4
w <- which(fracGenomeWgdComp[,"time.wgd"]<0.1 & fracGenomeWgdComp[,"nt.total"]/chrOffset["MT"] > 0.1 &  !isWgd & donor2type[sample2donor[names(finalBB)]]=="CNS-GBM")
#pdf(paste0(names(w[1]), ".pdf"), 4,4, pointsize=8)
p()
#dev.off()



#' ## Signatures
sigTable <- simplify2array(mclapply(finalSnv, function(vcf) table(classifyMutations(vcf, reclassify="none"), tncToPyrimidine(vcf)), mc.cores=2))
sigTable <- aperm(sigTable, c(2,3,1))

#' ## Real-time WGD & subclones
#' ### Prelim
wgdTime <- function(vcf, bb, clusters, purity){
	w <- which(info(vcf)$MajCN==2 & info(vcf)$MinCN==2& sapply(info(vcf)$CNID, length)==1 & isDeamination(vcf))
	v <- vcf[w]
	if(nrow(v)<=100) return(NULL)
	seqnames(rowRanges(v)) <- factor(rep(1, nrow(v)), levels=seqlevels(v))
	b <- GRanges(1, IRanges(1,max(end(v))), copy_number=4, major_cn=2, minor_cn=2, clonal_frequency=purity)
	computeMutCn(v, b, clusters, purity, isWgd=TRUE, n.boot=10)
}

nClones <- sapply(finalClusters, nrow)

finalWgdParam <- mclapply(names(finalSnv)[isWgd], function(ID){
			wgdTime(finalSnv[[ID]], finalBB[[ID]], clusters=finalClusters[[ID]], purity=finalPurity[ID])
		},  mc.cores=4)

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

correctAccelRand <- function(pi, ta=seq(0.66,1,0.01), a=seq(1,10,1)){
	x <- sapply(ta, function(taa) sapply(a, function(aa) correctAccel(pi, taa, aa)), simplify='array')
}

accel <- c(1,2.5,5,7.5,10)
names(accel) <- paste0(accel, "x")

#' ### WGD
foo <- apply(finalWgdPi[,,], 2:3,  function(x) correctAccelRand(x, a=accel))
finalWgdPiAdj <- array(foo, dim=c(2, length(accel), length(eval(formals(correctAccelRand)$ta)), dim(foo)[-1]))
dimnames(finalWgdPiAdj)[[1]] <- c('t.WGD','t.subclonal')
dimnames(finalWgdPiAdj)[[2]] <- paste0(accel, "x")
dimnames(finalWgdPiAdj)[[4]] <- dimnames(finalWgdPi)[[2]]
dimnames(finalWgdPiAdj)[[5]] <- names(finalBB)[isWgd][!sapply(finalWgdParam, is.null)]

#finalWgdPiAdj <- sapply(1:ncol(finalWgdPi), function(i) correctAccelRand(finalWgdPi[,i], a=accel), simplify='array')
#dimnames(finalWgdPiAdj)[[2]] <- paste0(accel, "x")

age <- clinicalData$donor_age_at_diagnosis
names(age) <- clinicalData$icgc_donor_id

n <- dimnames(finalWgdPiAdj)[[5]]
finalWgdTime <- finalWgdPiAdj[,,,,n] * rep(age[sample2donor[n]], each=2)

typeNa <- gsub("\t","",strsplit("Bone-Cart
				Breast-LobularCA
				Breast-DCIS
				Lymph-NOS
				Myeloid-MDS
				Cervix-AdenoCA", "\n")[[1]])

d <- droplevels(donor2type[sample2donor[n]])
s <- setdiff(levels(d), c(typeNa, names(which(table(d)<5))))
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

tissueBorder <- names(tissueColors) %in% c("Lung-SCC","Lung-AdenoCA")
names(tissueBorder) <- names(tissueColors)

#+ realTimeWgd, fig.height=3, fig.width=4
par( mar=c(7,3,1,1), mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
u <- names(finalSnv)[uniqueSamples]
qWgd <- sapply(timeWgd, function(x) apply(x[rownames(x) %in% u,"hat",], 2, quantile, c(0.25,0.5,0.75), na.rm=TRUE), simplify='array')
m <- qWgd[3,"5x",]#t[1,3,]
a <- "7.5x"
m <- qWgd[3,a,]#t[1,3,]
o <- order(m, na.last=NA)
plot(NA,NA, xlim=c(0.5,length(m[o])), ylim=c(0,max(2*qWgd[3,a,o], na.rm=TRUE)+.5), ylab="Years before diagnosis", xlab="", xaxt="n", yaxs="i")
x <- seq_along(m[o])
#rect(x-.45,qWgd[1,"10x",o] ,x+.45, qWgd[3,"1x",o],col=NA, border='grey') #mg14::colTrans(tissueColors[names(m)[o]], 3)
#segments(x-.45,t(q[3,1:3,o]) ,x+.45, t(q[3,1:3,o]),col='grey') #mg14::colTrans(tissueColors[names(m)[o]], 3)
#rect(x-.45,qWgd[1,a,o] ,x+.45, qWgd[3,a,o],col=mg14::colTrans(tissueColors[names(m)[o]],2))
mg14::rotatedLabel(x, labels=names(sort(m)))
y <- sapply(timeWgd, `[`, (quote(f(,)))[[2]], 1:3,a)
for(i in seq_along(o)){
	f <- function(x) x/max(abs(x))
	j <- f(mg14::violinJitter(na.omit(y[[o[i]]][,"hat"]))$y)/4 + i
	segments(j, na.omit(y[[o[i]]][,"up"]), j, na.omit(y[[o[i]]][,"lo"]), col='#DDDDDD')
	points(j, na.omit(y[[o[i]]][,"hat"]), pch=21, col=if(tissueColors[names(m)[o[i]]]=="#000000") "white" else "black", bg=tissueColors[names(m)[o[i]]], cex=1)
}
par(xpd=TRUE)
#s <- 12/8
#dev.copy2pdf(file="realTimeWgd.pdf", width=4*s, height=3.5*s, pointsize=8*s)

sapply(timeWgd, nrow)

#+ realTimeWgdAccel, fig.height=2, fig.width=2
par(mar=c(3,3,1,1), mgp=c(2,0.5,0), tcl=-0.5, bty="L")
plot(accel, qWgd["50%",,1], type='l', lty=0, ylim=c(0,30), xlab= "5meC rate acceleration", ylab="Median occurrence [years]", yaxs="i", xaxt="n")
axis(side=1, at=accel)
for(j in 1:dim(qWgd)[3]) lines(accel, qWgd["50%",,j], type='l', col=tissueColors[dimnames(qWgd)[[3]][j]])
#s <- 12/8
#dev.copy2pdf(file="realTimeWgdAccel.pdf", width=2*s, height=2*s, pointsize=8*s)


#' ### Subclones
effGenome <- unlist(mclapply(finalSnv, function(vcf) 2/avgWeights(vcf[na.omit(info(vcf)$CLS!="subclonal")], type="deam"), mc.cores=4))
names(effGenome) <- names(finalSnv)

subcloneDeam <- t(simplify2array(mclapply(finalSnv, function(vcf) {p <- info(vcf)$pSub[isDeamination(vcf)]; c(sum(p, na.rm=TRUE), sum(1-p, na.rm=TRUE))})))

d <- droplevels(donor2type[sample2donor[names(finalSnv)]])
s <- setdiff(levels(d), c(typeNa, names(which(table(d)<5))))

set.seed(42)
timeSubclones <- sapply(s, function(l) {
			i <- d==l
			tt0 <- subcloneDeam[i,]/cbind(finalPloidy[i], effGenome[i]) / cbind(nClones[i]-1, 1)
			resB <- sapply(1:1000, function(foo){
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
			r <- which(rowSums(subcloneDeam[i,]) < 50)
			arr[r,,] <- NA
			return(arr)
		})

t <- sapply(timeSubclones, colMeans, na.rm=TRUE)

m <- t[3,]

o <- order(m, na.last=NA)

timeSubclones0 <- sapply(s, function(l) {
			i <- d==l
 			tt <- subcloneDeam[i,]/cbind(finalPloidy[i], effGenome[i]) / cbind(nClones[i]-1, 1)
 			res <- sapply(accel, function(a)  tt[,1]/a/rowSums(tt/rep(c(a,1), each=nrow(tt)))) * age[sample2donor[names(finalSnv)[i]]]
 			colnames(res) <- paste0(accel, "x")
 			#res[res==0] <- NA
 			res})

#' Plot
#+ realTimeSubclone, fig.width=5, fig.height=4.667
u <- names(finalSnv)[uniqueSamples]
par( mar=c(7,3,1,1), mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
#qSubclone <- sapply(timeSubclones, function(x) apply(x[,], 2, quantile, c(0.25,0.5,0.75), na.rm=TRUE), simplify='array')
qSubclone <- sapply(timeSubclones, function(x) apply(x[rownames(x)%in%u,"hat",], 2, quantile, c(0.25,0.5,0.75), na.rm=TRUE), simplify='array')
a <- "7.5x"
m <- qSubclone[2,a,]#t[1,3,]
o <- order(m, na.last=NA)
plot(NA,NA, xlim=c(0.5,length(m[o])), ylab="Years before diagnosis", xlab="", xaxt="n", yaxs="i", ylim=c(0,30))
x <- seq_along(m[o])
#rect(x-.45,qWgd[1,"10x",o] ,x+.45, qWgd[3,"1x",o],col=NA, border='grey') #mg14::colTrans(tissueColors[names(m)[o]], 3)
#segments(x-.45,t(q[3,1:3,o]) ,x+.45, t(q[3,1:3,o]),col='grey') #mg14::colTrans(tissueColors[names(m)[o]], 3)
#rect(x-.45,qWgd[1,a,o] ,x+.45, qWgd[3,a,o],col=mg14::colTrans(tissueColors[names(m)[o]],2))
mg14::rotatedLabel(x, labels=names(sort(m)))
y <- sapply(timeSubclones, `[`, (quote(f(,)))[[2]], 1:3,a)
for(i in seq_along(o)){
	f <- function(x) x/max(abs(x))
	j <- f(mg14::violinJitter(na.omit(y[[o[i]]][,"hat"]))$y)/4 + i
	segments(j, na.omit(y[[o[i]]][,"97.5%"]), j, na.omit(y[[o[i]]][,"2.5%"]), col='#DDDDDD')
	points(j, na.omit(y[[o[i]]][,"hat"]), pch=21, col=if(tissueColors[names(m)[o[i]]]=="#000000") "white" else "black", bg=tissueColors[names(m)[o[i]]], cex=1)
}

#par(xpd=TRUE)
#s <- 12/8
#dev.copy2pdf(file="realTimeSubclone.pdf", width=5*s, height=3.5*3/4*s, pointsize=8*s)

sapply(timeSubclones, nrow)

#+ realTimeSubcloneWgdScatter
par( mar=c(4,3,1,1), mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
plot(qSubclone["50%",a,dimnames(qWgd)[[3]]], qWgd["50%",a,], col=tissueColors[dimnames(qWgd)[[3]]], pch=16, cex=2, xlab="Median time subclones", ylab="Median time WGD", xlim=c(0,5), ylim=c(0,10))

#+ realTimeSubcloneWgd, fig.width=2.5, fig.height=3.5
par( mar=c(3,3,3,10), mgp=c(2,.5,0), tcl=-0.25,cex=1, bty="n", xpd=FALSE, las=1)
w <- "50%"
plot(c(rep(1, dim(qSubclone)[3]), rep(2, each=dim(qWgd)[3])), c(qSubclone[w,a,],qWgd[w,a,]), bg=tissueColors[c(dimnames(qSubclone)[[3]], dimnames(qWgd)[[3]])], pch=21, cex=2, xaxt="n", ylab="Years before diagnosis", xlab="", xlim=c(0.5,2.5))
segments(rep(1, each=dim(qWgd)[3]), qSubclone[w,a,dimnames(qWgd)[[3]]], rep(2, each=dim(qWgd)[3]), qWgd[w,a,],col=tissueColors[dimnames(qWgd)[[3]]])
o <- order(qWgd[w,a,], na.last=NA)
y0 <- qWgd[w,a,o]
y1 <- mg14:::mindist(qWgd[w,a,o], diff(par('usr')[3:4])/30)
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

#' ## Session
#' ### Objects
l <- ls()
data.frame(variable=l, Reduce("rbind",lapply(l, function(x) data.frame(class=class(get(x)), size=format(object.size(get(x)), units="auto")))))
#' ### Packages
sessionInfo()
devtools::session_info()