#' ---
#' title: "Supplementary code: The evolutionary history of 2,658 cancers" 
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 4
#'     toc_float: true
#'     number_sections: true
#'     auto_identifiers: true
#'     table_captions: true
#' author: Moritz Gerstung and Santiago Gonzalez, on behalf of the PCAWG-11 Evolution and Heterogeneity Working Group
#' ---

#+ Preliminaries, echo=FALSE
options(width=120)
pdf.options(pointsize=8)
.par <- function() par(mar=c(3,3,1,1), bty="L", mgp=c(2,.5,0), tcl=-0.25, las=1)
require(knitr)
knit_hooks$set(smallMar = function(before, options, envir) {
			if (before) .par()
		})
opts_chunk$set(dev=c('my_png','pdf'), fig.ext=c('png','pdf'), fig.width=3, fig.height=3, smallMar=TRUE)
my_png <-  function(file, width, height, pointsize=12, ...) {
	png(file, width = 1.5*width, height = 1.5*height, units="in", res=72*1.5, pointsize=pointsize, ...)
}
options(mc.cores=as.numeric(Sys.getenv("LSB_MAX_NUM_PROCESSORS")))


#' # Prelim
#' ## Preprocessing
#' Prior to this script, `MutationTimeR` was run on 2x 2,778 VCF files for subs and indels each. This was done using the script `VCF-annotate.R`, attached at
#' the end of this vignette. 
#' 
#' Input files:
#' 
#' * VCF files: Folder `../final/final_consensus_12oct_passonly`
#' * Allele-specific copy number: Folder `../final/consensus.20170119.somatic.cna.annotated`
#' * Subclone sizes and cell frequencies: `../final/structure_weme_released_consensus_merged.txt`, `../final/wcc_consensus_values_9_12.tsv`
#' * Purity: `../final/consensus.20170218.purity.ploidy.txt` 
#' 
#' Output:
#+ files, ignore.stderr = FALSE
cat(paste0(system("for d in `find ../final/annotated_014 -maxdepth 1 -type d`
do echo $d; ls $d | head -6;
done", ignore.stderr=TRUE, intern=TRUE), collapse="\n"))
#' Annotated VCF files and copy number segments with timing information will be loaded as shown below. This requires about 100G memory.
#' 
#' ## Libraries
#' Load convenience R functions (code shown at the end of this vignette).
source("PCAWG-functions.R")

#' Prepare spreadsheets for figure data
library(xlsx)
Figure1 <- createWorkbook()
Figure2 <- createWorkbook()
Figure4 <- createWorkbook()
Figure5 <- createWorkbook()
ExtendedDataFigure3 <- createWorkbook()
ExtendedDataFigure6 <- createWorkbook()
ExtendedDataFigure8 <- createWorkbook()
ExtendedDataFigure9 <- createWorkbook()


#+ evalOff, echo=FALSE
dumpfile <- "2018-07-18-PCAWG-final.RData"
if(file.exists(dumpfile)){
	opts_chunk$set(eval=FALSE) # ie skip following steps
	load(dumpfile)
	source("PCAWG-functions.R")
}

#' # Load processed data from MutationTimeR
#' 
#' ## Whitelist
#' First 2,703 whitelisted samples
#' 
#' ### SNV and MNV
#' Load annotated VCF files for SNVs and MNVs into a list with `VariantAnnotation::VCF()` objects:
#+ loadSNV
p <- "../final/annotated_014/snv_mnv"
d <- dir(p, pattern="*.vcf.RData", full.names=TRUE)
finalSnv <- unlist(mclapply(split(d, seq_along(d) %/% 100), lapply, function(f) { # read in batches of 100
			e <- new.env()
			load(f, envir=e)
			e$vcf
		}, mc.preschedule=FALSE), recursive=FALSE)
names(finalSnv) <- sub(".conse.+","",dir(p, pattern="*.vcf.RData", full.names=FALSE))

#' ### Copy number profiles
#' Load copynumber profiles as a list of `GRanges()`:
#+ loadBB
p <- "../final/annotated_014/cn"
finalBB <- list()
for( f in dir(p, pattern="*.bb_granges.RData", full.names=TRUE)){
	load(f)
	colnames(mcols(bb)) <- sub("star.1","time.star",colnames(mcols(bb)) ) # Fix naming problem
	finalBB[[f]] <- bb
}
names(finalBB) <- sub(".conse.+","",dir(p, pattern="*.bb_granges.RData", full.names=FALSE))

#' ### Indels
#' Load annotated VCF files for indels into a list with `VariantAnnotation::VCF()` objects:
#+ loadIndel
p <- "../final/annotated_014/indel"
d <- dir(p, pattern="*.vcf.RData", full.names=TRUE)
finalIndel <- unlist(mclapply(split(d, seq_along(d) %/% 100), lapply, function(f) { # read in batches of 100
					e <- new.env()
					load(f, envir=e)
					e$vcfIndel
				}, mc.preschedule=FALSE), recursive=FALSE)
names(finalIndel) <- sub(".conse.+","",dir(p, pattern="*.vcf.RData", full.names=FALSE))

#' ### Clusters and purity
#' These were input to MutationTimeR. Loaded as a `list()` (clusters) and `data.frame` (purity/ploidy).
#+ loadClusters
finalClusters <- list()
finalPurity <- numeric()
p <- "../final/annotated_014/clusters"
for( f in dir(p, pattern="*.RData", full.names=TRUE)){
	load(f)
	finalClusters[[f]] <- clusters
	finalPurity[f] <- purity
}
names(finalClusters) <- names(finalPurity) <- sub(".conse.+","",dir(p, pattern="*.RData", full.names=FALSE))


#' ## Update drivers
#' `finalDrivers` is a `VariantAnnotation::VRanges()` object with driver point mutations (substitutions and indels), as provided by the PCAWG drivers 
#' working group. Here, we match these data with annotation from MutationTimeR. To this end, match samples and positions and copy annotation. 
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
	if(!is.na(j)){
		mcols(finalDriversAnnotated)[i,colnames(d)] <- info(v)[j, colnames(d)]
		refDepth(finalDriversAnnotated)[i] <- info(v)[j,"t_ref_count"]
		altDepth(finalDriversAnnotated)[i] <- info(v)[j,"t_alt_count"]
	}
	else
		mcols(finalDriversAnnotated)[i,colnames(d)] <- NA
}

#' ## Graylisted data
#' 75 graylisted samples were processed separately using MutationTimeR. We load these just like the whitelisted samples and concatenate the output.
#' 
#' ### SNV and MNV
#+ loadSnvGray
p <- "../final/annotated_014/graylist/snv_mnv"
finalSnvGray <- mclapply(dir(p, pattern="*.vcf.RData", full.names=TRUE), function(f) {
			e <- new.env()
			load(f, envir=e)
			e$vcf
		})
names(finalSnvGray) <- sub(".conse.+","",dir(p, pattern="*.vcf.RData", full.names=FALSE))
finalSnv[names(finalSnvGray)] <- finalSnvGray

#' ### Copy number profiles
#+ loadBBGray
p <- "../final/annotated_014/graylist/cn"
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
p <- "../final/annotated_014/graylist/indel"
finalIndelGray <- mclapply(dir(p, pattern="*.vcf.RData", full.names=TRUE), function(f) {
			e <- new.env()
			load(f, envir=e)
			e$vcfIndel
		})
names(finalIndelGray) <- sub(".conse.+","",dir(p, pattern="*.vcf.RData", full.names=FALSE))
finalIndel[names(finalIndelGray)] <- finalIndelGray


#' ### Clusters and purity
#+ loadClustersGray
finalClustersGray <- list()
finalPurityGray <- numeric()
p <- "../final/annotated_014/graylist/clusters"
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

#' ## Structural variants
#' Lastly, we load structural variant data. This is only available for a subset. We pad missing samples with NA
finalSv <- mclapply(dir("../final/pcawg_consensus_1.6.161116.somatic_svs", pattern='*.vcf.gz$', full.names=TRUE), function(x) {
			t <- try(readVcf(x))
			return(t)
		})
names(finalSv) <- sub("../final/pcawg_consensus_1.6.161116.somatic_svs/","", sub(".pcawg_consensus_1.6.161116.somatic.sv.vcf.gz","",dir("../final/pcawg_consensus_1.6.161116.somatic_svs", pattern='*.vcf.gz$', full.names=TRUE)))
finalSv <- finalSv[names(finalSnv)]


#' # QC
#' Here we calculate basic QC for MutationTimeR output. `info(vcf)$pMutCNTail` contains the beta-binomial tail probability `pbetabinom` for every variant
#' for its most likely clonal state. These should be roughly uniform, hence 1% outside [0.005,0.995] and about 5% outside [0.025,0.975]. Larger deviations
#' indicate that either purity, copy number or clusters were wrong in a given sample.
#+ QC
q1 <- sapply(finalSnv, function(vcf) mean(abs(0.5- info(vcf)$pMutCNTail) > 0.495 , na.rm=TRUE))
q5 <- sapply(finalSnv, function(vcf) mean(abs(0.5- info(vcf)$pMutCNTail) > 0.475 , na.rm=TRUE))

#' Most samples agree, on average:
#+ QCboxplot, eval=TRUE
par(mfrow=c(1,1))
boxplot(1-q5 ~ donor2type[sample2donor[names(finalSnv)]], las=2, ylab="Fraction of data inside theoretical 95% CI")
abline(h=0.95, lty=3)

#' QQ-plots show few outliers per sample:
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
#' Here we tabulate all driver point mutations and their clonal allele status.
#' 
#' ## MAP genotypes
#' These are the maximum a posteriori (MAP) assignments, used for the early/late/clonal/subclonal annotation output of MutationTimeR.
#+ finalGenotypes
finalGenotypesSnv <- simplify2array(mclapply(finalSnv[whiteList], getGenotype, useNA="always"))
finalGenotypesIndel <- simplify2array(mclapply(finalIndel[whiteList], getGenotype, useNA="always"))
finalGenotypes <- aperm(abind::abind(subs=finalGenotypesSnv,indels=finalGenotypesIndel, along=5), c(1,5,2,3,4))
rm(finalGenotypesSnv,finalGenotypesIndel)

#' ## Probabilistic genotypes
#' Alternatively to the hard MAP assignments one can calculate the actual probailities per state and assess their overall distribution. This could avoid
#' some biases of the MAP estimates, which are typically biased towards the more populous clonal states.
#+ finalGenotypesP
finalGenotypesSnvP <- simplify2array(mclapply(finalSnv[whiteList], probGenotype))
finalGenotypesIndelP <- simplify2array(mclapply(finalIndel[whiteList], probGenotype))
finalGenotypesP <- aperm(abind::abind(subs=finalGenotypesSnvP,indels=finalGenotypesIndelP, along=4), c(1,4,2,3))
rm(finalGenotypesSnvP,finalGenotypesIndelP)

#' ## Probabilistic genotypes - tail prob (QC)
#' As QC, we aggregate also the tail probabilities to find out whether there were annotation issues for some driver genes
#+ finalGenotypesQ
finalGenotypesSnvQ <- simplify2array(mclapply(finalSnv[whiteList], probGenotypeTail))
finalGenotypesIndelQ <- simplify2array(mclapply(finalIndel[whiteList], probGenotypeTail))
finalGenotypesQ <- aperm(abind::abind(subs=finalGenotypesSnvQ,indels=finalGenotypesIndelQ, along=3), c(1,3,2))
rm(finalGenotypesSnvQ,finalGenotypesIndelQ)

#' # Checkpoint
#' As the previous steps take rather long, introduce a checkpoint here and serialise the data, so it can be loaded later on.
#+ saveOut
save.image(file=dumpfile, compress=FALSE)
save(finalGenotypes, finalGenotypesP, finalGenotypesQ, file=paste0(Sys.Date(),"-FinalGenotypes.RData"))

#+ evalOn, eval=TRUE, echo=FALSE
opts_chunk$set(eval=TRUE)

#' # Timing of point mutations
#' Here we assess the timing of point mutations, including driver gene mutations.
#' 
#' ## Duplicated donors
#' First, find samples with the same donor to avoid double dipping for some computations.
w <- names(finalSnv)
n <- names(which(table(sample2donor[w]) > 1)) # donors
s <- w[w %in% names(sample2donor[sample2donor %in% n])] # multisamples
d <- setdiff(sample2donor[w], sample2donor[levels(finalDrivers$sample)]) # donors w/o driver
u <- sample2donor[s[sample2donor[s] %in% intersect(d,n)]]
selectedSamples <- !w %in% setdiff(s[!s %in% finalDrivers$sample ], names(u)[!duplicated(u)])
uniqueSamples <- !duplicated(sample2donor[names(finalSnv)])

#' ## Overall distribution
#' 
#' ### Subs or indels
#' First assess the overall distributions of early/late/clonal/subclonal variants per sample, separately for subs and indels.
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

#' The proportions of subs and indels are remarkably similar.
#' 
#' ### Figure 2a
#' For the manuscript only plot subs, together with the color code for each tumour type.
f <- function(x) unlist(sapply(seq_along(x), function(i) rep(i, x[i])))
d <- t(sapply(names(finalSnv)[whiteList &  selectedSamples], function(n) table(info(finalSnv[[n]])$CLS, useNA='a') + table(info(finalIndel[[n]])$CLS, useNA='a')))
o <- order(droplevels(donor2type[sample2donor[rownames(d)]]), -d[,1]/rowSums(d))
I <- t(apply(d/rowSums(d), 1, function(x) f(mg14:::roundProp(x * 100,p=100))))
s <- cumsum(table(droplevels(donor2type[sample2donor[rownames(d)]][o])))

col <- RColorBrewer::brewer.pal(9, "Set1")[c(3,4,2,1,9)] ## Colors for early-subclonal

#+ finalMutationPropAll, fig.width=7.25, fig.height=1.5
par(fig=c(0,1,0,1),mar=rep(0,4), bty="n", tcl=-0.25)
layout(matrix(1:5,nrow=5), height=c(1,4,1.2,0.8,4))
par(cex=1)
plot(NA,NA, xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0, nrow(I)), ylim=c(0,1), xaxs="i", yaxs='i')
legend("bottom", fill=col, legend=paste(dimnames(finalGenotypes)[[3]]), bty="n", horiz=TRUE, title="Mutation timing", cex=1)
par(mar=c(0,4,0,0)+.1, mgp=c(2,.5,0), bty="n", las=2, xpd=FALSE, cex=1)
image(z=I[o,], x=1:nrow(I), useRaster=TRUE, col=col, xaxt="n", ylab="Point mutations")
abline(v=s, col='white')
u <- par("usr")
par(bty="o")
barplot(t(t(table(droplevels(donor2type[sample2donor[rownames(d)]][o])))), horiz=TRUE, col=tissueColors[names(s)], border=NA, xlim=u[1:2], yaxs='i', xaxt="n")
par(bty="n", xpd=NA)
plot(NA,NA, xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0, nrow(I)), ylim=c(0,1), xaxs="i", yaxs='i')
d0 <- s - diff(c(0,s))/2
d1 <- mg14:::mindist(d0,60)
segments(d0, 0.66, d0, 1)
segments(d0, 0.66, d1, 0.33)
segments(d1, 0.33, d1, 0)
mg14::rotatedLabel(x0 = d1, names(s), y0=0)

#' Also dump the xlsx output.
Figure2a <- createSheet(Figure2, "Figure2a")
addDataFrame(data.frame(d,cancer=droplevels(donor2type[sample2donor[rownames(d)]]))[o,c(6,1:5)], Figure2a)

#t <- 12/8
#dev.copy2pdf(file="finalMutationPropAll.pdf", width=9*t, height=1.8*t, pointsize=8*t)


#' ### Barplot drivers
#' Not turn to driver point mutations. Extract from the large arrays calculated above.
p <- asum(finalGenotypesP[,,,selectedSamples[whiteList]], c(2,4))
g <- asum(finalGenotypes[,,,,selectedSamples[whiteList]], c(2,4:5))
g <- g[order(rowSums(g), decreasing=TRUE),]
colnames(g) <- paste(colnames(g))
rownames(g) <- paste(rownames(g)) 
rownames(p) <- paste(rownames(p))
p <- p[rownames(g),]
w <- rowSums(g) > 0
w[1] <- FALSE

#' Collapse into unique genes
l <- sapply(strsplit(rownames(p),"::"), function(x) paste(x[2:3], collapse=":"))
pu <- t(sapply(unique(l), function(u) asum(p[l==u,,drop=FALSE], 1)))
gu <- t(sapply(unique(l), function(u) asum(g[l==u,,drop=FALSE], 1)))
gu <- gu[order(rowSums(gu), decreasing=TRUE),]
pu <- pu[rownames(gu),]
wu <- rowSums(gu) > 0
wu[1] <- FALSE

#' Now plot as a barplot. All drivers, with the top 50 as an inset.
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

#' Also generate xlsx output
Figure2a2 <- createSheet(Figure2, "Figure2a2")
addDataFrame(gu[wu,], Figure2a2)

gt <- asum(finalGenotypes[,,,,selectedSamples[whiteList]], c(2:4))
t <- droplevels(donor2type[sample2donor[colnames(gt)]])
rownames(gt) <- paste(rownames(gt))
gt <- gt[rownames(g),] %*% model.matrix(~t-1)
colnames(gt) <- levels(t)
gtu <- t(sapply(unique(l), function(u) asum(gt[l==u,,drop=FALSE], 1)))
gtu <- gtu[rownames(gu),]

Figure2a3 <- createSheet(Figure2, "Figure2a3")
addDataFrame(gtu[wu,], Figure2a3)


#' Proportions may be more revealing than absolute numbers:
#+ finalDriversProp, fig.width=9, fig.height=2
par(fig=c(0,1,0,1),mar=c(3,4,1,1)+.1, mgp=c(3,.5,0))
w <- rowSums(pu) > 0
n <- 50
b <- barplot(t(pu /rowSums(pu))[,w], width=c(rep(2,n+1), rep(0.2,sum(w)-n-1)), space=c(0,2,rep(0.1,n), rep(0,sum(w)-n-2)), col=col, border=NA, ylab="Fraction of mutations", names.arg=rep("",sum(w)))
mg14::rotatedLabel(x=b[1:(n+1)],labels=c("Genome-wide", rownames(pu)[2:(n+1)]), cex=.5)

#s <- 12/8
#dev.copy2pdf(file="finalDriversProp.pdf", width=9*s, height=3*s, pointsize=8*s)

#' ### Cumulative distribution of driver mutations
#' Sort driver mutations by prevalence in each timing class. Clearly the MAP estimate is biased downwards for driver genes not observed in a final cohort size. 
#' Hence to this without and with a pseudocount of one (which will overestimate). Report the average of the two as the point estimate.
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

#' ### Figure 2d
#' Now plot the value where the number of unique genes contributes 50% of all driver mutations.
#+ finalGenes50, fig.width=3, fig.height=4
par(mar=c(4,3,2.5,1), mgp=c(2,.5,0), bty="L")
d50 <- apply((r[,,1]+r[,,2])/2 < 0.5, 2, which.min)[c(1,3,2,4)]
b <- barplot(d50,las=2, col=col[c(1,3,2,4)], border=NA, ylab="Genes contributing 50% of driver mutations")
lo <- apply(r[,,1] < 0.5, 2, which.min)[c(1,3,2,4)]
hi <- apply(r[,,2] < 0.5, 2, which.min)[c(1,3,2,4)]
segments(b,lo,b,hi)
mg14::rotatedLabel(x=b,labels=c("clonal [early]", "clonal [late]", "clonal [other]", "subclonal")[c(1,3,2,4)])
#dev.copy2pdf(file="finalGenes50.pdf", width=3,height=4)

#' xlxs output
#+ Figure2d
Figure2d <- createSheet(Figure2, "Figure2d")
addDataFrame(data.frame(d50, lo, hi, row.names=c("clonal [early]", "clonal [late]", "clonal [other]", "subclonal")[c(1,3,2,4)]), Figure2d)

#' # Whole-genome duplications
#' Now turn to the timing of copy number alterations. First assess the whole genome duplications (WGD).
#' 
#' ## Prelim
#' Final ploidy, weighted if subclonal CN
finalPloidy <- sapply(finalBB, averagePloidy)
names(finalPloidy) <- names(finalBB)

#' Final homozygousity, weighted if subclonal CN
finalHom <- sapply(finalBB, averageHom)
names(finalHom) <- names(finalBB)

#' ## WGD classification
#' WGD is primarily classified by the copy number data. Additionally we assess the concordance of timing. 
#' 
#' ### Based on ploidy and homozygousity
#' Using the quantities above, we classify WGD samples. This is the official WG output.
isWgd <- .classWgd(finalPloidy, finalHom)
table(isWgd)
#+ wdgHomPloidy, fig.widht=4, fig.height=4
plot(finalHom, finalPloidy, col=.classWgd( finalPloidy, finalHom)+1, xlim=c(0,1))

#' ### Based on timing
#' Further it's worthwile assessing whether the different copy number segments show the same timing, which would be expected for WGD.
#' We do this for both WGD and near diploid (ND) samples.
fracGenomeWgdComp <- t(sapply(finalBB, function(bb) {
					fgw <- try(fractionGenomeWgdCompatible(bb)); 
					if(class(fgw)!='try-error') fgw
					else rep(NA,10)}))
rownames(fracGenomeWgdComp) <- names(finalBB)

#' The temporal concordance defines out rating. Some samples are uninformative (especially ND) either due to too little point mutations (large timing 
#' CIs), or too few gained segments.
wgdStar <- factor(rep(1,nrow(fracGenomeWgdComp)), levels=0:3, labels=c("unlikely","uninformative","likely","very likely"))
wgdStar[fracGenomeWgdComp[,"avg.ci"]<=0.75 & fracGenomeWgdComp[,"nt.total"]/chrOffset["MT"] >= 0.33 ] <- "likely"
wgdStar[fracGenomeWgdComp[,"nt.wgd"]/fracGenomeWgdComp[,"nt.total"] < 0.66] <- "unlikely"
wgdStar[wgdStar=="likely" & fracGenomeWgdComp[,"nt.wgd"]/fracGenomeWgdComp[,"nt.total"] > 0.8 & fracGenomeWgdComp[,"sd.wgd"] < 0.1 &  fracGenomeWgdComp[,"nt.total"]/chrOffset["MT"] > 0.5] <- "very likely"
names(wgdStar) <-  names(finalBB)
prop.table(table(wgdStar[!isWgd]))

wgdPoss <- !isWgd & 2.5 - 1.5 * finalHom <= finalPloidy

wgdStat <- factor(wgdPoss + 2*isWgd - wgdPoss*isWgd, labels=c("absent","possible","present"))
table(wgdStat, wgdStar)

#' Overall, WGD samples display a high level of temporal concordance, as hoped.
#' 
#' # Temporal distribution of chromosomal gains
#' In this section we assess the temporal distribution of large-scale copy number gains.
#' 
#' ## Functions
#' This one aggregates individual segments by chromosome
aggregatePerChromosome <- function(bb, isWgd=FALSE){
	.aggregateSegments <- function(m){
		#m <- mcols(bb)
		t <- weighted.mean(m$time, m$n.snv_mnv, na.rm=TRUE)
		n <- sum(m$n.snv_mnv[!is.na(m$time)], na.rm=TRUE)
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

#' ## Aggregated timing per chromosome
#' This provides an average timing estimate per chromosome, per sample. 
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
#' Plot the timing histograms (deciles) per chromosome per cancer type, similar to Figure 1c in the final publication.
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
			abline(v=0, col='lightgrey', lty=3)
		if(i==2){
			abline(h=0.05*8/n, col='lightgrey', lty=1)
			axis(side=2, at=c(0,0.05*8/n), labels=c("",""), tcl=-.1)
		}
	}
#dev.copy2pdf(file="histTiming.pdf",width=6, height=6, pointsize=8)


vv <- v[dimnames(h)[[2]]]
vv <- vv/sum(vv)

hh <- matrix(matrix(aperm(h, c(1,3,2)), ncol=length(vv)) %*% vv, nrow=nrow(h))
rownames(hh) <- rownames(h)

#' ### Pan-Can histograms
#' It can also be instructive to study the timing histograms pan-cancer:
#+ histTimingPanCan, fig.height=2, fig.width=2
par(mar=c(3,3,1,1), mgp=c(2,.5,0), tcl=-0.5, bty="L", xpd=NA)
barplot(hh["WGD",], space=0, col=rev(col), xlab="Time [mutations]", ylab="Relative frequency", width=0.1, ylim=c(0,.065), yaxs='r', border=NA)
axis(side=1)
barplot(hh["All",], space=0, col=rev(col), xlab="Time [mutations]", ylab="Relative frequency", width=0.1, ylim=c(0,.065), yaxs='r', border=NA)
axis(side=1)
#dev.copy2pdf(file="histTimingPanCan.pdf",width=2, height=2, pointsize=8)

#' ## Glioblastoma
#' GBMs display very early timing on chromosomes 7, 20 and 21. Plot a few.
#+ timingExamplesGbm, fig.width=4, fig.height=4
w <- which(fracGenomeWgdComp[,"time.wgd"]<0.1 & fracGenomeWgdComp[,"nt.total"]/chrOffset["MT"] > 0.1 &  !isWgd & donor2type[sample2donor[names(finalBB)]]=="CNS-GBM")
plotSample(w[1])
plotSample(w[2])
plotSample(w[3])

#' ### Extended Data Figure 3a
#' Explore the early timing of chromosome 7 a bit more deeply.
#' Find all GBMs with +7
w <- which(donor2type[sample2donor[names(finalSnv)]]=="CNS-GBM") 
w <- w[sapply(finalBB[w], function(bb) sum(width(bb)[as.logical(seqnames(bb)==7) & bb$total_cn >= 3], na.rm=TRUE)/width(refLengths[7])>0.8)]

#' Plot the first 10
#+ GBM_tri7, fig.width=3.2, fig.height=3.2
for(ww in w[1:10]){ 
	finalSnv[[ww]] -> v
	v <- v[seqnames(v)==7]
	v <- v[which(info(v)$MajCN <= 4 & info(v)$MajCN > 1)]
	n <- sum(info(v)$MutCN==info(v)$MajCN, na.rm=TRUE)
	plotSample(ww, title=paste0(sample2icgc[names(finalSnv)[ww]], ", n=", n,"/",nrow(v), " SNVs pre +7"))
}

#' Tabulate the number of point mutations preceding +7 
t <- t(sapply(w, function(ww){
					finalSnv[[ww]] -> v
					v <- v[seqnames(v)==7]
					v <- v[which(info(v)$MajCN <= 4 & info(v)$MajCN > 1)]
					n <- sum(info(v)$MutCN==info(v)$MajCN, na.rm=TRUE)
					c(pre_7=n, total=nrow(v))
				}))
rownames(t) <- names(finalSnv)[w]

#' Less than 3 early
table(t[,1]<=3)

#' ### Extended Data Figure 3b
#' Plot the distributions of early and toal point mutations on chr7 across all GBM samples with +7:
#+ GBM_tri7_bee, fig.width=1.5, fig.height=2
.par()
par(bty="n")
beeswarm::beeswarm(as.numeric(pmax(t,0.5)) ~ factor(rep(c("pre +7","total"), each=nrow(t))), method='hex', pch=19, xlab="", ylab="Number of SNVs", cex=0.5, log=TRUE, yaxt='n', col=set1[c(3,9)])
axis(side=2, at=c(1,10,100,1000))
axis(side=2, at=0.5, label=0)
for( i in 0:2) axis(side=2, at=c(2,3,4,5,6,7,8,9)*10^i, labels=rep("",8), tcl=-0.1)
#dev.off()

#' Generate xlsx output:
#+ ExtendedDataFigure3b
ExtendedDataFigure3b <- createSheet(ExtendedDataFigure3, "ExtendedDataFigure3b")
addDataFrame(t, ExtendedDataFigure3b)

#' ## Medulloblastoma
#' Medulloblastoma showed very frequent and early i17q. Gather them all:
w <- which(donor2type[sample2donor[names(finalSnv)]]=="CNS-Medullo") 
w <- w[sapply(finalBB[w], function(bb) sum(width(bb)[as.logical(seqnames(bb)==17) & bb$total_cn >= 3], na.rm=TRUE)/width(refLengths[17])>0.5)]

#' ### Extended Data Figure 3c
#' Plot first 10 samples
#pdf("Medullo_i17q.pdf",3.2,3.2, pointsize=8)
for(ww in w[1:10]){
	finalSnv[[ww]] -> v
	v <- v[seqnames(v)==17]
	v <- v[which(info(v)$MajCN <= 4 & info(v)$MajCN > 1)]
	n <- sum(info(v)$MutCN==info(v)$MajCN, na.rm=TRUE)
	plotSample(ww, title=paste0(sample2icgc[names(finalSnv)[ww]], ", n=", n,"/",nrow(v), " SNVs pre i17q"))
}

#' Tabulate early and late substitutions. 
t <- t(sapply(w, function(ww){
					finalSnv[[ww]] -> v
					v <- v[seqnames(v)==17]
					v <- v[which(info(v)$MajCN <= 4 & info(v)$MajCN > 1)]
					n <- sum(info(v)$MutCN==info(v)$MajCN, na.rm=TRUE)
					c(pre_i17q=n, total=nrow(v))
				}))
rownames(t) <- names(finalSnv)[w]

#' Less than 1 early
table(t[,1]<=1)

#' ### Extended Data Figure 3d
#' Plot distribution across samples
#+ Medullo_i17q_bee, fig.width=1.5, fig.height=2
.par()
par(bty="n")
beeswarm::beeswarm(as.numeric(pmax(t,0.5) + runif(length(t))*0.25 -0.1) ~ factor(rep(c("pre i17q","total"), each=nrow(t))), method='hex', pch=19, xlab="", ylab="Number of SNVs", cex=0.5, log=TRUE, yaxt='n', col=set1[c(3,9)], corral="gutter", corralWidth=1.25)
axis(side=2, at=c(1,10,100,1000))
axis(side=2, at=0.5, label=0)
for( i in 0:2) axis(side=2, at=c(2,3,4,5,6,7,8,9)*10^i, labels=rep("",8), tcl=-0.1)

#' Dump xlsx output
#+ ExtendedDataFigure3d
ExtendedDataFigure3d <- createSheet(ExtendedDataFigure3, "ExtendedDataFigure3d")
addDataFrame(t, ExtendedDataFigure3d)
saveWorkbook(ExtendedDataFigure3,'ExtendedDataFigure3.xlsx')


#' # Synchronous gains
#' One thing observed in the temporal classification of WGD samples was that a surprisingly high fraction of near-diploid samples also 
#' display very concordant timings. Here we provide a basic assessment of the phenomenon. The final results were further refined, more accurately
#' considering the exact copy numeber configurations involved.
#'  
#' ## Synchronous/asynchronous classification
#' Note: Final figures have slightly deviated from this earlier version due to more elaborate handling of segments.
d <- fracGenomeWgdComp
i <- d[,"avg.ci"]<=0.5 & d[,"chr.all"] > 2 #&  fracGenomeWgdComp[,"nt.total"]/chrOffset["MT"] >= 0.1
timingClass <- paste(ifelse(isWgd,"WGD","ND"), ifelse(!i, "uninformative",""))
timingClass[i] <- paste0(timingClass[i], ifelse(d[i,"nt.wgd"]/d[i,"nt.total"] > 0.75,"sync","async"))
#timingClass[i] <- paste0(timingClass[i], cut(fracGenomeWgdComp[i,"nt.wgd"]/fracGenomeWgdComp[i,"nt.total"], c(0,0.5,0.8,1), include.lowest=TRUE))
timingClass <- factor(timingClass)

#' ### Figure 1f
#' First the pie chart with the timing classification
#+ timingClass, fig.width=4, fig.height=4
#pdf("TimingClass.pdf", 4,4)
colTime <- c("#A0C758","#6B8934","#BEC6AD","#CEB299","#CC6415","#EF7B00")
names(colTime) <- levels(timingClass)[c(4,5,6,3,2,1)]
c <- c(RColorBrewer::brewer.pal(9, "Pastel1"),"#DDDDDD")
t <- table(timingClass)[names(colTime)]
pie(t, init.angle=90, labels=paste0(names(t), ",\nn=", t), col=colTime)
#t <- table(isWgd)
par(new=TRUE)
symbols(x=0,y=0,circles=0.4, inches=FALSE, add=TRUE, bg="white")
#pie(t, labels=c("",""), col=NA, lwd=5, lty=1, init.angle=90)
#dev.off()

#' xlsx output
Figure1f <- createSheet(Figure1, "Figure1f")
addDataFrame(t, Figure1f)

colnames(d) <- c("ntCoamp","ntAmp","timeCoamp","segCoamp","segAmp","chrCoamp","chrAmp", "sdTimeCoamp","avgCiSeg","sdAllSeg")
timingInfo <- data.frame(avgPloidy=finalPloidy, avgHom=finalHom, isWgd=isWgd, d, informative=i, timingClass=timingClass)
#tab <- rbind(tab, data.frame(WGD_call=otherStat, WGD_timing=NA, ploidy=otherPloidy, hom=otherHom, nt.wgd=NA, nt.total=NA, time.wgd=NA, sd.wgd=NA,avg.ci=NA, sd.all=NA))
write.table(file=paste0(Sys.Date(),"-Timing-info.txt"), timingInfo, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")


#' ## Examples
#' Show 9 prototypical examples for the different timing classes.
#' 
#' ### Figure 1e
#+ timingExamples, fig.width=4, fig.height=4, warning=FALSE
w <- which(wgdStar=="likely" & !isWgd) # Garden variety near-diploid with concordant timing
plotSample(w[1]) 
plotSample(w[2])
plotSample(w[3])

w <- which(wgdStar=="very likely" & isWgd) # Now a-star WGD with highly concordant timing.
plotSample(w[1])
plotSample(w[2])
plotSample(w[9])

w <- which(wgdStar=="unlikely" & !isWgd & fracGenomeWgdComp[,"nt.total"]/chrOffset["MT"] > 0.25 & fracGenomeWgdComp[,"avg.ci"] < 0.5) # Near diploid with highly discordant timing
plotSample(w[1])
plotSample(w[2])
plotSample(w[3])


#' ## Relationship with mutation rates
#' As a sanity check for the molecular timing estimates, calculate number of substitutions and timing per tumour type. There shouldn't be a trend.
n <- nSub <- sapply(finalSnv, nrow)
n[timingInfo$timeCoamp==0] <- NA
q <- unlist(sapply(split(n, donor2type[sample2donor[names(finalSnv)]]), function(x) as.numeric(cut(x, {if(sum(!is.na(x))>1) quantile(x, seq(0,1,0.1), na.rm=TRUE) else 1:10}, include.lowest=TRUE))))
m <- match(names(finalSnv),unlist(split(names(finalSnv), donor2type[sample2donor[names(finalSnv)]])))
t <- timingInfo$timeCoamp
table(decSub=q[m], time=cut(t, seq(0,1,0.1)))

#' Also calculate deciles of timing per tumour type; this is an even stronger indication of independence, as hoped.
t[t==0] <- NA
r <- unlist(sapply(split(t, donor2type[sample2donor[names(finalSnv)]]), function(x) as.numeric(cut(x, {if(sum(!is.na(x))>1 & length(unique(x)) > 2) quantile(jitter(x), seq(0,1,0.1), na.rm=TRUE) else 1:10}, include.lowest=TRUE))))
table(decSub=q[m], decTime=r[m])

#' Plot the number of mutations per sample vs the average duplication time.
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


#' ## Secondary gains
#' We can also time secondary gains for certain configurations. This is always possible as long only one allele is gained, such as 3:1. If two alleles
#' are gained this required additional assumptions. For 3:2 the typical assumption would be 1:1 -> 2:2 -> 3:2. For 4:2 the assumption is 1:1 -> 2:1 -> 4:2.
#' Here we only consider gains of a single allele, which can be uniquely timed. These analyses were conducted by Lara Jerman. Here we load her data.
#' 
#' ### Load data
#' Load preprocessed data, aggregated by chromsome
load("two_gain_times.RData")
doubleGains <- as.data.frame(T.i.F)
m <- paste(doubleGains$sample, doubleGains$cnMaj, doubleGains$cnMin, doubleGains$chromosome, sep="_")
s <- split(doubleGains[,c("sample","tumor_type","T1_raw","T2_raw","n_mutations")], m)
doubleGainsAggregated <- Reduce("rbind",sapply(s, function(x) {
					data.frame(sample=x$sample[1], tumor_type=x$tumor_type[1], T1_raw=weighted.mean(x$T1_raw, x$n_mutations),T2_raw=weighted.mean(x$T2_raw, x$n_mutations), n_mutations=sum(x$n_mutations))
				}, simplify=FALSE))

#' ### Timing of primary and secondary gains
#' First plot all samples and chromosomes, first versus second gain. 
#+ multiGain, fig.height=2, fig.width=2
x <- doubleGainsAggregated$T1_raw/pmax(1, doubleGainsAggregated$T2_raw)
y <- doubleGainsAggregated$T2_raw/pmax(1, doubleGainsAggregated$T2_raw)
o <- order(doubleGainsAggregated$n_mutations, decreasing=TRUE)
plot(x[o], 
		y[o], 
		bg=tissueColors[as.character(donor2type[sample2donor[names(finalSnv)[doubleGainsAggregated$sample[o]]]])], pch=21,
		col=tissueBorder[as.character(donor2type[sample2donor[names(finalSnv)[doubleGainsAggregated$sample[o]]]])],
		xlab="Time [mutations], first gain",
		ylab="Time [mutations], second gain",
		cex=sqrt(doubleGainsAggregated$n_mutations[o]/500)+0.1,
		lwd=0.5)
t <- table(doubleGainsAggregated$sample)

#' The distribution is surprisingly independent, as one might have expected multiple copies of the same allele to be gained at the same time. 
#' 
#' Now split by timing class - no obvious trend.
#+ multiGainClass, fig.height=3, fig.width=3
names(timingClass) <- names(finalSnv)
par(mfrow=c(2,2))
for(l in grep("uninformative",levels(timingClass), invert=TRUE, value=TRUE)){
	w <- which(timingClass[doubleGains$sample]==l)
	o <- intersect(order(doubleGainsAggregated$n_mutations, decreasing=TRUE),w)
	plot(x[o], 
			y[o], 
			bg=tissueColors[as.character(donor2type[sample2donor[names(finalSnv)[doubleGainsAggregated$sample[o]]]])], pch=21,
			col=tissueBorder[as.character(donor2type[sample2donor[names(finalSnv)[doubleGainsAggregated$sample[o]]]])],
			xlab="Time [mutations], first gain",
			ylab="Time [mutations], second gain",
			cex=sqrt(doubleGainsAggregated$n_mutations[o]/500)+0.1)
	title(main=l, line=0, font.main=1)
}

#' Clusters are often caused by individual samples:
#+ multiGainSamples, fig.height=4, fig.width=4
par(mfrow=c(5,5))
for(i in as.numeric(names(t)[t>5])[1:25]){
	w <- which(doubleGainsAggregated$sample==i)
	plot(x[w],y[w], 
			col=tissueBorder[as.character(donor2type[sample2donor[names(finalSnv)[doubleGainsAggregated$sample[w]]]])], 
			bg=tissueColors[as.character(donor2type[sample2donor[names(finalSnv)[doubleGainsAggregated$sample[w]]]])], 
			type='p', xlim=c(0,1), ylim=c(0,1), 
			xlab="time 1",
			ylab="time 2",
			pch=21,
			cex=sqrt(doubleGainsAggregated$n_mutations[w]/500+0.1))
	
}

#' ### Relative latency
#' To formalise the observation that the second gains occur seemingly uniformly after the first, calculate and plot the relative latency.
#+ multiGainLatency, fig.height=1, fig.width=1.5
w <- y < 1 & x > 0
r <- ((y-x)/(1-x))
h <- hist(r[w], breaks=seq(0,1,0.025), plot=FALSE)
e <- d <- density(r[w], from=-1, to=2, bw="SJ")
i <- which(d$x < 0)
d$y[max(i) + seq_along(i)] <-  d$y[max(i) + seq_along(i)] + d$y[rev(i)]
i <- which(d$x > 1)
d$y[min(i) -  seq_along(i)] <-  d$y[min(i) -  seq_along(i)] + d$y[i]
i <- which(d$x >0 & d$x < 1)
d <- list(x=d$x[i], y=d$y[i])
plot(h$mids,h$counts/sum(h$counts),  pch=19, col='grey',ylim=c(0,max(h$counts/sum(h$counts))), xlab="Latency of second gain", ylab="Relative frequency", type='h')
#' The majority of secondary gains indeed appear to occur independently after the second. Only about 17% occur in close succession.

plot(d$x,cumsum(d$y * diff(d$x)[1]), xlim=c(0,1), type='l', ylim=c(0,1), xlab="Relative time of second gain", ylab="CDF")

#' #### Figure 1h
#' It is perhaps instructive to assess the relative latency by timing class. WGD appears underrepresented in cases with synchronous double gains.
#+ multiGainLatencyClass, fig.height=1, fig.width=1.5
c <- cut(r[w], 20)
t <- table(timingClass[doubleGainsAggregated$sample[w]],c)
barplot(t[names(colTime),]/sum(t), border=NA, col=colTime, width=1/24, space=0.2, names.arg=rep("",20, bty="L", yaxs="s"))
.par()
axis(side=1, line=0.2)

#' xlxs output
#+ Figure1h
Figure1h <- createSheet(Figure1, "Figure1h")
addDataFrame(t[names(colTime),], Figure1h)

#' Lastly, simulate higher order gains to cross-check whether the inference is credible.
n <- 100
c <- 40
purity=0.7
bb <- GRanges(seqnames=1, IRanges(1,1e8), major_cn=3, minor_cn=1, clonal_frequency=purity, n.snv_mnv=n)
bb$total_cn <- bb$major_cn+bb$minor_cn
t3 <- 0.8
t2 <- 0.8 # Simultaneous second amplification
d <- data.frame(cluster=1, n_ssms=1, proportion=purity)
cnStates <- defineMcnStates(bb,clusters=d, purity=purity)
pi <- matrix(c(4,2,1,0,1,0,0,0,1), byrow=TRUE, ncol=3) %*% c(1-t2,t2-t3,t3)
pi <- pi/sum(pi)
cnStates[[1]][,"P.m.sX"] <- pi
rho=0.01
cnStates[[1]][,"power.m.s"] <- 1-pbetabinom(2, size=c, prob=cnStates[[1]][,"f"], rho=rho )

bb$timing_param <- cnStates
s <- simplify2array(mclapply(1:50, function(foo){
					set.seed(foo)
					v <- simulateMutations(bb, n=40)
					bb0 <- bb
					bb0$timing_param <- NULL
					L <- computeMutCn(v, bb0, clusters=d, purity=purity, rho=rho, n.boot=0)
					L$P[[1]]
				}))

#+ simMultiGain, fig.height=2.5, fig.width=1.5
boxplot(t(s[2:3,"T.m.sX",]), at=3:2, xlab="Simulated time point", names=c("t2","t3"))
points(3:2,c(t2-t3, t3), col='red', pch=19)

#+ simMultiGainLatency, fig.height=1, fig.width=1.5
l <- s[2,"T.m.sX",]/(1-s[3,"T.m.sX",])
x <- seq(0,1,0.05)
plot(x[-1]+x[2]/2, as.numeric(prop.table(table(cut(l[l<1], x)))), xlab="Latency", ylab="frequency", type='h')
axis(side=1)

#' #### Figure 1g
#' We also calculate how often higher order gains occur in each timing class. To this end, to a tally of all copy number gains. 
cn <- do.call("rbind", sapply(names(finalBB), function(n){
					bb <- finalBB[[n]]
					data.frame(sample=n, chr=seqnames(bb), width=width(bb), M=bb$major_cn, m=bb$minor_cn)}, simplify=FALSE))

#' How does this relate to the length of the segments?
#+ distSegCopies, fig.width=2.5, fi.height=2.5
t <- table(pmin(cn$M,3) ,  pmax(3,round(log10(cn$width),1)), timingClass[cn$sample])
x <- as.numeric(colnames(t))
plot(NA,NA, type='p', col=colTime[1], pch=16, ylim=c(0,0.8), xlim=range(10^x), xlab="Length of segment", ylab="Proportion >2 allelic copies", log='x')
for(n in dimnames(t)[[3]]) {
	y <- as.numeric(t[4,,n]/colSums(t[3:4,,n]))
	lines(10^x,y, type='p', col=paste0(colTime[n],"44"), pch=16, cex=1)#sqrt(colSums(t[3:4,,i]/1000)))
	lines(10^x, predict(loess(y ~x)), col=colTime[n], lwd=2)
}

#' Plot the fraction of samples with double gains for each timing class. Asynchronous samples appear overrepresented.
#+ fracDoubleGains, fig.width=1.5, fig.height=1
tt <- mg14:::asum(t[,x>=7,],2)
o <- names(colTime)
p <- tt[4,o]/colSums(tt[3:4,o])
ci <- sapply(c(0.025, 0.975), qbeta, 0.025, shape1=tt[4,o]+1, shape2=tt[3,o]+1)
barplot(p, col=colTime, border=NA, las=2, ylab="Proportion >2 allelic copies", names=sub("ormative","",sub("near-diploid", "ND", names(colTime))), ylim=c(0,0.4)) -> b
segments(b, ci[,1], b, ci[,2])

#' Produce xlsx output
#+ Figure1g
Figure1g <- createSheet(Figure1, "Figure1g")
addDataFrame(t(tt[3:4,o]), Figure1g)




#' # Chronological WGD & MRCA
#' Here we calculate approximate real time estimates of WGD and MRCA using only CpG>TpG mutations, which are found in every tumoyr type 
#' and vary relatively little between samples.
#' 
#' ## CpG>TpG branch lengths and mutation rates
#' The approximate chronological timing is basd only CpG>TpG mutations, which universally occur in all tissues, are used as a precautionary measure.
#' We first study the burden of CpG>TpG mutations (RpCpG>RpTpG in Melanoma) in each sample across tissues and as a function of age. This enables us to 
#' remove some hypermutant samples and derive a range of possible mutation rate increase.
#' 
#' Branch lengths need to be adjusted for time-dependent copy number as well as the structure of the subclonal phyllogeny. 
#' 
#' ### Prelim
#' Get the age of every donor.
age <- clinicalData$donor_age_at_diagnosis
names(age) <- clinicalData$icgc_donor_id
#' Exclude these cancer types:
typeNa <- gsub("\t","",strsplit("Bone-Cart
						Breast-LobularCa
						Breast-DCIS
						Lymph-NOS
						Myeloid-MDS
						Cervix-AdenoCa", "\n")[[1]])
#' 
#' #### Effective (time-averaged) genome size
#' Calculate effective genome size, i.e. time-averaged ploidy from mutation copy numbers. This is useful to convert mutation counts into approximate rates.
#+ effGenome
effGenome <- unlist(mclapply(finalSnv, function(vcf) {
					w <- info(vcf)$CLS!="subclonal"
					if(donor2type[sample2donor[meta(header(vcf))$META["ID",]]]=="Skin-Melanoma")
						w <- w & isDeaminationNoUV(vcf)
					else
						w <- w & isDeamination(vcf)
					2/avgWeights(vcf[na.omit(w)])
				}))
names(effGenome) <- names(finalSnv)

#' #### Power per (sub)clone
#' Also calculate the power to detect variants for different subclones to extrapolate their true branch length.
#+ finalPower
finalPower <- sapply(names(finalBB), function(n) {
			x <- finalBB[[n]]
			f <- finalClusters[[n]]$proportion
			for(i in 1:length(x)){
				t <- x$timing_param[[i]]
				p <- t[match(f, t[,"cfi"]), "power.s"]
				if(!is.null(p)) if(all(!is.na(p))) break
			}
			if(is.null(p)) return(rep(NA, length(f)))
			return(p)
		})
#plot(unlist(lapply(wccClusters[names(finalSnv)], `[[`, "n_ssms")),unlist(lapply(finalClusters[names(finalSnv)], `[[`, "n_ssms"))/ unlist(finalPower), log='xy',
#		xlab="Cluster size WCC (consensus)", ylab="Cluster size MutationTime.R") 

#' #### Branch lengths
#' The following calculates the length of the trunk (clonal mutations) and the depth of the MRCA, scaled by power and using a densitly branching subclonal phylogeny.
#' Scaling branches by 1/f assumes the densest possible branching topology and this the shortest latency.
#+ branchDeam
branchDeam <- t(simplify2array(mclapply(finalSnv, function(vcf) {
							n <- meta(header(vcf))$META["ID",]
							if(donor2type[sample2donor[n]]=="Skin-Melanoma")
								w <- isDeaminationNoUV(vcf)
							else
								w <- isDeamination(vcf)
							if(sum(w)==0) return(c(0,0))
							p <- info(vcf)$pSub[w]; 
							n.subclonal <- aggregate(p, list(info(vcf)$CNF[w]), sum, na.rm=TRUE)
							m <- apply(abs(outer(n.subclonal$Group.1, finalClusters[[n]]$proportion, `-`)),1,which.min) # Match to subclones
							p.subclonal <- finalPower[[n]][m] # Power of subclones
							b.subclonal <- n.subclonal$x %*% (n.subclonal$Group.1 / p.subclonal) / max(n.subclonal$Group.1) # Subclonal branch, power adjusted & 1/f-scaled
							b.clonal <- sum(1-p, na.rm=TRUE)/finalPower[[n]][1] # Clonal branch (trunk), power adjusted & 1/f-scaled
							c(b.subclonal, b.clonal)})))

d <- droplevels(donor2type[sample2donor[names(finalSnv)]])
typesSubclones <- setdiff(levels(d), c(typeNa, names(which(table(d)<5))))

nClones <- sapply(finalClusters, nrow)

#' Comparison to linear branching in which the branch lengths of all subclonal populations is added. This constitutes the other extreme, corresponding to
#' the maximal latency.
#+ branchDeamLinear
branchDeamLinear <- t(simplify2array(mclapply(finalSnv, function(vcf) {
							if(donor2type[sample2donor[meta(header(vcf))$META["ID",]]]=="Skin-Melanoma")
								w <- isDeaminationNoUV(vcf)
							else
								w <- isDeamination(vcf)
							if(sum(w)==0) return(c(0,0))
							n <- meta(header(vcf))$META["ID",]
							if(donor2type[sample2donor[n]]=="Skin-Melanoma")
								w <- isDeaminationNoUV(vcf)
							else
								w <- isDeamination(vcf)
							if(sum(w)==0) return(c(0,0))
							p <- info(vcf)$pSub[w]; 
							n.subclonal <- aggregate(p, list(info(vcf)$CNF[w]), sum, na.rm=TRUE)
							m <- apply(abs(outer(n.subclonal$Group.1, finalClusters[[n]]$proportion, `-`)),1,which.min) # Match to subclones
							p.subclonal <- finalPower[[n]][m] # Power of subclones
							b.subclonal <- n.subclonal$x %*% (1 / p.subclonal)  # Subclonal branch, power adjusted 
							b.clonal <- sum(1-p, na.rm=TRUE)/finalPower[[n]][1] # Clonal branch (trunk), power adjusted
							c(b.subclonal, b.clonal)})))


#' Compare branching and linear topologies:
#+ branchingLinear
f <- (branchDeam[,1] / finalPloidy) / rowSums(branchDeam / cbind(finalPloidy, effGenome))
l <- (branchDeamLinear[,1]/ finalPloidy)/rowSums(branchDeamLinear / cbind(finalPloidy, effGenome))
t <- donor2type[sample2donor[names(finalSnv)]]
plot(f, l, xlab="Subclonal branch length (branching)", ylab="Subclonal branch length (linear)", pch=21, bg=tissueColors[t], col=tissueBorder[t], cex=tissueCex[t])
abline(0,1, lty=2)

quantile(l/f, na.rm=TRUE)
quantile(l, na.rm=TRUE)
quantile(f, na.rm=TRUE)

#' ### Average mutation rates
#' Here we convert branch lengths into average mutation rates (assuming no change, ie avgmu=n/age). 
#' First study relation of mutation burden to age, then calculate average mutation rates,
#' exclude hypermutators and also samples with tumour in normal 1% for further analyses. 
#' 
#' Plot a large overview:
#+ timeSubcloneAge, fig.width=10, fig.height=10
rateDeam <- cc <- list()
remove <- "8454fe53-869d-41c8-b0c8-a7929d00eec3" # a cell line, add more samples in the following
par(mfrow=c(6,6), mar=c(3,3,2,1),mgp=c(2,.5,0), tcl=-0.25,cex=1, bty="L", xpd=FALSE, las=1, xpd=FALSE)
for(n in typesSubclones){
	i <- d==n
	tt0 <- branchDeam[i,]/cbind(finalPloidy[i], effGenome[i]) / 3#cbind(nClones[i]-1, 1)/3 # 3Gb Haploid genome
	tt0[is.infinite(tt0)|is.nan(tt0)] <- 0
	yy <- rowSums(tt0)
	a <- age[sample2donor[names(finalSnv)[i]]]
	xx <- a
	r <- yy/xx 
	m <- median(r[TiN[names(xx)] <= 0.01 & ! is.na(TiN[names(xx)])],na.rm=TRUE)
	rateDeam[[n]] <- r
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
n <- names(rateDeam)
qRateDeam <- sapply(rateDeam, function(r){
			m <- median(r[TiN[sample2donor[names(r)]] <= 0.01 & ! is.na(TiN[sample2donor[names(r)]])],na.rm=TRUE)
			w <- (r-m)^2/m^2 <= 2^2 & TiN[sample2donor[names(r)]] <= 0.01 & ! is.na(TiN[sample2donor[names(r)]])
			quantile(r[w], na.rm=TRUE)})
plot(sapply(rateDeam, median, na.rm=TRUE), pch=NA , ylab="SNVs/Gb/yr", main="CpG>TpG rate", ylim=c(0, max(qRateDeam)), cex.main=1, xaxt='n', xlab="Tumour type")
segments(seq_along(rateDeam),qRateDeam["0%",],seq_along(rateDeam), qRateDeam["100%",], col=tissueLines[n], lty=1)
points(sapply(rateDeam, median, na.rm=TRUE), pch=21, col=tissueBorder[n], bg=tissueColors[n])

length(remove)

#' Plot the average rates as barplot across tissues. Relatively little variation.
#+deamRateBar, fig.width=4, fig.height=2
par(mar=c(6,3,1,1))
o <- order(qRateDeam["50%",])
barplot(qRateDeam["50%",][o], col=tissueColors[colnames(qRateDeam)][o], border=tissueLines[colnames(qRateDeam)][o], las=2,names.arg=rep("",ncol(qRateDeam)) , ylab="CpG>TpG rate [SNVs/Gb/yr]", ylim=c(0, max(qRateDeam))) -> b
mg14::rotatedLabel(b, labels=colnames(qRateDeam)[o])
segments(b, qRateDeam["50%",][o], b, qRateDeam["100%",][o], col=tissueLines[colnames(qRateDeam)][o], lwd=2)
segments(b, qRateDeam["0%",][o], b, qRateDeam["50%",][o], col=tissueBorder[colnames(qRateDeam)][o], lwd=2)


#' #### Extended Data Figure 8a
#' Plot the copynumber and ITH-adjusted mutation burden versus age and add a loess fit.
#+ timeSubcloneAgePancan, fig.width=2, fig.height=2
tt0 <- branchDeam/cbind(finalPloidy, effGenome) / cbind(nClones-1, 1)/3 # 3Gb Haploid genome
tt0[is.infinite(tt0)|is.nan(tt0)] <- 0
m <- sapply(rateDeam, function(r){m <- median(r[TiN[sample2donor[names(r)]] <= 0.01 & ! is.na(TiN[sample2donor[names(r)]])],na.rm=TRUE)})
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

#' xlsx
#+ ExtendedDataFigure8a
ExtendedDataFigure8a <- createSheet(ExtendedDataFigure8, "ExtendedDataFigure8a")
tab <- data.frame(age=x, `CpG>TpG/Gb`=s, tumour_type=t)
addDataFrame(tab, ExtendedDataFigure8a)

#' Is there a positive intercept in a linear mutation v age fit?
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

#' The two quantities above can be used to calculate a fraction of mutations due to linear accumulation
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
#' 
#' As the OLS fits per cancer type suggest a robust linear trend with small positive offsets, we repeat the analysis using a hierarchical Bayesian model.
#' This enables us to share some information across cancer types, implicitly encode the positivity constraints and also propagate uncertainty into 
#' functions of the inferred parameters, such as the fraction of mutations due to linear accumulation.
#' 
#' Prepare data for `rstan`
library(rstan)
y <- Reduce("c",rateDeam)
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

#' The actual stan model definition is defined in `./PCAWG-rates.stan`:
#+ PCAWG-rates.stan, eval=FALSE

#' Now fit the model using Hamiltonian Monte Carlo
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
s <- summary(fit, pars=c("offset","slope"))$summary
ab <- array(s, dim=c(ncol(tt),2,10), dimnames=list(levels(droplevels(t)), c("a","b"), colnames(s)))

#' Summary plot
#+ deamAgeBayes, fig.width=2, fig.height=2
plot(x,y, bg=tissueColors[t], pch=21, ylim=c(0,1000), col=tissueBorder[t], cex=tissueCex[t]*2/3, lwd=0.25, xlab="Age", ylab="SNVs/Gb")
for(i in 1:nrow(ab))
	abline(ab[i,1,"50%"], ab[i,2,"50%"], col=tissueLines[levels(droplevels(t))[i]], lty=tissueLty[levels(droplevels(t))[i]])

#' #### Extended Data Figure 8c
#' Plot inferred rates and offsets with 95% credible intervals.
#+rateOffsetBayes, fig.width=2, fig.height=2
plot(ab[,1,"50%"], ab[,2,"50%"], col=tissueColors[dimnames(ab)[[1]]], pch=NA, xlab="Offset", ylab="SNVs/Gb/yr", xlim=range(ab[,1,c("2.5%","97.5%")]), ylim=range(ab[,2,c("2.5%","97.5%")]))
segments(ab[,1,"50%"], ab[,2,"2.5%"],ab[,1,"50%"], ab[,2,"97.5%"], col=tissueLines[dimnames(ab)[[1]]], pch=19)
segments(ab[,1,"2.5%"], ab[,2,"50%"],ab[,1,"97.5%"], ab[,2,"50%"], col=tissueLines[dimnames(ab)[[1]]], pch=19)
points(ab[,1,"50%"], ab[,2,"50%"], pch=21, bg=tissueColors[dimnames(ab)[[1]]], col=tissueLines[dimnames(ab)[[1]]])
abline(h=0, lty=3)
abline(v=0, lty=3)

a <- extract(fit, pars="offset")$offset
b <- extract(fit, pars="slope")$slope
colnames(a) <- colnames(b) <- levels(droplevels(t))

#' xlsx ouput
#+ ExtendedDataFigure8c
ExtendedDataFigure8c <- createSheet(ExtendedDataFigure8, "ExtendedDataFigure8c")
addDataFrame(ab, ExtendedDataFigure8c)


#' #### Extended Data Figure 8b
#' Create a large overview panel with all data points used and the Bayesian linear fits with 95% CIs.
#+ timeSubcloneAgeBayes, fig.width=10, fig.height=10
par(mfrow=c(6,6), mar=c(3,3,2,1),mgp=c(2,.5,0), tcl=-0.25,cex=1, bty="L", xpd=FALSE, las=1, xpd=FALSE)
d <- droplevels(t)
for(n in levels(d)){
	i <- d==n
	yy <- y[i]
	xx <- x[i]
	m <- median(yy/xx, na.rm=TRUE)
	try({
				plot(xx, yy, bg=tissueColors[n], col=tissueBorder[n], pch=21, log='', xlab="Age at diagnosis", ylab="SNVs/Gb", main=n, ylim=c(0,pmin(1000,max(yy, na.rm=TRUE))), xlim=c(0,max(x, na.rm=TRUE)),  cex.main=1)
				#points(xx, yy, bg=tissueColors[n], col=ifelse(w,tissueBorder[n], tissueColors[n]), pch=ifelse(w,21,4))
				abline(0, m, lty=3)
				x0 <- seq(0,100,1)
				p <- apply(sapply(x0, function(x) a[,n] + b[,n]*x), 2, quantile, c(0.025,0.975), na.rm=TRUE)
				polygon(c(x0,rev(x0)), c(p["2.5%",], rev(p["97.5%",])), border=tissueBorder[n], col=paste0(tissueColors[n],"44"))
			})
}

#' xlsx outputs.
#+ ExtendedDataFigure8b
ExtendedDataFigure8b <- createSheet(ExtendedDataFigure8, "ExtendedDataFigure8b")
tab <- data.frame(age=x, `CpG>TpG/Gb`=y, tumour_type=d)
addDataFrame(tab, ExtendedDataFigure8b)

#' Calculate the fraction of mutations due to linear accumulation using the inferred parameters and propagate the error from the HMC samples.
q <- sapply(colnames(a), function(n){
			w <- which(t==n & !is.na(y))
			f <- sapply(w, function(j) x[j] * b[,n] / (a[,n] + x[j] * b[,n]))
			quantile(rowMeans(f), c(0.025, 0.25, .5,.75,.975))
		})*100

qPanCan=quantile(rowMeans(do.call("cbind",sapply(colnames(a), function(n){
									w <- which(d==n & !is.na(y))
									f <- sapply(w, function(j) x[j] * b[,n] / (a[,n] + x[j] * b[,n]))
								}))),
		c(0.025, 0.25, .5,.75,.975))*100



#' #### Extended Data Figure 8e
#' Plot the fraction of mutations due to linear accumulation using bar plots with 95% CIs, sorted.
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

#' xlsx output
#+ ExtendedDataFigure8e
ExtendedDataFigure8e <- createSheet(ExtendedDataFigure8, "ExtendedDataFigure8e")
addDataFrame(t(q), ExtendedDataFigure8e)


#' Quickly assess whether this fraction qualitative agrees with simulations using a 5x range increase occuring between 0-15yrs prior 
#' to sampling.
set.seed(42)
x <- df$x
a <-  runif(length(x), pmax(0.5, 1-15/x), 1) #
r <- rgamma(length(x), 10, 10)
y <- rpois(length(x), (a + (1-a)*5) * r * x * 6 * 0.2)/6
y[x < 40] <- NA

#+fracLinearSim, fig.width=2, fig.height=2
plot(x,y, ylim=c(0,max(y, na.rm=TRUE)), xlab="Age", ylab="SNVs/Gb", pch=21, bg="grey", lwd=0.5)
f <- lm(y~x)
summary(f)
c <- coef(f)
abline(c)

mean(c[1] / (c[1]+ x*c[2]), na.rm=TRUE)

#' ## CpG>TpG acceleration in relapse samples
#' Load a table with data from cancer relapse samples in 8 tissues. This data was collected from 8 separate studies and
#' compiled by Santiago Gonzalez.
relapse <- read.table("../data/2018_05_acceleration.txt", header=TRUE, sep="\t")
head(relapse)

#' ### Branch lengths
#' Now caculate branch lengths to visualise these data, adjusted for (time-dependent) ploidy. To this end first the ploidy in the primary branch.
p <- relapse$Primary_ploidy
p[is.na(p)] <- 2

#' And for the MRCA. Impute effective genome if timing is missing using typical scale factor for effective genome (this is just for graphical purposes).
m <- p
m[m>2.8 & is.na(relapse$Primary_effGenome)] <- m[m>2.8 & is.na(relapse$Primary_effGenome)] * 0.7 
m[m>2.8 & !is.na(relapse$Primary_effGenome)] <- relapse$Primary_effGenome[m>2.8 & !is.na(relapse$Primary_effGenome)]

#' Ploidy in the relapse branch
r <- relapse$Relapse_ploidy
r[is.na(r)] <- relapse$Primary_ploidy[is.na(r)]
r[is.na(r)] <- 2

#' Branch lengths, adjusted for ploidy
b_mrca <- relapse$Primary_CpG_muts * relapse$shared / 3 /m
b_relapse <- (relapse$Relapse_CpG_muts - relapse$Primary_CpG_muts * relapse$shared)/ 3/ r
b_primary <- (relapse$Primary_CpG_muts - relapse$Primary_CpG_muts * relapse$shared) / 3 /p

#' Corresponding time points
t_primary <- relapse$Primary_age
t_relapse <- t_primary + relapse$Relapse_time
t_mrca <- 0.5 * (t_primary * b_mrca/(b_mrca + b_primary) + pmin(t_primary, t_relapse * b_mrca/(b_mrca + b_primary) ))

#' #### Extended Data Figure 8f,g
#+ relapse_accel_branch, fig.width=4, fig.height=2
#pdf("relapse_accel_branch.pdf", pointsize=8, width=4, height=2)
.par()
par(mfrow=c(1,2))
for(t in c("Ovarian", "Breast")){
	w <- which(relapse$Ttype==t)
	r <- relapse[w,]
	plot(r$Primary_age+r$Relapse_time, r$Relapse_CpG_muts, ylim=c(0,150), xlim=c(0,80), pch=NA, xlab="Age", ylab="CpG>TpG/Gb")
	tt <- c(Ovarian="Ovary-AdenoCa",Breast="Breast-AdenoCa")[t]
	title(main=tt, line=0)
	c <- tissueColors[tt]
	#mrca <- r$Primary_CpG_muts*r$shared
#points(mrca_time, mrca)
	segments(0,0, t_mrca[w],b_mrca[w], col=c)
	segments(t_mrca[w], b_mrca[w], t_primary[w], (b_mrca+b_primary)[w], col=c)
	segments(t_mrca[w], b_mrca[w], t_relapse[w], (b_mrca+b_relapse)[w], col=c)
	points(t_primary[w], b_mrca[w]+b_primary[w], pch=21, bg=c, col="white")
	points(t_relapse[w], (b_mrca+b_relapse)[w], pch=21, bg=c, col="white")
}
#dev.off()

ExtendedDataFigure8f <- createSheet(ExtendedDataFigure8, "ExtendedDataFigure8f")
tab <- data.frame(t_mrca, t_primary, t_relapse, b_mrca, b_primary, b_relapse, type=relapse$Ttype, ID=relapse$ID)
addDataFrame(tab[which(relapse$Ttype=="Ovarian"),], ExtendedDataFigure8f)
ExtendedDataFigure8g <- createSheet(ExtendedDataFigure8, "ExtendedDataFigure8g")
addDataFrame(tab[which(relapse$Ttype=="Breast"),], ExtendedDataFigure8g)

#' ### Mutation rate increase
#' The acceleration - or mutation rate increase - is effectively the ratio of the ploidy adjusted CpG>TpG branch-length in the relase
#' sample versus the MRCA. See Supplementary Methods for details.
#' 
#' #### Extended Data Figure 8h
#' Plot the acceleration values inferred from CpG>TpG mutations in each branch
#+ accRelapse, fig.width=2, fig.height=2
#pdf('foo.pdf', 2,2, pointsize=8)
.par()
par(mar=c(5,3,1,1))
x <- relapse$Ttype
levels(x) <- c("Myeloid-AML","Breast-AdenoCa","CNS-LGG","Liver-HCC","Lung-AdenoCa","CNS-Medullo","Lymph-BNHL","Ovary-AdenoCa")
x <- as.character(x)
c <- tissueColors; c["CNS-LGG"] <- "black"
b <- tissueBorder; b["CNS-LGG"] <- "white"
t <- tissueCex; t["CNS-LGG"] <- 1
j <- jitter(as.numeric(factor(x)))
plot(NA,NA, xlim=c(.5,8), xaxt="n", xlab="", ylab="Acceleration", ylim=c(0,15))
segments(j, relapse$Acc_min, j, relapse$Acc_max, col=c[x])
points(Acc_CpG ~ j, data=relapse, bg=c[x], col=b[x], pch=21, cex=t[x])
abline(h=1, lty=3)
mg14::rotatedLabel(labels=levels(factor(x)))
#dev.off()

ExtendedDataFigure8h <- createSheet(ExtendedDataFigure8, "ExtendedDataFigure8h")
addDataFrame(relapse, ExtendedDataFigure8h)


#' The surplus of mutations in relapse samples (after adjusting for the depth of branching) shows a very good correlation with the time difference
#' between primary and relapse sample
#+ relapse_time, fig.width=2, fig.height=2 
#pdf("relapse_time.pdf", pointsize=8, width=2, height=2)
.par()
plot(t_relapse - t_primary, b_relapse-b_primary, xlim=c(0,13), xlab="Time to relapse", ylab="Surplus in relapse, CpG>TpG/Gb", pch=16)
w <- t_relapse - t_primary < 10
abline(lm((b_relapse - b_primary)[w] ~ (t_relapse-t_primary)[w]))
legend("topright", c(paste("rho =",round(cor(t_relapse - t_primary, b_relapse-b_primary, method='s'),2)),
				paste("p =",signif(cor.test(t_relapse - t_primary, b_relapse-b_primary, method='s')$p.value,1))), bty="n")
#dev.off()


#' ## Normal tissue mutation rates
#' Here we compare the infereed CpG>TpG baseline rates to data from four different studies on normal tissues.
#' 
#' ### Normal blood CSC
#' Data of blood colonies from a 59yr-old individual. Taken from Lee-Six et al. Nature 2018
hsc_data <- read.table(gzfile("../data/Shearwater_calls_FDR0.95_all_muts.txt.gz"), header=TRUE, sep="\t")
ct <- which((hsc_data$REF=="C" & hsc_data$ALT=="T") | (hsc_data$REF=="G" & hsc_data$ALT=="A" ))
v <- VRanges(seqnames = hsc_data$X.CHROM, ref=hsc_data$REF, alt=hsc_data$ALT, ranges=IRanges(hsc_data$POS, width=1))
tnc=scanFa(file=refFile, resize(granges(v), 3,fix="center"))
cpgtpg <- grepl("(A.CG)|(T..CG)", paste(alt(v),tnc))
n_cpgtpg <- colSums(hsc_data[cpgtpg,5:144], na.rm=TRUE)
normal_hsc_cpgtpg <- quantile(n_cpgtpg/59/6, c( 0.5, 0.025,0.975))

#' ### Normal colon
#' Colonic crypts of several individuals. Data taken from Lee-Six et al. bioRxiv 2018, https://www.biorxiv.org/content/10.1101/416800v1.
#+ normalColon
colon_sbs <- read.table("../data/model_input_with_CtoTatCpG.txt", header=TRUE, sep="\t")
foo <- as.data.frame(summary(lm(CtoTatCpG/6 ~ age-1, data=colon_sbs))$coef)
normal_colon_cpgtpg <- quantile(colon_sbs$CtoTatCpG/colon_sbs$age/6, c( 0.5, 0.025,0.975))#c(foo$Estimate, foo$Estimate - 2*foo$`Std. Error`, foo$Estimate + 2*foo$`Std. Error`)
plot(colon_sbs$age, colon_sbs$CtoTatCpG/6, xlab="Age", ylab="CpG>TpG/Gb", pch=16)

#' ### Endometrium
#' Normal endometrial glands. Data from Moore et al., bioRxiv 2018, https://www.biorxiv.org/content/10.1101/505685v1.
#+ normalEndometrium
tab <- read.table("../data/endom_subs.txt", sep="\t", header=TRUE)
quantile(tab$C.T.at.CpG/tab$Age, c(0.5, 0.025, 0.975))/6
normal_endometrium_cpgtpg <- quantile(tab$C.T.at.CpG/tab$Age, c(0.5, 0.025, 0.975))/6
plot(tab$Age, tab$C.T.at.CpG/6, xlab="Age", ylab="CpG>TpG/Gb", pch=16)

#' ### Normal skin
#' Data for a single normal skin biopsy from Martincorena et al., Science 2015
foo <- read.table("../ref/PD20399be_wg_caveman_annotated_with_dinucs_for_mg14.txt", header=TRUE, sep="\t")
is_cpgtpg <-  grepl("(A.CG[C,T])|(T.[A,G]CG)", paste(foo$mut,foo$trinuc_Ref))
normal_skin_cpgtpg <- sum(is_cpgtpg * foo$VAF)/0.375/55/6 # Adjust for YCG fraction, age and genome

normal_cpgtpg <- rbind(`Myeloid-MPN`=normal_hsc_cpgtpg,
		`Skin-Melanoma`=c(normal_skin_cpgtpg,NA,NA) ,
		`Uterus-AdenoCa`=normal_endometrium_cpgtpg,
		`ColoRect-AdenoCa`=normal_colon_cpgtpg)

x <- abind::abind(cancer=ab[rownames(normal_cpgtpg), "b",colnames(normal_cpgtpg)], normal=normal_cpgtpg[,], along=3)
x["Skin-Melanoma",,] <-x["Skin-Melanoma",,]


#' ### Extended Data Figure 8d
#' Side-by-side comparison of the median rates of normal tissues with linear (ie presumed baseline) rates from cancer samples with 95% CI.
#+cancer_normal_cpgtpg.pdf, fig.width=2.5, fig.height=2.5
#pdf("cancer_normal_cpgtpg.pdf", 2.5,2.5, pointsize=8)
par(mar=c(3,3,1,1), bty="L", mgp=c(2,.5,0), tcl=-0.25, las=1)
t(barplot(t(x[,"50%",]), beside=TRUE, col=rep(tissueColors[rownames(normal_cpgtpg)], each=2), density=rep(c(NA,36), nrow(normal_cpgtpg)), ylim=c(0,max(x, na.rm=TRUE)),
				ylab="Mutation rate [CpG>TpG/Gb/yr]", names.arg=c("Blood","Skin","Uterus","Colon"))) -> b
legend("topleft", c("Cancer","Normal"), fill="black", density=c(NA,36), bty="n")
#foo <- do.call("rbind", list(
#				data.frame(cpgtpg=n_cpgtpg/59/6, x=b[1,2]), 
#				data.frame(cpgtpg=normal_skin_cpgtpg, x=b[2,2]),
#				data.frame(cpgtpg=tab$C.T.at.CpG/tab$Age/6, x=b[3,2]),
#				data.frame(cpgtpg=colon_sbs$CtoTatCpG/colon_sbs$age/6, x=b[4,2])))
#points(cpgtpg ~ jitter(x,.5), data=foo, cex=.5, pch=16, col="#00000044")
segments(b,x[,"2.5%",],b,x[,"97.5%",])
#dev.off()

#' xlsx output
#+ ExtendedDataFigure8d
ExtendedDataFigure8d <- createSheet(ExtendedDataFigure8, "ExtendedDataFigure8d")
addDataFrame(x, ExtendedDataFigure8d)

#' ## MRCA
#' Now we time the latency of the most recent common ancestor (MRCA), or more accurately we compute the time between MRCA and the last observable subclone.
#' An implicit assumption is that that the latency between the emergence of the last detectable subclone and diagnosis is short with respect
#' to the time between fertilisation and diagnosis. 
#' 
#' ### Functions
#' Acceleration values to simulate
accel <- c(1,2.5,5,7.5,10,20)
names(accel) <- paste0(accel, "x")

#' The actual timing using the age of each individual. This assumes that the branch length of the subclonal period is to be devided by the 
#' acceleration (ie mutation rate increase). We propagate Poisson error to obtain CIs. 
#+ timeSubclones, warning=FALSE
set.seed(42)
d <- droplevels(donor2type[sample2donor[names(finalSnv)]])

computeSubclonesTimeAbs <- function(l, b) {
	i <- d==l
	tt0 <- b[i,]/cbind(finalPloidy[i], effGenome[i]) #/ cbind(nClones[i]-1, 1)
	resB <- sapply(1:1000, function(foo){ ## Assess the impact of Poisson fluctuations on numbers
				tt <- matrix(rpois(length(tt0), lambda=tt0), ncol=ncol(tt0))
				res <- sapply(accel, function(a)  tt[,1]/a/rowSums(tt/rep(c(a,1), each=nrow(tt)))) * age[sample2donor[names(finalSnv)[i]]]
				colnames(res) <- paste0(accel, "x")
				#res[res==0] <- NA
				res}, simplify='array')
	res <- sapply(accel, function(a)  tt0[,1]/a/rowSums(tt0/rep(c(a,1), each=nrow(tt0)))) * age[sample2donor[names(finalSnv)[i]]]
	colnames(res) <- paste0(accel, "x")	
	resCI <- apply(resB,1:2, quantile, c(0.1,0.9), na.rm=TRUE)
	arr <- abind::abind(res, resCI, along=1)
	rownames(arr)[1] <- "hat"
	arr <- aperm(arr, c(2,1,3))
	tt0[is.infinite(tt0)|is.nan(tt0)] <- 0
	r <- which(rowSums(b[i,]) < 50 ) ## Exclude samples with less than 50 subs 
	arr[r,,] <- NA
	return(arr)
}

#' ### Timing
#' Calculate for linear and branching topologies
subclonesTimeAbs <- mclapply(typesSubclones, computeSubclonesTimeAbs, b=branchDeam)
subclonesTimeAbsLinear <- mclapply(typesSubclones, computeSubclonesTimeAbs, b=branchDeamLinear)
names(subclonesTimeAbsLinear) <- names(subclonesTimeAbs) <- typesSubclones

#' Buest guess acceleration values (usually 5x, except Ovarian and CNS).
guessAccel <- sapply(subclonesTimeAbs, function(x) "5x")
guessAccel["Ovary-AdenoCa"] <- guessAccel["Liver-HCC"] <- "7.5x"
guessAccel[grep('CNS', names(guessAccel))] <- "2.5x"

#' #### Extended Data Figure 9f
#' Plot the inferred latencies per sample with overlaid boxplots.
#+ realTimeSubclone, fig.width=6, fig.height=2.225
u <- setdiff(names(finalSnv)[uniqueSamples], remove)
par( mar=c(7,3,1,1), mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
qSubclone <- sapply(subclonesTimeAbs, function(x) apply(x[,"hat",][rownames(x)%in%u,,drop=FALSE], 2, quantile, c(0.05,0.25,0.5,0.75,0.95), na.rm=TRUE), simplify='array')
a <- "5x"
subclonesTimeAbsType <- sapply(names(subclonesTimeAbs), function(n) {x <- subclonesTimeAbs[[n]]; x[,,guessAccel[n]][setdiff(rownames(x),remove), 1:3, drop=FALSE]})
m <- diag(qSubclone["50%",guessAccel[dimnames(qSubclone)[[3]]],])#t[1,3,]
names(m) <- dimnames(qSubclone)[[3]]
nSubclones <- sapply(subclonesTimeAbsType, function(x) sum(!is.na(x[,1])))
m[nSubclones < 5] <- NA
o <- order(m, na.last=NA)
plot(NA,NA, xlim=c(0.5,length(m[o])), ylab="Years before diagnosis", xlab="", xaxt="n", yaxs="i", ylim=c(0,30))
abline(h=seq(10,20,10), col="#DDDDDD", lty=3)
x <- seq_along(m[o])
mg14::rotatedLabel(x, labels=names(sort(m)))
for(i in seq_along(o))try({
				n <- names(m)[o[i]]
				f <- function(x) x/max(abs(x))
				a <- guessAccel[n]
				bwd <- 0.8/2
				j <- if(length(na.omit(subclonesTimeAbsType[[o[i]]][,"hat"]))>1) f(mg14::violinJitter(na.omit(subclonesTimeAbsType[[o[i]]][,"hat"]))$y)/4 + i else i
				tpy <- 2
				segments(j, na.omit(subclonesTimeAbsType[[o[i]]][,"90%"]), j, na.omit(subclonesTimeAbsType[[o[i]]][,"10%"]), col=mg14::colTrans(tissueLines[n],tpy))
				points(j, na.omit(subclonesTimeAbsType[[o[i]]][,"hat"]), pch=21, col=mg14::colTrans(tissueBorder[n], tpy), bg=mg14::colTrans(tissueColors[n],tpy), 
						cex=tissueCex[n]*2/3, lwd=1)
				rect(i-bwd,qSubclone["25%",a,n],i+bwd,qSubclone["75%",a,n], border=tissueLines[n],  col=paste(tissueColors[n],"44", sep=""))
				segments(i-bwd,qSubclone["50%",a,n],i+bwd,qSubclone["50%",a,n],col=tissueLines[n], lwd=2)
				segments(i,qSubclone["75%",a,n],i,qSubclone["95%",a,n],col=tissueLines[n], lwd=1.5)
				segments(i,qSubclone["5%",a,n],i,qSubclone["25%",a,n],col=tissueLines[n], lwd=1.5)
				f <- function(x) x/max(abs(x))
			})

#par(xpd=TRUE)
#s <- 12/8
#dev.copy2pdf(file="realTimeSubclone.pdf", width=6*s, height=3.5*3/5*s, pointsize=8*s)#' xlsx

#' xlsx output
#+ ExtendedDataFigure9f
ExtendedDataFigure9f <- createSheet(ExtendedDataFigure9, "ExtendedDataFigure9f")
tab <- data.frame(do.call("rbind", subclonesTimeAbsType))
tab$tumour_type <- donor2type[sample2donor[rownames(tab)]]
addDataFrame(tab, ExtendedDataFigure9f)



#' #### Figure 5c
#' Alternative representation of median acceleration across the range of simulated rate increases, per cancer type:
#+ realTimeSubcloneMed, fig.width=6, fig.height=2.225
#pdf("realTimeSubcloneMed.pdf", width=6, height=2.225, pointsize=8)
u <- setdiff(names(finalSnv)[uniqueSamples], remove)
m <- qSubclone["50%","5x",]#t[1,3,]
names(m) <- dimnames(qSubclone)[[3]]
m[nSubclones < 5] <- NA
o <- order(m, na.last=NA)
n <- dimnames(qSubclone)[[3]]
par( mar=c(7,3,1,1), mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
plot(NA,NA, xlim=c(0.5,length(m[o])), ylab="Latency [yr]", xlab="", xaxt="n", yaxs="i", ylim=c(0,15))
abline(h=seq(10,20,10), col="#DDDDDD", lty=3)
x <- seq_along(m[o])
mg14::rotatedLabel(x, labels=names(sort(m)))
b <- .3
rect(seq_along(o)-b,qSubclone["50%","1x",o],seq_along(o)+b,qSubclone["50%","10x",o], col=paste(tissueColors[n[o]],"88", sep=""), border=1)
rect(seq_along(o)-b,qSubclone["50%","2.5x",o],seq_along(o)+b,qSubclone["50%","7.5x",o], col=paste(tissueColors[n[o]],"FF", sep=""), border=1)
rect(seq_along(o)-b,qSubclone["50%","20x",o],seq_along(o)+b,qSubclone["50%","10x",o], col=paste(tissueColors[n[o]],"22", sep=""), border=1)
segments(seq_along(o)-b,qSubclone["50%","5x",o],seq_along(o)+b,qSubclone["50%","5x",o])

#par(xpd=TRUE)
#s <- 12/8
#dev.copy2pdf(file="realTimeSubclone.pdf", width=6*s, height=3.5*3/5*s, pointsize=8*s)

sapply(subclonesTimeAbs, nrow)

#' xlsx output
#+ Figure5c
Figure5c <- createSheet(Figure5, "Figure5c")
addDataFrame(aperm(qSubclone, c(3,1,2)), Figure5c)

#' #### Extended Data Figure 9g
#' Comparison of branching v linear
#+ realTimeBranchLinear, fig.width=2.5, fig.height=3.5
subclonesTimeAbsTypeLinear <- sapply(names(subclonesTimeAbsLinear), function(n) {x <- subclonesTimeAbsLinear[[n]]; x[,,guessAccel[n]][setdiff(rownames(x),remove), 1:3, drop=FALSE]})
qSubcloneLinear <- sapply(subclonesTimeAbsLinear, function(x) apply(x[,"hat",][rownames(x)%in%u,,drop=FALSE], 2, quantile, c(0.05,0.25,0.5,0.75,0.95), na.rm=TRUE), simplify='array')
n <- diag(qSubcloneLinear["50%",guessAccel[dimnames(qSubcloneLinear)[[3]]],])#t[1,3,]
names(n) <- dimnames(qSubcloneLinear)[[3]]
par( mar=c(5,3,3,10), mgp=c(2,.5,0), tcl=-0.25,cex=1, bty="n", xpd=FALSE, las=1)
plot(c(rep(1, length(m)), rep(2, each=length(n))), c(m,n), bg=tissueColors[c(names(m), names(n))], pch=21, cex=1, xaxt="n", ylab="Years after MRCA", xlab="", xlim=c(0.5,2.5), ylim=c(0, max(n, na.rm=TRUE)))
segments(rep(1, each=length(m)), m, rep(2, each=length(n)), n,col=tissueLines[names(m)], lty= ifelse(sapply(subclonesTimeAbsType, nrow) <= 5, 3, tissueLty[names(m)]))
o <- order(n, na.last=NA)
y0 <- n[o]
y1 <- mg14:::mindist(n[o], diff(par('usr')[3:4])/30)
par(xpd=NA)
mtext(names(m)[o], at=y1, side=4 )
segments(2.1,y0,2.2,y0)
segments(2.2,y0,2.3,y1)
segments(2.3,y1,2.4,y1)
mg14::rotatedLabel(1:2, labels=c("Branching","Linear"))


#' xlsx output
#+ ExtendedDataFigure9g
ExtendedDataFigure9g <- createSheet(ExtendedDataFigure9, "ExtendedDataFigure9g")
addDataFrame(data.frame(branching=m, linear=n), ExtendedDataFigure9g)

#' Numbers per decade
yy <- do.call("rbind",subclonesTimeAbsType)
yy <- yy[setdiff(rownames(yy), remove),"hat"]
table(cut(yy, seq(0,60,10)))

#' ## WGD
#' The other landmark event that can be timed in approximately 30% of samples are whole-genome doublings.
#' 
#' ### Functions
#' Calculate relative timing estimates based on deaminations. The following function calculates 
#' molecular timing estimates by only evaluating CpG>TpG mutations, separately for aggreagated 
#' 2:0, 2:1, and 2:2 regions in each sample.

computeWgdParamDeam <- function(vcf, bb, clusters, purity){
	# 1. Find segments compatible with WGD
	min.dist <- 0.05
	m <- findMainCluster(bb)
	l <- pmin(bb$time.lo, bb$time - min.dist)
	u <- pmax(bb$time.up, bb$time + min.dist)
	o <- which(l <= m & u >= m)
	
	# 2. Find deaminations in compatible segments
	w <- which(info(vcf)$MajCN==2 & sapply(info(vcf)$CNID, length)==1 & isDeamination(vcf) & vcf %over% bb[o])
	if(donor2type[sample2donor[meta(header(vcf))$META["ID",]]]=="Skin-Melanoma")
		w <- intersect(w, which(isDeaminationNoUV(vcf)))
	v <- vcf[w]
	if(nrow(v)<=100) return(NULL) # At least 100 SNVs
	seqnames(rowRanges(v)) <- factor(3-info(v)$MinCN, levels=seqlevels(v))
	v <- sort(v)
	
	# 3. Merged CN segments
	b <- GRanges(1:3, IRanges(rep(1,3),rep(max(end(v)),3)), copy_number=4:2, major_cn=2, minor_cn=2:0, clonal_frequency=as.numeric(purity))
	
	# 4. Calculate times
	l <- computeMutCn(v, b, clusters, purity, isWgd=TRUE, n.boot=200, rho=0.01, xmin=3)
	b$n.snv_mnv <- l$n <- table(factor(info(v)$MinCN, levels=2:0))
	l$time <- bbToTime(b, l$P)
	return(l)
}

#' ### Timing
#' #### Molecular timing estimates
#' Takes ~1h.
#+ finalWgdParam, eval=FALSE
wgdParamDeam <- mclapply(names(finalSnv)[isWgd], function(ID){
			try(computeWgdParamDeam(finalSnv[[ID]], finalBB[[ID]], clusters=finalClusters[[ID]], purity=finalPurity[ID]))
		})
names(wgdParamDeam) <- names(finalSnv)[isWgd]

#+ finalWgdParamLoad, echo=FALSE
wgdParamDeam <- readRDS("2018-06-13-finalWgdParam.rds")

#' Samples with insufficient data
void <- sapply(wgdParamDeam, function(x) is.null(x) | class(x)=="try-error")

#' Some checks
t <- sapply(wgdParamDeam[!void], function(x) {r <- as.matrix(x$time[,2:4]); rownames(r) <- x$time[,1];r}, simplify='array')
pairs(t(t[,"time",]))

#' #### Adjust for mutation rate increase
#' Now transform the molecular timing estimates (effectively branch lengths, scaled to 1) accounting for possinle late mutation 
#' rate increase of up to 20x. This returns a large array.
wgdTimeDeamAcc <- simplify2array(mclapply(names(wgdParamDeam[!void]), function(n) {
					x <- wgdParamDeam[!void][[n]]
					
					T.clonal <- as.matrix(x$time[,2:4]) # Time of WGD as fraction of clonal					
					
					n.subclonal <- aggregate(x$D[,"pSub"], list(x$D[,"CNF"]), sum)
					m <- match(n.subclonal$Group.1, finalClusters[[n]]$proportion)
					p.subclonal <- x$power.c[m] # Power of subclones
					b.subclonal <- n.subclonal$x %*% (n.subclonal$Group.1 / p.subclonal) / max(n.subclonal$Group.1) # Subclonal branch, power adjusted & 1/f-scaled
					b.clonal <- sum(1-x$D[,"pSub"])/p.subclonal['1'] # Clonal branch (trunk), power adjusted & 1/f-scaled
					f.subclonal <- b.subclonal / (b.subclonal + b.clonal)
					G.clonal <- sum (1-x$D$pSub)/sum((1-x$D$pSub)*x$D$MutCN/(x$D$MajCN + x$D$MinCN)) # Effective ploidy clonal, adjusted for timing
					G.subclonal <- sum(x$D$pSub*(x$D$MajCN + x$D$MinCN))/ sum (x$D$pSub) # Final ploidy
					if(is.nan(G.subclonal)) G.subclonal <- mean(x$D$MajCN + x$D$MinCN)
					
					ag <- age[sample2donor[names(finalBB)[isWgd][!void][j]]]
					tmin <- max(0.5, 1-15/ag) # 15yrs or 50%, whatever smaller (median ~ 0.75 mutation time)
					if(is.na(tmin)) tmin <- 0.8
					ta=seq(tmin,1,l=20)
					
					.correctAccel <- function(T.clonal, f.subclonal, G.clonal, G.subclonal, ta, a){ # Helper function to correct accel a at clonal time ta
						t1 <- T.clonal + (1-T.clonal) *(a-1)/a*ta #acc before dup
						t2 <- T.clonal * (ta + a*(1-ta)) ## after
						T.clonal.adj <- pmin(t1, t2) # as fraction of clonal
						a.clonal <- ta + (1-ta)*a # effective rate, avg over clonal
						T.subclonal.abs <- f.subclonal / G.subclonal / a
						T.clonal.abs <- (1-f.subclonal) / G.clonal/ a.clonal
						T.clonal.abs <- T.clonal.abs / (T.clonal.abs + T.subclonal.abs) # as fraction of all mutations
						return(c(T.WGD=T.clonal.adj * T.clonal.abs, T.MRCA=T.clonal.abs))
					}
					
					.correctAccelRand <- function(T.clonal, f.subclonal, G.clonal, G.subclonal, ta=seq(0.8,1,0.01), a=seq(1,10,1)){ # Helper to calculate range of accel a and times
						sapply(ta, function(taa) sapply(a, function(aa) .correctAccel(T.clonal, f.subclonal, G.clonal, G.subclonal, taa, aa)), simplify='array')
					}
					
					res <- apply(T.clonal, 1:2, .correctAccelRand, f.subclonal, G.clonal, G.subclonal, a=accel, ta=seq(tmin,1,l=20))
					dim(res) <- c(2, length(accel), length(ta), dim(res)[-1])
					return(res)
					
				}))
dimnames(wgdTimeDeamAcc)[1:2] <- list(c("T.WGD","T.MRCA"), names(accel))
dimnames(wgdTimeDeamAcc)[[5]] <- colnames(wgdParamDeam[[1]]$time)[2:4] 
dimnames(wgdTimeDeamAcc)[[4]] <- levels(finalBB[[1]]$type)[c(3,1,2)]
dimnames(wgdTimeDeamAcc)[[6]] <- names(wgdParamDeam[!void])

n <- dimnames(wgdTimeDeamAcc)[[6]]
d <- droplevels(donor2type[sample2donor[n]])
s <- setdiff(levels(d), c(typeNa, names(which(table(d)<3))))

#' #### Convert to chronological time
#' Lastly calculate the corresponding real time by scaling with age at diagnosis
#+ wgdTimeAbs, warning=FALSE
f <- function(n, mu, a, b){ ## asymmetric normal to interpolate CIs of the timing estimates
	r <- rnorm(n, sd=a)
	w <- which(r>0)
	r[w] <- r[w]*(b/a)[w]
	return(r + mu)
}
wgdTimeAbs <- sapply(s, function(l) {
			set.seed(42)
			i <- d==l & ! n %in% c(rownames(purityPloidy)[purityPloidy$wgd_uncertain])
			
			## absolute time by multiplying with age at diagnosis
			absTimeSeg <- aperm((1-wgdTimeDeamAcc["T.WGD",,,,,i]) * rep(age[sample2donor[n]][i], each = prod(dim(wgdTimeDeamAcc)[c(2,3,4,5)])))
			w <- t(sapply(wgdParamDeam[n[i]], `[[`, "n")) #number of SNV as weights
			
			## remove NA due to zero mutations
			for(c in 1:3)
				absTimeSeg[,,c,,][is.na(absTimeSeg[,,c,,]) & w[,c]==0] <- 0
			
			## weighted average over 2+0, 2+1 and 2+2 segments
			absTime <- (absTimeSeg[,,1,,] * w[,1] + absTimeSeg[,,2,,] * w[,2] + absTimeSeg[,,3,,] * w[,3]) / rowSums(w) 
			rownames(absTime) <- n[i]
			
			## Median over acceleration onset
			absTimeMed <- apply(absTime, c(1,2,4), median, na.rm=TRUE) 
			colnames(absTimeMed) <- c("hat","up","lo")
			
			## Simulate distribution sampling from timing onset and mutation time CI 
			ts <- sapply(1:1000, function(foo) {ax <- sample(1:20,1); matrix(f(length(absTime[,"time",ax,]),a=abs(absTime[,"time",ax,] - absTime[,"time.up",ax,])/2, b=abs(absTime[,"time.lo",ax,]-absTime[,"time",ax,])/2, mu=absTime[,"time",ax,]), nrow=dim(absTime)[1])}, simplify='array')
			me <- apply(ts, 1:2, quantile, c(0.1, 0.8), na.rm=TRUE) # 80% CIs
			absTimeMed[,"lo",] <- (me[1,,])
			absTimeMed[,"up",] <- (me[2,,])
			absTimeMed[,c(1,3,2),]
		}, simplify=FALSE)


#' #### Extended Data Figure 9a
#' First plot the inferred latency of WGD for each sample, per cancer type, overlaid by boxplots, using buest-guess acceleration.
#+ realTimeWgd, fig.height=3, fig.width=4.5
par( mar=c(7,3,1,1), mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
u <- setdiff(names(finalSnv)[uniqueSamples], remove)
qWgd <- sapply(wgdTimeAbs, function(x) apply(x[rownames(x) %in% u,"hat",], 2, quantile, c(0.05,0.25,0.5,0.75,0.95), na.rm=TRUE), simplify='array')
nWgd <- sapply(wgdTimeAbs, function(x) sum(x[rownames(x) %in% u,"hat","1x"]!=0, na.rm=TRUE))
wgdTimeAbsType <- sapply(names(wgdTimeAbs), function(n) {x <- wgdTimeAbs[[n]]; x[,,guessAccel[n]][setdiff(rownames(x),remove), 1:3, drop=FALSE]})
m <- diag(qWgd["75%",guessAccel[dimnames(qWgd)[[3]]],])#t[1,3,]
names(m) <- dimnames(qWgd)[[3]]
o <- order(m, na.last=NA)
x <- seq_along(m[o])
plot(NA,NA, xlim=c(0.5,length(m[o])), ylim=c(0,max(do.call('rbind',wgdTimeAbsType)[,1], na.rm=TRUE)+5), ylab="Years before diagnosis", xlab="", xaxt="n", yaxs="i")
abline(h=seq(10,60,10), col="#DDDDDD", lty=3)
mg14::rotatedLabel(x, labels=names(sort(m)))
for(i in seq_along(o)){
	n <- names(m)[o[i]]
	f <- function(x) x/max(abs(x))
	a <- guessAccel[n]
	j <- f(mg14::violinJitter(na.omit(wgdTimeAbsType[[o[i]]][,"hat"]))$y)/4 + i #rank(na.omit(tWgdByType[[o[i]]][,"hat"]))/2/length(na.omit(tWgdByType[[o[i]]][,"hat"]))-0.25+i #
	tpy <- if(grepl("Skin|Lung", n)) 4 else 2
	segments(j, na.omit(wgdTimeAbsType[[o[i]]][,"up"]), j, na.omit(wgdTimeAbsType[[o[i]]][,"lo"]), col=mg14::colTrans(tissueLines[n],tpy), lty=tissueLty[n])
	points(j, na.omit(wgdTimeAbsType[[o[i]]][,"hat"]), pch=21, col=mg14::colTrans(tissueBorder[n], tpy), bg=mg14::colTrans(tissueColors[n],tpy), 
			cex=tissueCex[n]*2/3, lwd=1)
	bwd <- 0.8/2
	rect(i-bwd,qWgd["25%",a,n],i+bwd,qWgd["75%",a,n], border=tissueLines[n],  col=paste0(tissueColors[n],"44"))
	segments(i-bwd,qWgd["50%",a,n],i+bwd,qWgd["50%",a,n],col=tissueLines[n], lwd=2)
	segments(i,qWgd["75%",a,n],i,qWgd["95%",a,n],col=tissueLines[n], lwd=1.5)
	segments(i,qWgd["5%",a,n],i,qWgd["25%",a,n],col=tissueLines[n], lwd=1.5)
}
par(xpd=TRUE)
#s <- 12/8
#dev.copy2pdf(file="realTimeWgd.pdf", width=4*s, height=3.5*s, pointsize=8*s)

#' xlsx output
#+ ExtendedDataFigure9a
ExtendedDataFigure9a <- createSheet(ExtendedDataFigure9, "ExtendedDataFigure9a")
tab <- data.frame(do.call("rbind", wgdTimeAbsType))
tab$tumour_type <- donor2type[sample2donor[rownames(tab)]]
addDataFrame(tab, ExtendedDataFigure9a)

#' #### Figure 5b
#' Plot the median timing value for WGD, as a function of the unknown rate increase.
#+ realTimeWgdMed, fig.height=3, fig.width=4.5
#pdf("realTimeWgdMed.pdf", width=4.5, height=3, pointsize=8)
u <- setdiff(names(finalSnv)[uniqueSamples], remove)
m <- qWgd["50%","5x",]#t[1,3,]
names(m) <- dimnames(qWgd)[[3]]
m[nWgd < 5] <- NA
o <- order(m, na.last=NA)
n <- dimnames(qWgd)[[3]]
par( mar=c(7,3,1,1), mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
plot(NA,NA, xlim=c(0.5,length(m[o])), ylab="Latency [yr]", xlab="", xaxt="n", yaxs="i", ylim=c(0,40))
abline(h=seq(10,20,10), col="#DDDDDD", lty=3)
x <- seq_along(m[o])
mg14::rotatedLabel(x, labels=names(sort(m)))
b <- .3
rect(seq_along(o)-b,qWgd["50%","1x",o],seq_along(o)+b,qWgd["50%","10x",o], col=paste(tissueColors[n[o]],"88", sep=""), border=1)
rect(seq_along(o)-b,qWgd["50%","2.5x",o],seq_along(o)+b,qWgd["50%","7.5x",o], col=paste(tissueColors[n[o]],"FF", sep=""), border=1)
rect(seq_along(o)-b,qWgd["50%","20x",o],seq_along(o)+b,qWgd["50%","10x",o], col=paste(tissueColors[n[o]],"22", sep=""), border=1)
segments(seq_along(o)-b,qWgd["50%","5x",o],seq_along(o)+b,qWgd["50%","5x",o])

#' xlsx output
#+ Figure5b
Figure5b <- createSheet(Figure5, "Figure5b")
addDataFrame(aperm(qWgd, c(3,1,2)), Figure5b)

#' As a quality control, plot a range of extremely early samples
#+ earlyWgdExamples, fig.width=4, fig.height=4
t <- do.call("rbind", wgdTimeAbsType)
o <- order(t[,"hat"], na.last=NA)
for(n in rownames(t)[tail(o, 20)])
	plotSample(n, title=paste0(sub("-.+","",n),", ", donor2type[sample2donor[n]], ", ",round(t[n,"hat"]),"yr"))

#' Numbers per decade
yy <- do.call("rbind",wgdTimeAbsType)
yy <- yy[setdiff(rownames(yy), remove),"hat"]
table(cut(yy, seq(0,60,10)))

#' #### Extended Data Figure 9b
#' A panel of scatter plots showing WGD time as a function of age at diagnosis. 
#+ wgdAge, fig.width=10, fig.height=10
par(mfrow=c(6,6), mar=c(3,3,2,1),mgp=c(2,.5,0), tcl=-0.25,cex=1, bty="L", xpd=FALSE, las=1)
for(i in seq_along(wgdTimeAbsType)){
	n <- names(wgdTimeAbsType)[i]
	y <- wgdTimeAbsType[[n]][,"hat"]
	x <- age[sample2donor[names(y)]]
	plot(x,x-y, pch=NA, bg=tissueColors[n], col=tissueBorder[n], xlab="Age [yr]", ylab="WGD [yr]", cex=tissueCex[n], xlim=c(0,90), ylim=c(0,90))
	segments(x, y0=x-wgdTimeAbsType[[n]][,"up"],y1=x-wgdTimeAbsType[[n]][,"lo"],col=tissueLines[n])
	points(x,x-y, pch=21, bg=tissueColors[n], col=tissueBorder[n], cex=tissueCex[n])
#	d <- density(na.omit(x), bw="SJ", from=0)
#	lines(d$x,d$y*100,col=tissueLines[n], lty=tissueLty[n])
#	d <- density(na.omit(x-y), bw="SJ", from=0)
#	lines(d$y*100, d$x,col=tissueLines[n], lty=tissueLty[n])
	rug(x, col=tissueLines[n],)
	rug(x-y, side=2, col=tissueLines[n])
	title(main=n, line=0, font.main=1, cex.main=1)
	abline(0,1, lty=3)
}

#' xlsx output
#+ ExtendedDataFigure9b
ExtendedDataFigure9b <- createSheet(ExtendedDataFigure9, "ExtendedDataFigure9b")
tab <- data.frame(do.call("rbind", wgdTimeAbsType))
tab$age <- age[sample2donor[rownames(tab)]]
tab$tumour_type <- donor2type[sample2donor[rownames(tab)]]
addDataFrame(tab, ExtendedDataFigure9b)

#' ### Deaminations v all mutations
#' For a rough comparison take the average co-amplification molecular time estimate calculated above against the estimate based on CpG>TpG mutations.
#+ wgdDeamAll, fig.width=2, fig.height=2
t <- sapply(wgdParamDeam[!void], function(x) x$time$time)
plot(timingInfo[colnames(t),'timeCoamp'], 
		colMeans(t), 
		cex=sqrt(nSub[colnames(t)]/10000), 
		pch=21,bg=tissueColors[donor2type[sample2donor[colnames(t)]]], 
		col=tissueBorder[donor2type[sample2donor[colnames(t)]]],
		xlab="time [all mutations]",
		ylab="time [CpG>TpG]")
abline(0,1, lty=3)

#' Plot difference as boxplots:
#+ wgdDeamAllBox, fig.width=6, fig.height=3
par(mar=c(6,4,1,1), xaxs='i')
boxplot(timingInfo[colnames(t),'timeCoamp'] - colMeans(t) ~ droplevels(donor2type[sample2donor[colnames(t)]]), 
		col=tissueColors[levels(droplevels(donor2type[sample2donor[colnames(t)]]))],
		na.rm=TRUE, las=2,
		lty=1,
		staplewex=0,
		pch=16,
		cex=0.8,
		ylab='time [all] - time [CpG>TpG]',
		names=NA
)
mg14::rotatedLabel(labels=levels(droplevels(donor2type[sample2donor[colnames(t)]])))
abline(h=0, lty=3)

#' ### Assessment of absolute mutation counts
#' To gain further insighe assess the distribution of the actual mutation numbers used for WGD timing across (2:0, 2:1, 2:2) regions
nDeam <- sapply(wgdParamDeam[!void], function(x) if(!is.null(x$n)) sum(x$n, na.rm=TRUE) else NA)

#' Number of deaminations only in 2:2 regions
nDeam22 <- sapply(wgdParamDeam[!void], function(x) if(!is.null(x$n)) as.numeric(x$n[1]) else NA)
w22 <- sapply(finalBB[isWgd][!void], function(bb) {
			w <- bb$major_cn==2 & bb$minor_cn==2 & !duplicated(bb)
			sum(as.numeric(width(bb)[w]), na.rm=TRUE)})
nDeam22 <- nDeam22/w22*1e9

#' Fraction of deam on 1 and 2 copies
d <- nDeam22 *  t(sapply(wgdParamDeam[!void], function(x) if(!is.null(x$P)) x$P[[1]][1:2,"P.m.sX"] else c(NA,NA)))
d[rownames(d) %in% remove] <- NA

#' Unadjusted time (inc. of subclonal mutations)
t0 <- colMeans(wgdTimeDeamAcc["T.WGD","1x",1,,"time",],na.rm=TRUE) # WGD
names(t0) <- dimnames(wgdTimeDeamAcc)[[6]]
m0 <- colMeans(wgdTimeDeamAcc["T.MRCA","1x",1,,"time",],na.rm=TRUE) # MRCA
names(m0) <- dimnames(wgdTimeDeamAcc)[[6]]

#' #### Figure 5a
#' Conceptual plot
#+ concept, fig.height=2, fig.width=2
#pdf("chronRateIncrease.pdf",2,2,pointsize=8)
.par()
a <- c(1,2.5,5,7.5,10, 20)
x <-  age[sample2donor["052665d1-ab75-4f40-be5a-b88154c8beed"]]
y <- nDeam["052665d1-ab75-4f40-be5a-b88154c8beed"]
t <- seq(0,x, by=0.1)
f <- function(t, ta, a){y <- t; y[t>ta] <- y[t>ta] + (y[t>ta] - min(y[t>ta]))*(a-1); y <- y/max(y)}
r <- sapply(a, function(aa){
			apply(sapply(x - seq(15,0,l=100), function(ta) f(t, ta, aa)),1,mean)
		})
plot(NA, NA, xlim=c(0,x), ylim=c(0,y), xlab="Age", ylab="CpG>TpG mutations")
for(i in seq_along(a)) lines(t, r[,i]*y)
w <- t0["052665d1-ab75-4f40-be5a-b88154c8beed"]
v <- m0["052665d1-ab75-4f40-be5a-b88154c8beed"]
segments(0,0, 0,w*y, col=set1[3]) 
text(10, w*y, "WGD", pos=3)
segments(0,w*y, 0, v*y, col=set1[4], lwd=2) # (WGD, MRCA)
text(0, v*y, "MRCA", pos=4)
segments(0,v*y, 0, y, col=set1[1], lwd=3) # (MRCA, diag)
for(i in seq_along(a)){
	x0 = t[which.min(abs(w - r[,i]))]
	segments(x0,0,x0,w*y, lty=3)
}
segments(0,w*y, x0, w*y, lty=3)
#dev.off()


#' #### Extended Data Figure 9c&d
#' Plot time v early and total number of mutations 
#+ nDeam22Time, fig.width=2, fig.height=2
y <- d[names(t0),]/6
x <- t0
t <- donor2type[sample2donor[names(t0)]]
plot(x,y[,2], bg=tissueColors[t], pch=21,  col=tissueBorder[t], cex=tissueCex[t]*1, lwd=0.5, xlab="Time", ylab="Early SNVs/Gb", log='')
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
#' Reassuringly, the molecular timing estimate appears independent of the total mutation count and largely driven by a depletion of early variants 
#' (rather than an excess of late variants).
#' 
#' xlsx outputs
#+ ExtendedDataFigure9c
ExtendedDataFigure9c <- createSheet(ExtendedDataFigure9, "ExtendedDataFigure9c")
tab <- data.frame(`Molecular Time`=x, `Early CpG>TpG/Gb`=y[,2])
addDataFrame(tab, ExtendedDataFigure9c)

#+ ExtendedDataFigure9d
ExtendedDataFigure9d <- createSheet(ExtendedDataFigure9, "ExtendedDataFigure9d")
tab <- data.frame(`Molecular Time`=x, `Total CpG>TpG/Gb`=rowSums(y))
addDataFrame(tab, ExtendedDataFigure9d)


#' ### Acceleration adjustment relative to lowest quintile.
#' So far we used a fixed rate increase across samples. If some samples experienced a higher increase, this should manifest as a relative excess of mutations.
#' Hence, calculate timing estimates using greater rate increase for samples with a higher CpG>TpG burden (per cancer type).
#+ accelRelWgd, fig.height=2, fig.width=2
accelRelWgd <- sapply(names(wgdTimeAbs), function(n) {
			r <- rateDeam[[n]]
			x <- wgdTimeAbs[[n]]
			r0 <- quantile(r[!names(r) %in% remove], 0.20, na.rm=TRUE)
			a <- cut(pmin(pmax(1,(r/r0-0.9)/0.1),10), c(0,1.5,3.75,6.25,8.75,20), labels=c("1x","2.5x","5x","7.5x","10x"))
			print(table(a))
			names(a) <- names(r)
			ta <- sapply(rownames(x), function(ss) x[ss, "hat",a[ss]])
			ta
		})
par(mar=c(3,3,1,1),mgp=c(2,.5,0), tcl=-0.25,cex=1, bty="L", xpd=FALSE, las=1)
x <- Reduce("c",accelRelWgd)
y <- Reduce("c",sapply(wgdTimeAbs, function(x) x[,"hat","5x"]))

t <- donor2type[sample2donor[names(x)]]
plot(y+runif(length(y)),x+runif(length(y)), pch=21, bg=tissueColors[t], col=tissueBorder[t], cex=tissueCex[t]*1, lwd=0.5, xlab="Time (constant acceleration)", ylab="Time (sample-specific accel.)", log='')
#s <- 12/8; dev.copy2pdf(file="accelRelWgd.pdf", width=2*s, height=2*s, pointsize=8*s)

#' #### Extended Data Figure 9e
#' Compare quartiles of fixed and sample-specific timing approaches:
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

#' xlsx output
#+ ExtendedDataFigure9e
ExtendedDataFigure9e <- createSheet(ExtendedDataFigure9, "ExtendedDataFigure9e")
tab <- data.frame(fixed_acc=t(qWgd[2:4,"5x",]),sample_acc=t(qAccelRelWgd[2:4,]))
addDataFrame(tab, ExtendedDataFigure9e)


#' ### Average rate v time
#' Mutations per year vs time. If the mutation rate was constant, there should be a proportional increase due to the double opportunity to mutate after WGD. 
#+ mutYearTime, fig.height=10, fig.width=10
par(mfrow=c(6,6), mar=c(3,3,2,1),mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
for(n in names(wgdTimeAbsType)){
	a <- age[sample2donor[rownames(wgdTimeAbsType[[n]])]]
	yy <- nDeam22[rownames(wgdTimeAbsType[[n]])]/a
	xx <- 1-t0[rownames(wgdTimeAbsType[[n]])]#y[[n]][,"hat"]/a
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


#' ## MRCA v WGD
#' Compare timing estimates of MRCA and WGD. First, scatter of median
#+ realTimeSubcloneWgdScatter
par( mar=c(4,3,1,1), mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
a <- "5x"
plot(qSubclone["50%",a,dimnames(qWgd)[[3]]], qWgd["50%",a,], col=tissueColors[dimnames(qWgd)[[3]]], pch=16, cex=2, xlab="Median time subclones", ylab="Median time WGD", xlim=c(0,5), ylim=c(0,10))

#' A dotplot may be more informative 
#+ realTimeSubcloneWgd, fig.width=2.5, fig.height=3.5
par( mar=c(3,3,3,10), mgp=c(2,.5,0), tcl=-0.25,cex=1, bty="n", xpd=FALSE, las=1)
w <- "50%"

x <- diag(qSubclone[w,guessAccel[dimnames(qSubclone)[[3]]],])#t[1,3,]
names(x) <- dimnames(qSubclone)[[3]]
y <- diag(qWgd[w,guessAccel[dimnames(qWgd)[[3]]],])#t[1,3,]
names(y) <- dimnames(qWgd)[[3]]

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

#' ## Multi-sample cases
#' Lastly demonstrate concordance of our timing inferences based on multi-sample cases.
donor2sample <- split(names(finalBB),sample2donor[names(finalBB)])
donor2sample <- donor2sample[sapply(donor2sample, length)>1]

#' ### WGD
#' WGD - samples with * removed
#+ wgdMulti, fig.width=6, fig.height=4.5
par(mfrow=c(2,3), mar=c(4,3,2,1), cex=1)
for(n in names(donor2sample)){
	t <- as.character(donor2type[n])
	try({
				s <- wgdTimeAbs[[t]][,,guessAccel[t]]
				s <- s[rownames(s) %in% donor2sample[[n]],,drop=FALSE]
				if(length(s) ==0  ) next 
				if(nrow(s) ==1) next
				plot(s[,"hat"], ylim=c(0, 30), main = paste(n, t), ylab="Time WGD [yr]", xlab="", xaxt="n", xlim=c(0,nrow(s)+1), pch=NA)
				mg14::rotatedLabel(labels=paste0(sample2icgc[rownames(s)],ifelse(rownames(s) %in% remove, "*","")))
				segments(x0=1:nrow(s),y0=s[,'lo'], y1=s[,'up'], col=tissueLines[t],lty=tissueLty[t])
				points(s[,"hat"], pch=21, bg=tissueColors[t], cex=tissueCex[t], col=tissueBorder[t])
				abline(h= mean(s[,'hat']), col=tissueLines[t], lty=3)
			})
}

#' ### MRCA
#' MRCA - samples with * removed
#+ mrcaMulti, fig.width=20, fig.height=20
par(mfrow=c(8,8), mar=c(4,3,2,1), cex=1)
for(n in names(donor2sample)){
	t <- as.character(donor2type[n])
	try({
				s <- subclonesTimeAbs[[t]][,,guessAccel[t]]
				s <- s[rownames(s) %in% donor2sample[[n]],,drop=FALSE]
				if(length(s) ==0  ) next 
				if(nrow(s) ==1) next
				plot(s[,"hat"], ylim=c(0, 20), main = paste(n, t), ylab="Time WGD [yr]", xlab="", xaxt="n", xlim=c(0,nrow(s)+1), pch=NA)
				mg14::rotatedLabel(labels=paste0(sample2icgc[rownames(s)],ifelse(rownames(s) %in% remove, "*","")))
				par(xpd=NA)
				segments(x0=1:nrow(s),y0=s[,'10%'], y1=s[,'90%'], col=tissueLines[t],lty=tissueLty[t])
				points(s[,"hat"], pch=21, bg=tissueColors[t], cex=tissueCex[t], col=tissueBorder[t])
				par(xpd=TRUE)
				abline(h= mean(s[,'hat']), col=tissueLines[t], lty=3)
			})
}



#' # Mutation spectra
#' In this section we analyse changes of mutation spectra across the clonal timing categories.
#' 
#' ## Timed mutation spectra 
#' Caclulate trinucleotide substitution spectra across all samples and all timing categories.
tncTime <- simplify2array(mclapply(finalSnv, function(vcf) table(tncToPyrimidine(vcf), info(vcf)$CLS)))

#' ## Early to late clonal changes
#' ### Distance-based
#' Calculate cosine distance between clonal and subclonal, plot 10 samples with strongest change and 1000 SNVs each
w <- which(apply(mg14:::asum(tncTime[,c(1,2),], 1)>1000, 2, all))
dEarlyLate <- sapply(w, function(i) {x <- tncTime[,,i]; 1-cosineDist(x[,1,drop=FALSE],x[,2,drop=FALSE])})

makeTitle <- function(n){
	d <- sample2donor[n]
	paste0(sample2icgc[n], ", ", donor2type[d], ", ", age[d], "yr")
}

#' ### Figure 4a and Extended Data Figure 6a
#' As examples show the 10 samples with the most striking changes.
#+ sigChangeEarlyLate, fig.width=2.5, fig.height=2
n <- names(tail(sort(dEarlyLate),10))
sigCol <- sapply(c("#2196F3","#212121","#f44336","#BDBDBD","#8BC34A","#FFAB91"), rep,16)
ttl <- function(){
	l <- gsub("\\[|\\]||>.","", dimnames(tncTime)[[1]])
	mtext(text=l, at=b, las=2, cex=0.25, side=1)
	mtext(at=b[8+16*(0:5)], side=1,text=c("C>A","C>G","C>T","T>A","T>C","T>G"), las=1, line=1)
}
for(s in n){
	.par()
	par(mfrow=c(2,1), las=2, mar=c(2,3,2,1))	
	b <- barplot(as.numeric(tncTime[,1,s]), ylab="SNVs", main=makeTitle(s),col=sigCol,border=NA, font.main=1, las=2, cex.main=1)
	text(x=b[1],y=par("usr")[4], labels="early clonal", xpd=TRUE, adj=c(0,1))
	ttl()
	barplot(as.numeric(tncTime[,2,s]), ylab="SNVs", col=sigCol, border=NA,las=2)
	text(x=b[1],y=par("usr")[4], labels="late clonal", xpd=TRUE, adj=c(0,1))
	ttl()
}

#' xlsx output
Figure4a <- createSheet(Figure4, "Figure4a")
addDataFrame(as.data.frame(tncTime[,1:2,names(which(sample2icgc[n]=="SA434348"))]), Figure4a)
ExtendedDataFigure6a <- createSheet(ExtendedDataFigure6, "ExtendedDataFigure6a")
x <- tncTime[,1:2,n]
dimnames(x)[[3]] <- sample2icgc[n]
addDataFrame(as.data.frame(x), ExtendedDataFigure6a)


#' ### Test for spectral changes
#' To formatlise what constitutes a measureable change, calculate an LRT
tncLrt <- function(x){
	l1 <- dmultinom(x=x[,1], prob=x[,1], log=TRUE) + dmultinom(x=x[,2], prob=x[,2], log=TRUE)
	l0 <- dmultinom(x=x[,1], prob=x[,1]+x[,2], log=TRUE) + dmultinom(x=x[,2], prob=x[,1]+x[,2], log=TRUE)
	chisq <- 2*(l1-l0)
	df <- nrow(x)-1
	return(c(chisq=chisq, df=df, p=pchisq(chisq, df=df, lower.tail=FALSE)))
}

w <- which(colSums(tncTime[,1,])>0 & colSums(tncTime[,2,])>0) # Samples with at mutations in both categories
pEarlyLate <- apply(tncTime[,1:2,w],3, tncLrt)
table(p.adjust(pEarlyLate["p",], "bonf")<0.05)

#' Refit some mutational signatures, which have apparently been missed in the time-agnostic signature attribution
sbs <- read.csv("../ref/sigProfiler_SBS_signatures_2018_03_28.csv") # Official table

n <- names(which(sample2icgc=="SA434348"))
e <- nmSolve(tncTime[,,n], P = as.matrix(sbs[,grep("(SBS7[a-d])|SBS38|SBS40", colnames(sbs))]), maxIter=10000, tol=1e-6)
a <- rbind(UV=colSums(e[1:4,]),e[5:6,])[,1:3]
apply(t(t(a)/colSums(a))*100,2,mg14:::roundProp)

n <- names(which(sample2icgc=="SA309724"))
e <- nmSolve(tncTime[,,n], P = as.matrix(sbs[,grep("SBS(2|13|4|5)$", colnames(sbs))]), maxIter=10000)
a <- rbind(APOBEC=colSums(e[c(1,4),]), e[2:3,])
apply(t(t(a)/colSums(a))*100,2,mg14:::roundProp)

n <- names(which(sample2icgc=="SA100712"))
e <- nmSolve(tncTime[,,n], P = as.matrix(sbs[,grep("SBS(1|17.|18|5|40)$", colnames(sbs))]), maxIter=10000)
a <- rbind(`17`=colSums(e[c(3:4),]), e[-(3:4),])
apply(t(t(a)/colSums(a))*100,2,mg14:::roundProp)


#' ## Clonal to subclonal changes
#' ### Distance
#' Calculate cosine distance between clonal and subclonal, plot 10 samples with strongest change and 1000 SNVs each
w <- which(mg14:::asum(tncTime[,1:3,], c(1,2))>1000 & mg14:::asum(tncTime[,4,], 1)>1000)
dClonalSubclonal <- sapply(w, function(i) {x <- tncTime[,,i]; 1-cosineDist(x[,1:3] %*% rep(1,3),x[,4,drop=FALSE])})

#' #### Figure 4b and Extended Data Figure 6b
#+ sigChangeClonalSubclonal, fig.width=2.5, fig.height=2
n <- names(tail(sort(dClonalSubclonal),10))
for(s in n){
	.par()
	par(mfrow=c(2,1), las=2, mar=c(2,3,2,1))	
	b <- barplot(as.numeric(rowSums(tncTime[,1:3,s])), ylab="SNVs", main=makeTitle(s),col=sigCol,border=NA, font.main=1, las=2, cex.main=1)
	text(x=b[1],y=par("usr")[4], labels="clonal", xpd=TRUE, adj=c(0,1))
	ttl()
	barplot(as.numeric(tncTime[,4,s]), ylab="SNVs", col=sigCol, border=NA,las=2)
	text(x=b[1],y=par("usr")[4], labels="subclonal", xpd=TRUE, adj=c(0,1))
	ttl()
}

#' xlsx output
x <- tncTime[,c(1,4),n]
x[,1,] <- mg14:::asum(tncTime[,1:3,n], 2)
dimnames(x)[[3]] <- sample2icgc[n]
dimnames(x)[[2]][1] <- "clonal [early/late/NA]"
Figure4b <- createSheet(Figure4, "Figure4b")
addDataFrame(as.data.frame(x[,1:2,"SA557034"]), Figure4b)
ExtendedDataFigure6b <- createSheet(ExtendedDataFigure6, "ExtendedDataFigure6b")
addDataFrame(as.data.frame(x), ExtendedDataFigure6b)

#' ### Test
#' TNC LRT
w <- which(mg14:::asum(tncTime[,1:3,], c(1,2))>0 & mg14:::asum(tncTime[,4,], 1)>0)
pClonalSubclonal <- apply(tncTime[,,w],3, function(x) tncLrt(x %*% matrix(c(1,1,1,0,0,0,0,1), ncol=2)))
table(p.adjust(pClonalSubclonal["p",], "bonf")<0.05)

#' Mean absolute differences
n <- names(which(p.adjust(pEarlyLate["p",], "bonf")<0.05))
madEarlyLate <- rowSums(abs(t(tncTime[,1,n])/colSums(tncTime[,1,n]) - t(tncTime[,2,n])/colSums(tncTime[,2,n])))/2
n <- names(which(p.adjust(pClonalSubclonal["p",], "bonf")<0.05))
madClonalSubclonal <- rowSums(abs(t(mg14:::asum(tncTime[,1:3,n],2))/colSums(mg14:::asum(tncTime[,1:3,n],2)) - t(tncTime[,4,n])/colSums(tncTime[,4,n])))/2

#' More signature changes
n <- names(which(sample2icgc=="SA557034"))
e <- nmSolve(tncTime[,,n], P = as.matrix(sbs[,grep("SBS(1|2|13|3|5)$", colnames(sbs))]), maxIter=10000)
e <- cbind(clonal=rowSums(e[,1:3]), subclonal=e[,4])
a <- rbind(`APOBEC`=colSums(e[c(2,5),]), e[-c(2,5),])
apply(t(t(a)/colSums(a))*100,2,mg14:::roundProp)

#' ## Signature analysis
#' Now assess how the signature attribution changes in samples with measurable spectral shifts. These data were
#' calculated by Clemency Jolly.
sfc <- read.table("../ref/2019-01-03-allSignatureChanges.txt", header=TRUE, sep="\t")

v <- sfc$samplename %in% names(which(p.adjust(pEarlyLate["p",])<0.05))
#mg14:::ggPlot(pmin(100,pmax(0.01,exp(sfc$log2fc_earlyLate[v]*log(2)))),sfc$signature[v], log='y')

#' ### Early to late clonal change
#' #### Figure 4c
#+ sigFcEarlyLate, fig.width=6, fig.height=3
#cairo_pdf("sigFcEarlyLate.pdf", 6, 3, pointsize=8)
.par()
par(mar=c(9,3,1,1), bty="n")
s <- names(sort(sapply(split(sfc$log2fc_earlyLate[v]*log(2)/log(10), sfc$signature[v]), median, na.rm=TRUE)))
s <- setdiff(s, c(names(which(table(sfc$signature[v])<10)),"n_unassigned"))
t <- donor2type[sample2donor[as.character(sfc$samplename)[v]]]
f <- function(x) pmax(-2,pmin(2,x*log(2)/log(10)))
x <- jitter(as.numeric(factor(as.character(sfc$signature[v]), levels=s)))
plot(f(sfc$log2fc_earlyLate[v]) ~ x, pch=NA, xaxt='n', xlab="", ylab="fold change", ylim=c(-2.2,2.2), xaxt="n", yaxt="n", lwd=0.5)
#segments(x, y0=f(sfc$lCI_earlyLate[v]), y1=f(sfc$uCI_earlyLate[v]), col='grey')
points(f(sfc$log2fc_earlyLate[v]) ~ x, bg=tissueColors[t], pch=21, col=tissueBorder[t], cex=tissueCex[t])
boxplot(sfc$log2fc_earlyLate[v]*log(2)/log(10) ~factor(as.character(sfc$signature[v]), levels=s), pch=NA, staplewex=0, lty=1, add=TRUE, col=NA, xaxt="n", yaxt="n")
mg14::rotatedLabel(x0=seq_along(s), labels=s)
axis(side=2, at=-2:2, labels=paste0(c("\u2265","","","","\u2264"),10^(-2:2)))
m <- log10(2:9 %o% 10^(-2:1))
axis(side=2, at=m, tcl=-0.1, labels=rep("", length(m)))
axis(side=1, at=seq_along(s), labels=rep("", length(lengths(s))))
text(0,2.2, "late", pos=4)
text(0,-2.2, "early", pos=4)
abline(h=0, lty=3)
#dev.off()

#' xlsx output
Figure4c <- createSheet(Figure4, "Figure4c")
addDataFrame(data.frame(sfc[v,c("samplename","signature","log2fc_earlyLate")], tumour_type=t), Figure4c)


#' Some numbers
t(signif(2^sapply(split(sfc$log2fc_earlyLate[v] , sfc$signature[v])[s], quantile, na.rm=TRUE),2))

#' ### Clonal to subclonal changes
#' #### Figure 4d
w <- sfc$samplename %in% names(which(p.adjust(pClonalSubclonal["p",])<0.05))

#+ sigFcClonalSubclonal, fig.width=6, fig.height=3
#cairo_pdf("sigFcClonalSubclonal.pdf", 6, 3, pointsize=8)
.par()
par(mar=c(9,3,1,1), bty="n")
s <- names(sort(sapply(split(sfc$log2fc_clonalSubclonal[w]*log(2)/log(10), sfc$signature[w]), median, na.rm=TRUE)))
s <- setdiff(s, c(names(which(table(sfc$signature[w])<10)),"n_unassigned"))
t <- donor2type[sample2donor[as.character(sfc$samplename)[w]]]
f <- function(x) pmax(-2,pmin(2,x*log(2)/log(10)))
x <- jitter(as.numeric(factor(as.character(sfc$signature[w]), levels=s)))
plot(f(sfc$log2fc_clonalSubclonal[w]) ~ x, pch=NA, xaxt='n', xlab="", ylab="fold change", ylim=c(-2.2,2.2), xaxt="n", yaxt="n", lwd=0.5)
#segments(x, y0=f(sfc$lCI_clonalSubclonal[v]), y1=f(sfc$uCI_clonalSubclonal[v]), col='grey')
points(f(sfc$log2fc_clonalSubclonal[w]) ~ x, bg=tissueColors[t], pch=21, col=tissueBorder[t], cex=tissueCex[t])
boxplot(f(sfc$log2fc_clonalSubclonal[w]) ~factor(as.character(sfc$signature[w]), levels=s), pch=NA, staplewex=0, lty=1, add=TRUE, col=NA, xaxt="n", yaxt="n")
mg14::rotatedLabel(x0=seq_along(s), labels=s)
axis(side=2, at=-2:2, labels=paste0(c("\u2265","","","","\u2264"),10^(-2:2)))
m <- log10(2:9 %o% 10^(-2:1))
axis(side=2, at=m, tcl=-0.1, labels=rep("", length(m)))
axis(side=1, at=seq_along(s), labels=rep("", length(lengths(s))))
text(0,2.2, "subclonal", pos=4)
text(0,-2.2, "clonal", pos=4)
abline(h=0, lty=3)
#dev.off()

#' xlsx output
Figure4d <- createSheet(Figure4, "Figure4d")
addDataFrame(data.frame(sfc[w,c("samplename","signature","log2fc_clonalSubclonal")], tumour_type=t), Figure4d)


#' Some numbers
t(signif(2^sapply(split(sfc$log2fc_clonalSubclonal[w] , sfc$signature[w])[s], quantile, na.rm=TRUE),2))
 
#' Legend for figure
#+ sigFcLegend, fig.width=2, fig.height=4
#pdf("legend.pdf", 2,4)
.par()
par(mar=c(0,0,0,0), bty="n")
plot(NA,NA, xlim=c(-1,1), ylim=c(-1,1), xaxt="n",yaxt="n", xlab='', ylab="")
l <- as.character(sort(unique(donor2type[sample2donor[as.character(sfc$samplename)[w | v]]])))
n <- unique(c(names(which(p.adjust(pClonalSubclonal["p",])<0.05)),names(which(p.adjust(pEarlyLate["p",])<0.05))))
m <- unique(c(names(pClonalSubclonal["p",]),names(pEarlyLate["p",])))
s <- table(donor2type[sample2donor[n]])[l]
t <- table(donor2type[sample2donor[m]])[l]
#l <- l[order(-s/t)]
legend("top", legend=paste0(l, " (",s[l],"/",t[l],")"), pt.bg=tissueColors[l], pch=21, col=tissueBorder[l], pt.cex=tissueCex[l], bty="n")
dev.off()
 


#' # Outputs
#' ## Figure data
#' Save output underlying figures.
saveWorkbook(Figure1,'Figure1.xlsx')
saveWorkbook(Figure2,'Figure2.xlsx')
saveWorkbook(Figure4,'Figure4.xlsx')
saveWorkbook(Figure5,'Figure5.xlsx')
saveWorkbook(ExtendedDataFigure3,'ExtendedDataFigure3.xlsx')
saveWorkbook(ExtendedDataFigure6,'ExtendedDataFigure6.xlsx')
saveWorkbook(ExtendedDataFigure8,'ExtendedDataFigure8.xlsx')
saveWorkbook(ExtendedDataFigure9,'ExtendedDataFigure9.xlsx')

#' ## Real time WGD & MRCA
#' Other working group output, to be released on synapse https://www.synapse.org/#!Synapse:syn14193595.
#' These data may be more useful than the specific subsets underlying each published figure.
#+ wgdMrcaOut
n <- names(finalSnv)
wgdMrcaTimingData <- data.frame(
		uuid=n,
		icgc_sample_id=sample2icgc[n], 
		icgc_donor_id=sample2donor[n], 
		tissue=donor2type[sample2donor[n]], 
		WGD=isWgd, 
		ploidy=finalPloidy,
		eff_ploidy=effGenome,
		purity=finalPurity,
		age=round(age[sample2donor[n]]),
		n_snv_mnv=sapply(finalSnv, nrow),
		CpG_TpG_trunk_pwradj=round(branchDeam[,2],2),
		CpG_TpG_subclonal_branch_pwradj=round(branchDeam[,1],2),
		CpG_TpG_subclonal_linear_pwradj=round(branchDeamLinear[,1],2),
		TiN = TiN[sample2donor[n]],
		remove = n %in% remove,
		accel = ifelse(n %in% remove,NA,guessAccel[as.character(donor2type[sample2donor[n]])]), 
		row.names=n)

t <- as.data.frame(round(Reduce("rbind", wgdTimeAbsType),2))
colnames(t) <- c("time", "time.lo","time.up")
w <- data.frame(row.names=rownames(t), t, `CpG_TpG_total`=nDeam[rownames(t)])
colnames(w) <- sub("lo", "10%", sub("up","90%", colnames(w)))

t <- cbind(branching=as.data.frame(round(Reduce("rbind", subclonesTimeAbsType),2)), linear=as.data.frame(round(Reduce("rbind", subclonesTimeAbsTypeLinear),2)))
colnames(t) <- sub("\\.hat","", colnames(t))

wgdMrcaTimingData <- cbind(wgdMrcaTimingData, 
		WGD=w[rownames(wgdMrcaTimingData),],
		MRCA.time=t[rownames(wgdMrcaTimingData),]
)

write.table(wgdMrcaTimingData, file=paste0(Sys.Date(),"-wgdMrcaTiming.txt"), sep="\t", quote=FALSE, row.names=FALSE)

#' ## All segments, MutationTime.R raw values
#+ segOut, eval=FALSE
t <- do.call("rbind", mclapply(finalBB, function(bb) as.data.frame(bb[,c(1:6,38:47)])))
n <- rownames(t) 
t <- as.data.frame(lapply(t, function(x) if(class(x)=="numeric") round(x,3) else x)) 
t$sample <- sub("\\..+","",n)
write.table(t, file=paste0(Sys.Date(),"-allSegmentsTimeRaw.txt"), quote=FALSE, sep="\t")

#' ## Drivers
#+ driverOut
t <- as.data.frame(finalDriversAnnotated)[c(1,2,6,7,9,10,11,12,13,14,31:44)]
t <- as.data.frame(lapply(t, function(x) if(class(x)=="numeric") round(x,3) else x));
write.table(t, file=paste0(Sys.Date(),"-driversTiming.txt"), quote=FALSE, sep="\t")


#' # Session
#' ## Commit
cat(system("git log -n1", intern=TRUE), sep="\n")
#' ## Objects
l <- ls()
data.frame(variable=l, Reduce("rbind",lapply(l, function(x) data.frame(class=class(get(x)), size=format(object.size(get(x)), units="auto")))))
#' ## Packages
sessionInfo()
devtools::session_info()

#' # Other code
#' All code is available at github.com/gerstung-lab/PCAWG-11
#+ additionalCode, cache=FALSE, echo=FALSE, eval=TRUE
read_chunk('../modules/MutationTime.R/MutationTime.R', labels="MutationTimer")
read_chunk('./PCAWG-functions.R', labels="PCAWG-functions")
read_chunk('./VCF-annotate.R', labels="VCF-annotate")

#' ## MutationTime.R
#' See https://github.com/gerstung-lab/MutationTime.R
#+ MutationTimer, eval=FALSE

#' ## PCAWG-functions.R
#' All basic functions
#+ PCAWG-functions, eval=FALSE

#' ## VCF-annotate.R
#+ VCF-annotate, eval=FALSE
