#' # Pilot 63 data analysis
#' Sanger pipeline only

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


#' ### Libraries
library(VariantAnnotation)
library(Matrix)
library(CoxHD)
library(igraph)

vcfPath <- 'pilot63'
dpPath <- '/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/pilot_64/sanger_analysis/deliverables/dirichlet_clustering_002_calibration_final/sanger/data/'
mutsigDrivers <- read.table('/nfs/team78pc10/im3/Reference_data/putative_cancer_genes_MutSigCV_5000.txt')$V1
purityPloidy <- read.table('/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/pilot_64/sanger_analysis/deliverables/battenberg_002_calibration_final/purity_ploidy.txt', header=TRUE, row.names=1)

sampleIds <- sub("_DPou.+","",dir(dpPath))

loadClusters <- function(ID){
 	file <- paste0(dpPath,"/",ID,"_DPoutput_1250iters_250burnin/",ID,"_optimaInfo.txt")
	read.table(file, header=TRUE, sep="\t")
}

loadVcf <- function(ID){
	file <- dir(vcfPath, pattern=paste0(ID, ".+somatic.snv_mnv.TNC.vcf.bgz$"), full.names=TRUE)
	pos <- loadPositions(ID)
	vcf <- readVcf(file, genome="GRCh37") #, param=ScanVcfParam(which=pos))
	vcf <- vcf[findOverlaps(pos, vcf, select="first")]
	vcf <- addDriver(vcf, mutsigDrivers)
	i = header(vcf)@header$INFO
	exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="Numeric",Description="DP cluster", row.names="DPC"))
	i = header(vcf)@header$INFO
	exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="Numeric",Description="DP cluster probability", row.names="DPP"))
	info(vcf)$DPC <- pos$cluster
	info(vcf)$DPP <- pos$likelihood
	vcf
}

isDeamination <- function(vcf) grepl("(A.CG)|(T..CG)", paste(as.character(unlist(alt(vcf))),vcf@info$TNC))

testNonsyn <- function(vcf) sapply(info(vcf)$VC, function(x) if(length(x) ==0) FALSE  else x %in% c('nonsense','missense','ess_splice'))

addDriver <- function(vcf, mutsigDrivers){
	r <- paste(paste0("^",mutsigDrivers,"(?=\\|)"), collapse="|")
	rowIdx <- grepl(r, info(vcf)$VD, perl=TRUE) & testNonsyn(vcf)
	g <- sapply(info(vcf)$VD, function(x) sub("\\|.+","", x))
	d <- ifelse(rowIdx,g,character(0))
	d[is.na(d)] <- chararc
	i = header(vcf)@header$INFO
	exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="String",Description="Driver gene", row.names="DG"))
	info(vcf)$DG <- as(d,"CharacterList")
	vcf
}

loadAssignment <- function(ID){
	file <- paste0(dpPath,"/",ID,"_DPoutput_1250iters_250burnin/",ID,"_DP_and_cluster_info.txt")
	read.table(file, header=TRUE, sep="\t")
}

addAssignment <- function(vcf, ID){
	a <- loadAssignment(ID)
	i = header(vcf)@header$INFO
	exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=ncol(a),Type="Numeric",Description="DP probability", row.names="DPP"))
	info(vcf)$DPP <- as.matrix(a)
	vcf
}

loadPositions <- function(ID){
	file <- paste0(dpPath,"/",ID,"_DPoutput_1250iters_250burnin/",ID,"_1250iters_250burnin_bestConsensusAssignments.bed")
	tmp <- read.table(file, header=TRUE, sep="\t")
	GRanges(tmp$chr, IRanges(tmp$start+1, tmp$end), cluster=tmp$cluster, likelihood=tmp$likelihood)
}

applyPigeonHole <- function(ID){
	c <- loadClusters(ID)
	p <- purityPloidy[ID,"purity"]
	mcf <- c$location*p
	l <- sapply(1:length(mcf), function(i) mcf[i] > pmax(0,1-mcf))
	which(l & upper.tri(l), arr.ind=TRUE)
}

reduceToCoverRelations <- function(rel){
	if(nrow(rel) ==0) return(rel)
	n <- max(rel)
	P <- matrix(FALSE,n,n)
	for(i in 1:nrow(rel))
		P[rel[i,1],rel[i,2]] <- TRUE
	for(i in 1:n){
		q <- list()
		visit <- logical(n)
		for(j in 1:n)
			if(P[i,j])
				q <- c(q,j)
		while(length(q)>0){
			j <- q[[1]]
			q[[1]] <- NULL
			for(k in 1:n){
				if(P[j,k] & !visit[k]){
					visit[k] <- TRUE
					q <- c(q,k)
					if(P[i,k])
						P[i,k] <- FALSE
					if(P[k,i])
						stop("Error.")
						
				}
			}
		}
	}
	which(P, arr.ind=TRUE)
}

branchLengths <- function(vcf, type=c("all","deamination")){
	type <- match.arg(type)
	if(type=="deamination")
		w <- isDeamination(vcf)
	else
		w <- rep(TRUE, nrow(vcf))
	t <- table(info(vcf)$DPC, w)[,"TRUE"]
	names(t) <- sort(na.omit(unique(info(vcf)$DPC)))
	t
}

allPosets <- sapply(sampleIds, function(ID){
			reduceToCoverRelations(applyPigeonHole(ID))
		})

#' ### Load all to VCF and annotate
#+ allVcf, cache=TRUE
allVcf <- sapply(sampleIds, loadVcf)

toGraph <- function(edgelist, branch.length, edge.labels, node.labels=1:max(edgelist)){
	g <- graph.edgelist(edgelist)
	E(g)$weight <- branch.length
	E(g)$name <- edge.labels
	V(g)$name <- node.labels
	return(g)	
}

allGraphs <- lapply(sampleIds, function(ID){
			edgelist <- reduceToCoverRelations(applyPigeonHole(ID))
			branch.length <- branchLengths(allVcf[[ID]], type='deam')
			b <- as.numeric(names(branch.length))
			edgelist <- edgelist[edgelist[,2] %in% b & edgelist[,1] %in% b,,drop=FALSE]
			edgelist <- rbind(edgelist, cbind(setdiff(b, edgelist[,1]),max(b)+1))
			edgelist <- edgelist[,2:1, drop=FALSE]
			labels <- sapply(split(unlist(info(allVcf[[ID]])$DG), info(allVcf[[ID]])$DPC), function(x) paste(f(x), collapse=","))
			r <- max(edgelist):1
			branch.length <- branch.length[edgelist[,2]]
			labels <- labels[edgelist[,2]]
			edgelist <- cbind(r[edgelist[,1]], r[edgelist[,2]]) ## reverse node ids
			toGraph(edgelist, branch.length, edge.labels=labels, node.labels=r[1:max(edgelist)])
		})
names(allGraphs) <- sampleIds

#' Plot all posets
#+ allGraphs, fig.width=2, fig.height=2
par(mar=c(0,0,0,0))
i <- 1; for(g in allGraphs){
	plot(g, layout=layout.reingold.tilford(g), main=names(allGraphs)[i], edge.label=E(g)$name); i<-i+1}

		
posetDist <- function(graph) {
	e <- get.edgelist(graph)
	w <-  c(0,E(graph)$weight[order(e[,2])]/2)
	d <- shortest.paths(graph, mode='out') - rep(w, each=length(w))
	diag(d) <- NA
	d
}

allDist <- sapply(allGraphs, function(g){ #all
			posetDist(graph = g)
		})

f <- function(x) x[!is.na(x)]
allDrivers <- sapply(allVcf, function(x) {sapply(split(unlist(info(x)$DG), info(x)$DPC), f)})
t <- table(unlist(allDrivers))
sort(t)

allMutations <- names(t)[which(t >=3)]
allMutationsGermline <- c("germline",allMutations)

getMutationCluster <- function(allMutations, vcf){
	m <- match(allMutations, unlist(info(vcf)$DG))
	info(vcf)$DPC[m]
}

allDistMutations <- array(0, c(rep(length(allMutationsGermline),2), length(sampleIds)), dimnames=list(allMutationsGermline, allMutationsGermline, sampleIds))
for(n in sampleIds){
	c <- getMutationCluster(allMutations, allVcf[[n]])
	d <- allDist[[n]]
	m <- c(germline=max(as.numeric(colnames(d))), c)
	l <- match(paste(1:max(as.numeric(colnames(d)))), colnames(d))
	d <- d[l, l]
	#d[is.infinite(d)] <- NA
	allDistMutations[,,n] <- d[m,m]
}

weights <- sapply(sampleIds, function(x){
			m <- getMutationCluster(allMutations, allVcf[[x]])
			l <- branchLengths(allVcf[[x]], type="deam")
			1/l[m]
		})
rownames(weights) <- allMutations

#' # Fit model
#+ fit, results='asis'
observed <- !is.na(allDistMutations) & !is.infinite(allDistMutations) & c(TRUE,rep(FALSE, length(allMutations))) # only distance to root for the moment..
w <- which(observed, arr.ind = TRUE)
y <- allDistMutations[observed]
x <- CoxHD:::MakeInteger(factor(w[,2], levels=seq_along(allMutationsGermline), labels=allMutationsGermline)) - MakeInteger(factor(w[,1], levels=seq_along(allMutationsGermline), labels=allMutationsGermline))

fit <- lm(y ~ . -1 -germline, data=x , weights=1/mapply(function(i,j) weights[i,j], w[,2]-1,w[,3]))
s <- summary(fit)
kable(s$coefficients)

#' # Plot
r <- rank(coef(fit), ties.method = 'random')
c <- coef(fit)
u <- unlist(allDrivers)
par(bty="n", xpd=NA, mar=c(5,3,1,5))
v1 <- s$coefficients[,2]^2
d1 <- colSums(abs(x[,-1]))
d0 <- 0 # prior df
v0 <- 1 # prior variance, i.e., any time
v <- (d1*v1 + d0 * v0)/(d0+d1)  
sd <- sqrt(v)
plot(c, r, xlab="relative time", yaxt='n', pch=19, cex=sqrt(table(u)[gsub("`","",names(r))]/2), xlim=range(c(c+sd, c-sd)))
segments(pmax(0,c-sd),r, pmax(1,c+sd),r)
text(pmax(1,c+sd),r, paste0(names(c), ", n=",table(u)[names(c)]), pos=4)

#' Pairwise precendences
allPrecMutations <- rowSums(!is.na(allDistMutations) & !is.infinite(allDistMutations), dims=2)
Matrix(allPrecMutations)