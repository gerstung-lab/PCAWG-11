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

vcfPath <- '../pilot63'
dpPath <- '/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/pilot_64/sanger_analysis/deliverables/dirichlet_clustering_002_calibration_final/sanger/data/'
mutsigDrivers <- read.table('/nfs/team78pc10/im3/Reference_data/putative_cancer_genes_MutSigCV_5000.txt')$V1
purityPloidy <- read.table('/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/pilot_64/sanger_analysis/deliverables/battenberg_002_calibration_final/purity_ploidy.txt', header=TRUE, row.names=1)
cnPath <- '/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/for_moritz/mle_estimates_14Apr15'
bbPath <- '/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/pilot_64/sanger_analysis/deliverables/battenberg_002_calibration_final/output'

sampleIds <- sub("_DPou.+","",dir(dpPath))

loadClusters <- function(ID){
 	file <- paste0(dpPath,"/",ID,"_DPoutput_1250iters_250burnin/",ID,"_optimaInfo.txt")
	read.table(file, header=TRUE, sep="\t")
}

parseRegion <- function(regions){
	chr <- regmatches(regions, regexpr("^.+(?=:)",regions,perl=TRUE))
	start <- as.numeric(regmatches(regions, regexpr("(?<=:)[0-9]+(?=-)",regions,perl=TRUE)))
	end <- as.numeric(regmatches(regions, regexpr("(?<=-)[0-9]+$",regions,perl=TRUE)))
	data.frame(chr,start,end)
}

loadCn <- function(ID){
	file <- paste0(cnPath, "/",ID,"_timingsMle.txt")
	tab <- read.table(file, header=TRUE, sep=" ")
	reg <- parseRegion(tab$ascatId)
	GRanges(tab$chr, IRanges(reg$start,reg$end), strand="*", tab[-1])
}

loadBB <- function(ID){
	file <- paste0(bbPath, "/",ID,"_segments.txt")
	tab <- read.table(file, header=TRUE, sep='\t')
	GRanges(tab$chromosome, IRanges(tab$start, tab$end), strand="*", tab[-3:-1])
}

allBB <- sapply(sampleIds, loadBB) 

getTumorCounts <- function(vcf){
	sapply(grep("(F|R).Z", names(geno(vcf)), value=TRUE), function(i) geno(vcf)[[i]][,"TUMOUR"])
}

getTumorDepth <- function(vcf)
	rowSums(getTumorCounts(vcf))

getAltCount <- function(vcf){
	c <- getTumorCounts(vcf)
	t <- c[,1:4] + c[,5:8]
	colnames(t) <- substring(colnames(t),2,2)
	a <- as.character(unlist(alt(vcf)))
	sapply(seq_along(a), function(i) t[i, a[i]])
}

computeMutCn <- function(vcf, bb){
	f <- findOverlaps(vcf, bb)
	cn <- split(bb$copy_number[subjectHits(f)], queryHits(f))
	cf <- split(bb$clonal_frequency[subjectHits(f)], queryHits(f))
	p <- purityPloidy[ID, 'purity']
	a <- getAltCount(vcf)
	d <- getTumorDepth(vcf)
	d <- t(sapply(seq_along(a), function(i){
						if(!i %in% names(cn)) return(rep(NA,4))
						cni <- cn[[as.character(i)]]
						cfi <- cf[[as.character(i)]]
						totalCnTumor <- sum(cni * cfi)
						totalCnNormal <- 2
						cnStates <- unlist(sapply(seq_along(cni), 
										function(j){s <- 0:cni[j] * cfi[j]
											names(s) <- paste0(j,':',0:cni[j]) 
											return(s)}, simplify=FALSE)) * p / (totalCnTumor * p + (1-p) * totalCnNormal)
						prob <- sapply(cnStates, function(s) dbinom(a[i],d[i],s)) + .Machine$double.eps
						prob <- prob/sum(prob)
						w <- which.max(prob)
						idx <- as.numeric(strsplit(names(prob[w]), ":")[[1]])
						c(MCN=idx[2], TCN=cni[idx[1]], CNF=cfi[idx[1]], PMCN=prob[w])
					}))
	d <- DataFrame(d, CNID = as(f,"List"))
	colnames(d) <- gsub("[[:punct:]]+[0-9]+","", colnames(d))
	return(d)
}

addMutCn <- function(vcf, bb){
	i = header(vcf)@header$INFO
	exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=c(1,1,1,1,"."),Type=rep("Numeric",5), Description=c("Mutation copy number","Total copy number","Copy number frequency (relative to all cancer cells)", "MCN probability","BB segment ids"), row.names=c("MCN","TCN","CNF","PMCN","CNID")))
	info(vcf) <- cbind(info(vcf), computeMutCn(vcf, bb))
	return(vcf)
}

t <- read.table("ref/zack2013_140.txt", sep="\t", header=TRUE)[-141,]
r <- parseRegion(t$Peak.region)
gisticCn <- GRanges(sub("chr","",r$chr), IRanges(r$start,r$end), strand="*",t)

tryExceptNull <- function(x) if(class(x)=="try-error") GRanges() else x

allCn <- do.call("GRangesList",sapply(sampleIds, function(x) tryExceptNull(try(loadCn(x)))))

c <- coverage(allCn)
c[c >= 3]

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
	#d[is.na(d)] <- character(0)
	i = header(vcf)@header$INFO
	exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="String",Description="Driver gene", row.names="DG"))
	info(vcf)$DG <- as(d,"CharacterList")
	vcf
}

addMcn <- function(vcf, ID){
	
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

tncToPyrimidine <- function(vcf){
	a <- unlist(alt(vcf))
	r <- ref(vcf)
	tnc <- DNAStringSet(info(vcf)$TNC)
	rc <- grepl("A|G", r)
	tnc[rc] <- reverseComplement(tnc[rc])
	a[rc] <- reverseComplement(a[rc])
	t <- paste0(substr(tnc,1,1), "[",substr(tnc,2,2), ">",a, "]", substr(tnc,3,3))
	n <- c("A","C","G","T")
	f <- paste0(rep(n, each=4), "[", rep(c("C","T"), each=96/2), ">", c(rep(c("A","G","T"), each=48/3), rep(c("A","C","G"), each=48/3)), "]", n)
	factor(t, levels=f)
} 

applyPigeonHole <- function(ID){
	c <- loadClusters(ID)
	p <- purityPloidy[ID,"purity"]
	mcf <- c$location#*p
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

#' ### Add mutation copy numbers
#+ allVcfMCN, cache=TRUE
allVcf <- sapply(sampleIds, function(ID) addMutCn(allVcf[[ID]], allBB[[ID]]))

#' ### Compute graphs (posets)
toGraph <- function(edgelist, branch.length, edge.labels, node.labels=1:max(edgelist)){
	g <- graph.edgelist(edgelist)
	E(g)$weight <- branch.length
	E(g)$name <- edge.labels
	V(g)$name <- node.labels
	return(g)	
}

na.rm <- function(x) x[!is.na(x)]

allGraphs <- lapply(sampleIds, function(ID){
			edgelist <- reduceToCoverRelations(applyPigeonHole(ID))
			branch.length <- branchLengths(allVcf[[ID]], type='deam')
			b <- as.numeric(names(branch.length))
			edgelist <- edgelist[edgelist[,2] %in% b & edgelist[,1] %in% b,,drop=FALSE]
			edgelist <- rbind(edgelist, cbind(setdiff(b, edgelist[,1]),max(b)+1))
			edgelist <- edgelist[,2:1, drop=FALSE]
			labels <- sapply(split(unlist(info(allVcf[[ID]])$DG), info(allVcf[[ID]])$DPC), function(x) paste(na.rm(x), collapse=","))
			r <- max(edgelist):1
			branch.length <- branch.length[edgelist[,2]]
			labels <- labels[edgelist[,2]]
			edgelist <- cbind(r[edgelist[,1]], r[edgelist[,2]]) ## reverse node ids
			node.size <- branchLengths(allVcf[[ID]])
			alleleFrequency <- loadClusters(ID)$location
			node.labels <- c("Germline", paste0(r[1:max(edgelist)][-1], ", n=", node.size[r[-1]], ", f=", round(alleleFrequency[r[-1]],2)))
			g <- toGraph(edgelist, branch.length, edge.labels=labels, node.labels=node.labels)
			V(g)$size <- c(1-purityPloidy[ID,"purity"],alleleFrequency[r[-1]])
			return(g)
		})
names(allGraphs) <- sampleIds

plotPoset <- function(g){
	c <- colorRampPalette(brewer.pal(9, "Spectral"))(10)
	plot(g, layout=layout.reingold.tilford(g), edge.label=E(g)$name, vertex.size = 25*pmin(1,V(g)$size), edge.width=E(g)$weight/100)
}

#' Plot all posets
#+ allGraphs, fig.width=2, fig.height=2
par(mar=c(1,1,1,1))
i <- 1; for(g in allGraphs){
	plotPoset(g);
	title(main=names(allGraphs)[i])
	i<-i+1
}

#' Distance in poset		
posetDist <- function(graph) {
	e <- get.edgelist(graph)
	w <-  c(0,E(graph)$weight)
	names(w) <- c("Germline", e[,2])
	d <- shortest.paths(graph, mode='out')
	d <- d - rep(w[colnames(d)], each=length(w))/2
	diag(d) <- NA
	d
}

allDist <- sapply(allGraphs, function(g){ #all
			posetDist(graph = g)
		})

#' Map to drivers
allDrivers <- sapply(allVcf, function(x) {sapply(split(unlist(info(x)$DG), info(x)$DPC), na.rm)})
t <- table(unlist(allDrivers))
sort(t)

minDriver <- 2
allMutations <- names(t)[which(t >= minDriver)]
allMutations <- allMutations[allMutations!="ARID1A"]
allMutationsGermline <- c("germline",allMutations)

getMutationCluster <- function(allMutations, vcf){
	m <- match(allMutations, unlist(info(vcf)$DG))
	info(vcf)$DPC[m]
}

#' Combine into array
allDistMutations <- array(0, c(rep(length(allMutationsGermline),2), length(sampleIds)), dimnames=list(allMutationsGermline, allMutationsGermline, sampleIds))
for(n in sampleIds){
	c <- getMutationCluster(allMutations, allVcf[[n]])
	d <- allDist[[n]]
	colnames(d) <- sub(",.+","", colnames(d))
	colnames(d) <- sub("Germline", max(as.numeric(sub("Germline","0",colnames(d))))+1, colnames(d))
	m <- c(germline=max(as.numeric(colnames(d))), c)
	l <- match(paste(1:max(as.numeric(colnames(d)))), colnames(d))
	d <- d[l, l]
	#d[is.infinite(d)] <- NA
	allDistMutations[,,n] <- d[m,m]
}

weights <- sapply(sampleIds, function(ID){
			m <- getMutationCluster(allMutations, allVcf[[ID]])
			l <- branchLengths(allVcf[[ID]], type="deam")
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
#+ btTiming, fig.width=5, fig.height=3
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
plot(c, r, xlab="Time [C>T@CpG]", yaxt='n', pch=19, cex=sqrt(table(u)[gsub("`","",names(r))]/2), xlim=pmax(0,range(c(pmin(5000,c+sd), c-sd))), ylab="")
segments(pmax(0,c-sd),r, pmin(5000,c+sd),r)
text(pmax(1,c+sd),r, paste0(names(c), ", n=",table(u)[names(c)]), pos=4)


#' Distributions
boxplot(t(allDistMutations[1,-1,][order(c),]), las=2, horizontal=TRUE)
#' Pairwise precendences
allPrecMutations <- rowSums(!is.na(allDistMutations) & !is.infinite(allDistMutations), dims=2)
Matrix(allPrecMutations)

#' NMF
snmf <- new.env()
load("ref/snmf.RData", envir=snmf)
predictRealTime <- function(x, signatures=snmf$snmfFit$P[-1,]){
	snmf$snmfFit$P[1,1] / snmf$weight * snmf$nmSolve(x, signatures, maxIter=100)[1,]
}

set.seed(42)
sigData <- sapply(allVcf, function(x) table(tncToPyrimidine(x)))
n <- 6
nmfFit <- mynmf(sigData,n)
f <- factor(substr(rownames(d),3,5))
b <- barplot(t(nmfFit$P), las=2, col=brewer.pal(8,"Set1"), border=NA,names.arg=rep("", nrow(sigData)))
mtext(at=b,rownames(sigData), side=1, las=2, cex=.7, col=brewer.pal(6,"Set1")[f])

cor(snmf$snmfFit$P[-1,], nmfFit$P[,])
newSig <- nmfFit$P
w <- which.max(cor(snmf$snmfFit$P[-1,1], nmfFit$P[,]))
newSig <- newSig[,c(w, setdiff(1:n,w))]

cor(predictRealTime(d, newSig),age)

summary(predictRealTime(sigData, newSig))

#' WGD
wgd <- sapply(sampleIds, function(ID){
			if(purityPloidy[ID,2]<3)
				return(rep(NA,2))
			vcf <- allVcf[[ID]]
			t <- table(info(vcf)$MCN, info(vcf)$TCN, isDeamination(vcf))
			c(t["2","4",2], rowSums(t["2",,2,drop=FALSE]))
		})

#+ wdgTime, 3,3
o <- order(colMeans(wgd))
plot(wgd[,o[1:25]], rep(1:25,each=2), xlab="Time [C>T@CpG]", ylab="", pch=rep(c(1,19), 25))
segments(wgd[1,o[1:25]], 1:25, wgd[2,o[1:25]], 1:25)
legend("bottomright", c("MCN=2; CN=4", "MCN=2"), pch=c(1,19))

#' Fix allCN
allCn <- GRangesList(sapply(allCn, function(x){
			reg <- parseRegion(x$ascatId)
			GRanges(reg$chr, ranges(x), strand(x), mcols(x))
		}))

#' Driver CN to obs
driverCnTime <- sapply(allCn, function(cn){
			gamp <- gisticCn[gisticCn$Type=="Amp"]
			gdel <- gisticCn[gisticCn$Type=="Del"]
			c(cn[grep("ains",cn$type)]$pi0[findOverlaps(gamp, cn[grep("ains",cn$type)], select='first')],
			cn[grep("CNLOH",cn$type)]$pi0[findOverlaps(gdel, cn[grep("CNLOH",cn$type)], select='first')])
		})
rownames(driverCnTime) <- gisticCn$Peak.Name
#+ driverCnTime, 4,4
o <- order(rowMeans(driverCnTime, na.rm=TRUE))
w <- rowSums(!is.na(driverCnTime)) >= 10
boxplot(t(driverCnTime[o,][w[o],]), horizontal=TRUE, las=2, col=gisticCn$Type[o][w[o]], xlab="pi0")

#+ driverCnHist, fig.width=6, fig.height=3
barplot(sort(rowSums(!is.na(driverCnTime)), decreasing=TRUE), las=2, border=NA, cex.names=.5)

allClusters <- sapply(sampleIds, loadClusters)

driverStatus <- sapply(sampleIds, function(ID){
			m <- match(mutsigDrivers,unlist(info(allVcf[[ID]])$DG))
			allClusters["location",ID][[1]][info(allVcf[[ID]])$DPC[m]] > 0.95
		})

driverStatus <- apply(driverStatus,1, function(x) table(factor(x, levels=c(TRUE,FALSE), labels=c("clonal","subclonal"))))

#+ driverSubHist, fig.width=6, fig.height=3
w <- colSums(driverStatus)>0
o <- order(colSums(driverStatus[,w]), decreasing=TRUE)
barplot(driverStatus[,w][,o], las=2, legend=TRUE, border=NA, col=brewer.pal(3,"Set1"), cex.names=.5)
dev.off()