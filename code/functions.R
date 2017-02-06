# TODO: Add comment
# 
# Author: mg14
###############################################################################

library(Rsamtools)

vcfPath <- '/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/final_consensus_12oct_passonly/snv_mnv'
basePath <-  '/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/dp/20161213_vanloo_wedge_consSNV_prelimConsCNAallStar'
dpPath <- paste0(basePath,'/2_subclones/')
CANCERGENES <- read.table('/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/ref/cancer_genes.txt')$V1
purityPloidy <- read.table( '/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/consensus.20170119.purity.ploidy.annotated.txt', header=TRUE, row.names=1)
#colnames(purityPloidy) <- c("purity","ploidy")
cnPath <- paste0(basePath,'/4_copynumber/')
bbPath <- paste0(basePath,'/4_copynumber/')

allGender <- read.table('/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/gender/2016_12_09_inferred_sex_all_samples.txt', header=TRUE, sep='\t')
allGender <- allGender[!duplicated(allGender$tumourid) & allGender$tumourid != 'tumourid',]
#write.table(gender,'/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/gender/2016_12_09_inferred_sex_all_samples_CORRECTED_MG.txt', row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
rownames(allGender) <- allGender$tumourid

addTNC <- function(vcf){	
	r = "/lustre/scratch112/sanger/cgppipe/PanCancerReference/genome.fa.gz" #meta(header(v))["reference",]
	if(!"TNC" %in% rownames(header(vcf)@header$INFO)){
		tnc=scanFa(file=r, resize(granges(vcf), 3,fix="center"))
		i = header(vcf)@header$INFO
		exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="String",Description="Trinucleotide context", row.names="TNC"))
		info(vcf)$TNC <- as.character(tnc)
	}
	return(vcf)
}

dpFiles <- dir(dpPath, pattern="_subclonal_structure.txt", recursive=TRUE)
	
bbFiles <- dir(bbPath, pattern="_segments.txt", recursive=TRUE)

loadClusters <- function(ID){
	file <- paste0(dpPath,"/",grep(paste0(ID,"[[:punct:]]"), dpFiles, value=TRUE, perl=TRUE))
	if(grepl(".gz", file))
		file <- gzfile(file)
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
	t <- try({
				file <- grep(paste0(ID,"[[:punct:]]"), dir(bbPath, pattern="segments.txt", recursive=TRUE, full.names=TRUE), value=TRUE)
				if(grepl(".gz", file))
					file <- gzfile(file)
				tab <- read.table(file, header=TRUE, sep='\t')
				GRanges(tab$chromosome, IRanges(tab$start, tab$end), strand="*", tab[-3:-1])
			})
	if(class(t)=='try-error') GRanges(copy_number=numeric(), major_cn=numeric(), minor_cn=numeric(), clonal_frequency=numeric()) else t
}

loadConsensusCNA <- function(ID, purity=1, path="/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/mg14/final/consensus.20170119.somatic.cna.annotated"){
	file <- grep(paste0(ID,"[[:punct:]]"), dir(path, pattern="cna.annotated.txt", recursive=TRUE, full.names=TRUE), value=TRUE)
	if(grepl(".gz", file))
		file <- gzfile(file)
	tab <- read.table(file, header=TRUE, sep='\t')
	subclonalIndex <- !is.na(tab$total_cn) & !is.na(tab$battenberg_nMaj2_A) & !is.na(tab$battenberg_nMin2_A) & !is.na(tab$battenberg_frac2_A) & (tab$battenberg_nMaj1_A == tab$major_cn & tab$battenberg_nMin1_A == tab$minor_cn | tab$battenberg_nMaj2_A == tab$major_cn & tab$battenberg_nMin2_A == tab$minor_cn)
	ix <- c(1:nrow(tab), which(subclonalIndex))
	gr <- GRanges(tab$chromosome, IRanges(tab$start, tab$end), strand="*", clonal_frequency=purity, tab[-3:-1])[ix]
	if(any(subclonalIndex)){
		gr$clonal_frequency[which(subclonalIndex)] <- tab$battenberg_frac1_A[subclonalIndex] * purity
		gr$major_cn[which(subclonalIndex)] <- tab$battenberg_nMaj1_A[subclonalIndex]
		gr$minor_cn[which(subclonalIndex)] <- tab$battenberg_nMin1_A[subclonalIndex]
		gr$clonal_frequency[-(1:nrow(tab))] <- tab$battenberg_frac2_A[subclonalIndex] * purity
		gr$major_cn[-(1:nrow(tab))] <- tab$battenberg_nMaj2_A[subclonalIndex]
		gr$minor_cn[-(1:nrow(tab))] <- tab$battenberg_nMin2_A[subclonalIndex]
	}
	sort(gr)
}

refFile = "/lustre/scratch112/sanger/cgppipe/PanCancerReference/genome.fa.gz" #meta(header(v))["reference",]
refLengths <- scanFaIndex(file=refFile)
chrOffset <- cumsum(c(0,as.numeric(width(refLengths))))
names(chrOffset) <- c(seqlevels(refLengths), "NA")

averagePloidy <- function(bb) {
	sum(width(bb) * bb$copy_number * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
}

averageEvenPloidy <- function(bb){
	sum(width(bb) * (!bb$major_cn %% 2 & !bb$minor_cn %% 2) * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
}

getTumorCounts <- function(vcf){
	sapply(grep("(F|R).Z", names(geno(vcf)), value=TRUE), function(i) geno(vcf)[[i]][,"TUMOUR"])
}

getTumorDepth <- function(vcf){
	if("t_alt_count" %in% colnames(info(vcf))){ ## consensus data, snv and indel
		info(vcf)$t_alt_count + info(vcf)$t_ref_count
	}else{ ## older data
		if("FAZ" %in% rownames(geno(header(vcf)))){
			rowSums(getTumorCounts(vcf))
		}else{
			geno(vcf)$DEP[,2]
		}
	}
}

getAltCount <- function(vcf){
	if("t_alt_count" %in% colnames(info(vcf))){ ## consensus data, snv and indel
		info(vcf)$t_alt_count
	}else{ ## older formats
		if("FAZ" %in% rownames(geno(header(vcf)))){ ## ie subs
			c <- getTumorCounts(vcf)
			t <- c[,1:4] + c[,5:8]
			colnames(t) <- substring(colnames(t),2,2)
			a <- as.character(unlist(alt(vcf)))
			a[!a%in%c('A','T','C','G')] <- NA
			sapply(seq_along(a), function(i) if(is.na(a[i])) NA else t[i, a[i]])
		}
		else{ ## ie indel
			#(geno(vcf)$PP + geno(vcf)$NP + geno(vcf)$PB + geno(vcf)$NB)[,"TUMOUR"]
			geno(vcf)$MTR[,2]
		}
	}
}



mergeClusters <- function(clusters, deltaFreq=0.05){
	if(nrow(clusters) <= 1) return(clusters)
	h <- hclust(dist(clusters$proportion))
	ct <- cutree(h, h=deltaFreq)
	cp <- as.matrix(cophenetic(h))
	Reduce("rbind",lapply(unique(ct), function(i) {
						n_ssms <- sum(clusters$n_ssms[ct==i])
						w <- max(cp[ct==i,ct==i])
						data.frame(new.cluster=i, n_ssms=n_ssms, proportion=sum(clusters$proportion[ct==i]*clusters$n_ssms[ct==i])/n_ssms, width=w)
					}
	))
}

clustersFromBB <- function(bb){
	w <- bb$clonal_frequency == max(bb$clonal_frequency, na.rm=TRUE) | bb$clonal_frequency < 0.5 * max(bb$clonal_frequency, na.rm=TRUE)
	t <- table(bb$clonal_frequency[w])
	cl <- data.frame(cluster=seq_along(t), n_ssms=as.numeric(t), proportion=as.numeric(names(t)))
	mergeClusters(cl, deltaFreq=0.2)
}



probGenotype <- function(vcf){
	dg <- factor(paste(unlist(info(vcf)$DG)), levels=c("NA",as.character(CANCERGENES)))
	P <- pGainToTime(vcf)
	G <- matrix(0, ncol=5, nrow=nlevels(dg), dimnames=list(levels(dg), c(colnames(P),"NA")))
	t <- table(dg)
	for(g in names(t[t>0]))
		G[g,] <- c(colSums(P[dg==g,,drop=FALSE],na.rm=TRUE), "NA"=sum(is.na(P[dg==g,1])))
	return(G)
}

getGenotype <- function(vcf, reclassify='missing', ...){
	cls <- classifyMutations(vcf = vcf, reclassify=reclassify)
	t <- info(vcf)$TCN
	if(is.null(t))
		t <- info(vcf)$MinCN + info(vcf)$MajCN
	hom <- factor(info(vcf)$MutCN==t, levels=c(TRUE,FALSE))
	dg <- factor(unlist(info(vcf)$DG), levels=as.character(CANCERGENES))
	table(gene=dg, class=cls, homozygous=hom, ...)
}

tryExceptNull <- function(x) if(class(x)=="try-error") GRanges() else x

loadVcf <- function(ID){
	file <- dir(vcfPath, pattern=paste0(ID, ".+somatic.snv_mnv.TNC.vcf.bgz$"), full.names=TRUE)
	pos <- loadPositions(ID)
	vcf <- readVcf(file, genome="GRCh37") #, param=ScanVcfParam(which=pos))
	f <- findOverlaps(pos, vcf, select="first")
	vcf <- vcf[na.omit(f)]
	vcf <- addDriver(vcf, CANCERGENES)
	i = header(vcf)@header$INFO
	exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="Numeric",Description="DP cluster", row.names="DPC"))
	i = header(vcf)@header$INFO
	exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="Numeric",Description="DP cluster probability", row.names="DPP"))
	info(vcf)$DPC <- pos$cluster[!is.na(f)]
	info(vcf)$DPP <- pos$likelihood[!is.na(f)]	
	vcf
}

isDeamination <- function(vcf) grepl("(A.CG)|(T..CG)", paste(as.character(unlist(alt(vcf))),vcf@info$TNC))

testDriver <- function(vcf) sapply(info(vcf)$VC, function(x) if(length(x) ==0) FALSE  else any(x %in% c('nonsense','missense','ess_splice','frameshift','inframe','cds_distrupted')))

addDriver <- function(vcf, mutsigDrivers){
	i = header(vcf)@header$INFO
	exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="String",Description="Driver gene", row.names="DG"))
	if(nrow(vcf)==0){
		info(vcf)$DG <- CharacterList()
		return(vcf)
	}
	
	r <- paste(paste0("^",mutsigDrivers,"(?=\\|)"), collapse="|")
	rowIdx <- grepl(r, info(vcf)$VD, perl=TRUE) & testDriver(vcf)
	g <- sapply(info(vcf)$VD, function(x) sub("\\|.+","", x))
	d <- ifelse(rowIdx,g,character(0))
	#d[is.na(d)] <- character(0)
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
	file <- gzfile(paste0(dpPath,"/",ID,"_mutation_assignments.txt.gz"))
	tmp <- read.table(file, header=TRUE, sep="\t")
	#GRanges(tmp$chr, IRanges(tmp$start+1, tmp$end), cluster=tmp$cluster, likelihood=tmp$likelihood)
	GRanges(tmp$chr, IRanges(tmp$pos, width=1), cluster=tmp$cluster)
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
	mcf <- c$proportion#*p
	l <- sapply(1:length(mcf), function(i) mcf[i] > pmax(0,1-mcf))
	w <- which(l & upper.tri(l), arr.ind=TRUE)
	cbind(c$cluster[w[,1]], c$cluster[w[,2]])
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

cnWeights <- function(vcf){
	t <- if(is.null(info(vcf)$TCN)) (info(vcf)$MajCN + info(vcf)$MinCN) else info(vcf)$TCN
	info(vcf)$MutCN / t
}

branchLengths <- function(vcf, type=c("all","deamination")){
	type <- match.arg(type)
	if(type=="deamination")
		w <- isDeamination(vcf)
	else
		w <- rep(TRUE, nrow(vcf))
	cnw <- cnWeights(vcf)
	u <- allClusters[[meta(header(vcf))["ID",1]]]$cluster
	b <- sapply(u, function(c) sum(2*cnw * (info(vcf)$DPC==c & w), na.rm=TRUE))
	#if(length(b)==0) b <- rep(0, length(u))
	names(b) <- u
	return(b)
}

avgWeights <- function(vcf, type=c("all","deamination")){
	type <- match.arg(type)
	if(type=="deamination")
		w <- isDeamination(vcf)
	else
		w <- rep(TRUE, nrow(vcf))
	cnw <- cnWeights(vcf)
	mean(2*cnw[w], na.rm=TRUE)
}

predictRealTime <- function(x, signatures=snmf$snmfFit$P[-1,]){
	snmf$snmfFit$P[1,1] / snmf$weight * snmf$nmSolve(x, signatures, maxIter=100)[1,]
}

#' ### Compute graphs (posets)
toGraph <- function(edgelist, branch.length, edge.labels, node.labels=1:max(edgelist)){
	g <- graph.edgelist(edgelist)
	E(g)$weight <- branch.length
	E(g)$name <- edge.labels
	V(g)$name <- node.labels
	return(g)	
}

na.rm <- function(x) x[!is.na(x)]

plotPoset <- function(g){
	c <- colorRampPalette(brewer.pal(9, "Spectral"))(10)
	plot(g, layout=layout.reingold.tilford(g), edge.label=E(g)$name, vertex.size = 25*pmin(1,V(g)$size), edge.width=E(g)$weight/100)
}

#' Distance in poset		
posetDist <- function(g) {
	e <- get.edgelist(g)
	w <-  c(0,E(g)$weight)
	names(w) <- c("Germline", e[,2])
	d <- shortest.paths(g, mode='out')
	d <- d - rep(w[colnames(d)], each=length(w))/2
	diag(d) <- NA
	d
}

getMutationCluster <- function(allMutations, vcf){
	m <- match(allMutations, unlist(info(vcf)$DG))
	info(vcf)$DPC[m]
}

distAsRange <- function(g){
	e <- get.edgelist(g)
	w <-  c(0,E(g)$weight)
	names(w) <- c("Germline", e[,2])
	d <- shortest.paths(g, mode='out')
	y <- shift(IRanges(-w[colnames(d)],0), d["Germline", ])
	names(y) <- paste0(colnames(d), ", genes=",E(g)$name[match(colnames(d), e[,2])])
	y
}


nmSolve <- function(D, P, maxIter = 500, tol=1e-3) {
	n <- nrow(D)
	m <- ncol(D)
	s <- ncol(P)
	tP <- t(P)
	rP <- rep(colSums(P), m)
	D <- as.matrix(D)
	E1 <- E2 <- matrix(runif(s * m, 1e-3, 1), ncol = m)
	err <- 2*tol
	
	iter <- 1
	while (iter < maxIter & err > tol) {
		E1 <- E2
		E2 <- E1 * (tP %*% (D/(P %*% (E1 + .Machine$double.eps))))/rP
		iter <- iter + 1
		err <- mean(abs(E2 - E1)/(E1+.Machine$double.eps), na.rm=TRUE)
		if(iter %% 100 == 0) cat(round(-log10(err)))
	}
	cat("\n")
	if(iter == maxIter) warning(paste("No convergence after",iter, "iterations."))
	E2
}

wnmSolve <- function(D, P, weights =  rep(0, ncol(P)), maxIter = 500, tol=1e-3) {
	D <- as.matrix(D)
	D <- rbind(D, matrix(weights, ncol=ncol(D), nrow=ncol(P)))
	P <- rbind(P, diag(rep(1,ncol(P))))
	n <- nrow(D)
	m <- ncol(D)
	s <- ncol(P)
	tP <- t(P)
	rP <- rep(colSums(P), m)
	E1 <- E2 <- matrix(runif(s * m, 1e-3, 1), ncol = m)
	err <- 2*tol
	
	iter <- 1
	while (iter < maxIter & err > tol) {
		E1 <- E2
		E2 <- E1 * (tP %*% (D/(P %*% (E1 + .Machine$double.eps))))/rP
		iter <- iter + 1
		err <- mean(abs(E2 - E1)/(E1+.Machine$double.eps), na.rm=TRUE)
		if(iter %% 100 == 0) cat(round(-log10(err)))
	}
	cat("\n")
	if(iter == maxIter) warning(paste("No convergence after",iter, "iterations."))
	E2
}



wgdTest <- function(vcf){
	id <- meta(header(vcf))["ID",1]
	bb <- allBB[[id]]
	ix <- which(bb$copy_number==4 & bb$minor_cn==2)
	v <- vcf[vcf %over% bb[ix]]
	#w <- sum(as.numeric(width(reduce(bb[ix]))))
	t <- table(info(v)$MutCN, info(v)$TCN, as.character(seqnames(v)), info(v)$DPC)
}

#' Power
power <- function(f,n, theta=6.3, err=1e-4) if(any(is.na(c(f,n,theta,err)))) NA else sum((log10(dbinom(0:n, n, 0:n/n) / dbinom(0:n,n,err)) > theta) * dbinom(0:n,n,f))

testIndel <- function(vcf) sapply(info(vcf)$VC, function(x) if(length(x) ==0) FALSE  else any(x %in% c('frameshift','inframe','ess_splice','SO:0001577:complex_change_in_transcript', 'SO:0001582:initiator_codon_change', 'splice_region')))

asum <- function(x, dim) apply(x, setdiff(seq_along(dim(x)), dim), sum)

#' official driver file
#library(VariantAnnotation)
#drivers <- read.table("/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/driver/pcawg_whitelist_coding_drivers_v1_sep302016.txt", header=TRUE, sep="\t")
finalData <- read.table("/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/ref/release_may2016.v1.4.tsv", header=TRUE, sep="\t")
#r <- DNAStringSet(drivers$ref)
#a <- DNAStringSet(drivers$alt)
#m <- sapply(levels(drivers$sample_id), function(x) grep(x, finalData$sanger_variant_calling_file_name_prefix))
#driVers <- VRanges(seqnames = drivers$chr, ranges=IRanges(drivers$pos, width =  width(r)), ref=r, alt=a, sampleNames = finalData$icgc_donor_id[m[drivers$sample_id]])
#mcols(driVers) <- cbind(samples=finalData$sanger_variant_calling_file_name_prefix[m[drivers$sample_id]], drivers[,-c(1,3,4,5,6)])
#save(driVers, file = "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/driver/pcawg_whitelist_coding_drivers_v1_sep302016.RData")
load("/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/driver/pcawg_whitelist_coding_drivers_v1_sep302016.RData")
CANCERGENES <- levels(driVers$gene)

addFinalDriver <- function(vcf, driVers){
	i = header(vcf)@header$INFO
	exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="String",Description="Driver gene", row.names="DG"))
	info(vcf)$DG <- factor(rep(NA,nrow(vcf)), levels = levels(driVers$gene))
	if(nrow(vcf)==0)
		return(vcf)
	ID <- meta(header(vcf))["ID",1]
	d <- driVers[grep(ID, driVers$samples)]
	if(length(d)==0)
		return(vcf)
	overlaps <- findOverlaps(vcf, d)
	g <- factor(rep(NA,nrow(vcf)), levels = levels(d$gene))
	g[queryHits(overlaps)] <- d$gene[subjectHits(overlaps)]
	info(vcf)$DG <- g
	return(vcf)
}


t <- read.table("../ref/PCAWG_colors.tsv", sep='\t', header=FALSE, comment="")
tissueColors <- as.character(t$V4)
names(tissueColors) <- t$V2
clinicalData <- read.table("../ref/pcawg_donor_clinical_August2016_v7-2.tsv", header=TRUE, sep="\t", comment="", quote="")
specimenData <- read.table("../ref/pcawg_specimen_histology_August2016_v6.tsv", header=TRUE, sep="\t", comment="", quote="")

s <- strsplit(as.character(finalData$sanger_variant_calling_file_name_prefix),",")
sample2donor <- as.character(finalData$icgc_donor_id[unlist(sapply(seq_along(s), function(i) rep(i, length(s[[i]]))))])
names(sample2donor) <- unlist(s)

donor2type <- factor(specimenData$histology_abbreviation, levels=c(levels(specimenData$histology_abbreviation)[-1], ""))
names(donor2type) <- specimenData$icgc_donor_id
levels(donor2type)[levels(donor2type)==""] <- "Other/NA"


piToTime <- function(timing_param, type=c("Single Gain","CN-LOH", "WGD")){
	type <- match.arg(type)
	pi <-  timing_param[timing_param[,"state"]==1,c("P.m.sX","P.m.sX.lo","P.m.sX.up")]
	pi[1,2:3] <- pi[1,3:2]
	t <- if(type=="Single Gain"){
				3*pi[2,]/(2*pi[2,] + pi[1,])
			}else if(type=="CN-LOH"){
				2*pi[2,]/(2*pi[2,] + pi[1,])
			}else if(type=="WGD"){
				2*pi[2,]/(2*pi[2,] + pi[1,])
			}
	names(t) <- c("","lo","up")
	t[2] <- min(t[1],t[2])
	t[3] <- max(t[1],t[3])
	return(pmin(t,1))
}

bbToTime <- function(bb){
	uniqueBB <- unique(bb)
	sub <- countQueryHits(findOverlaps(bb, bb))
	maj <- bb$major_cn
	min <- bb$minor_cn
	type <- sapply(seq_along(bb), function(i){
				if(maj[i] != 2 | is.na(maj[i])) return(NA)
				type <- if(min[i]==1){ "Single Gain" 
						}else if(min[i]==0){"CN-LOH"}
						else "WGD"
				return(type)
			})
	time <- t(sapply(seq_along(bb), function(i){
				if(sub[i] == 2 | is.na(type[i])) return(c(NA,NA,NA)) # Exclude segments with subclonal CN
				else piToTime(bb$timing_param[[i]],type[i])
			}))
	data.frame(type=factor(type, levels=c("Single Gain","CN-LOH","WGD")), time=time)
}

averageHom <- function(bb){
	sum(width(bb) * (bb$minor_cn == 0) * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
}

.classWgd <- function(ploidy, hom) 2.9 -2*hom <= ploidy

classWgd <- function(bb) .classWgd(averagePloidy(bb), averageHom(bb))

plotBB <- function(bb, ylim=c(0,max(max(bb$total_cn, na.rm=TRUE)))){
	col=RColorBrewer::brewer.pal(4,"Set2")
	s <- c(1:22, "X","Y")
	l <- as.numeric(width(refLengths[seqnames(refLengths) %in% s]))
	names(l) <- s
	plot(NA,NA, ylab="Copy number",xlab="",xlim=c(0,sum(l)), ylim=ylim, xaxt="n")
	c <- cumsum(l)[-length(l)]
	axis(side=1, at=c, labels=rep('', length(l)-1))
	mtext(side=1, at= cumsum(l) - l/2, text=names(l), line=1)
	abline(v=c, lty=3)
	x0 <- start(bb) + cumsum(l)[as.character(seqnames(bb))] - l[as.character(seqnames(bb))]
	x1 <- end(bb) + cumsum(l)[as.character(seqnames(bb))] - l[as.character(seqnames(bb))]
	segments(x0=x0, bb$major_cn ,x1, bb$major_cn, col=col[1], lwd=5* bb$clonal_frequency)
	segments(x0=x0, bb$minor_cn -.125,x1, bb$minor_cn-.125, col=col[2], lwd=5* bb$clonal_frequency)
	segments(x0=x0, bb$total_cn+.125,x1, bb$total_cn+.125, col=1, lwd=5* bb$clonal_frequency)
#	cv <- coverage(bb)
#	cv <- cv[s[s%in%names(cv)]]
#	par(xpd=NA)
#	for(n in names(cv)){
#		cc <- cv[[n]]
#		segments(start(cc) + cumsum(l)[n] - l[n] ,-runValue(cc)/2,end(cc)+ cumsum(l)[n] - l[n], -runValue(cc)/2, col=4)
#	}
	abline(v = chrOffset[1:25], lty=3)
	mtext(side=1, line=1, at=chrOffset[1:24] + diff(chrOffset[1:25])/2, text=names(chrOffset[1:24]))
	legend("topleft", c("Total CN","Major CN","Minor CN"), col=c("black", col[1:2]), lty=1, lwd=2)
}

timeToBeta <- function(time){
	mu <- time[1]
	if(is.na(mu)) return(c(NA,NA))
	if(mu == 1) mu <- 1 - 1e-3
	if(mu == 0) mu <- 1e-3
	v <- (0.5 * (time[3]-time[2]))^2
	alpha <- mu * (mu * (1-mu) / v - 1)
	beta <- (1-mu) *  (mu * (1-mu) / v - 1)
	return(c(alpha, beta))
}

plotTiming <- function(bb, time, col=RColorBrewer::brewer.pal(5,"Set2")[c(3:5)]){
	plot(NA,NA, xlab='', ylab="Time", ylim=c(0,1), xlim=c(0,chrOffset["MT"]), xaxt="n")
	for(i in seq_along(bb)) try({
					s <- start(bb)[i]
					e <- end(bb)[i]
					x <- chrOffset[as.character(seqnames(bb)[i])]
					y <- time[i,2]
					rect(s+x,time[i,3],e+x,time[i,4], border=NA, col=paste0(col[time[i,1]],"88"))
					segments(s+x,y,e+x,y)
				}, silent=TRUE)
	abline(v = chrOffset[1:25], lty=3)
	mtext(side=1, line=1, at=chrOffset[1:24] + diff(chrOffset[1:25])/2, text=names(chrOffset[1:24]))
	legend("topleft", levels(t[,1]), fill=col)
}

source("ComputeMCN.R")
