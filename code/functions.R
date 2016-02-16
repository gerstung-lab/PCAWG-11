# TODO: Add comment
# 
# Author: mg14
###############################################################################

vcfPath <- '/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/subs/annotated'
basePath <-  '/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/dp/2015_11_23_vanloo_wedge_batch06070809_sanger'
dpPath <- paste0(basePath,'/2_subclones/')
mutsigDrivers <- read.table('/nfs/team78pc10/im3/Reference_data/putative_cancer_genes_MutSigCV_5000.txt')$V1
purityPloidy <- read.table(paste0(basePath,'/1_purity_ploidy/purity_ploidy.txt'), header=TRUE, row.names=1)
cnPath <- paste0(basePath,'/4_copynumber/')
bbPath <- paste0(basePath,'/4_copynumber/')

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
	
loadClusters <- function(ID){
	file <- paste0(dpPath,"/",grep(ID, dpFiles, value=TRUE))
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
	file <- gzfile(paste0(bbPath, "/",ID,"_segments.txt.gz"))
	tab <- read.table(file, header=TRUE, sep='\t')
	GRanges(tab$chromosome, IRanges(tab$start, tab$end), strand="*", tab[-3:-1])
}

getTumorCounts <- function(vcf){
	sapply(grep("(F|R).Z", names(geno(vcf)), value=TRUE), function(i) geno(vcf)[[i]][,"TUMOUR"])
}

getTumorDepth <- function(vcf){
	if("cavemanVersion" %in% rownames(meta(header(vcf)))){
		rowSums(getTumorCounts(vcf))
	}else{
		(geno(vcf)$PR + geno(vcf)$NR)[,"TUMOUR"]
	}
}

getAltCount <- function(vcf){
	if("cavemanVersion" %in% rownames(meta(header(vcf)))){ ## ie subs
		c <- getTumorCounts(vcf)
		t <- c[,1:4] + c[,5:8]
		colnames(t) <- substring(colnames(t),2,2)
		a <- as.character(unlist(alt(vcf)))
		sapply(seq_along(a), function(i) t[i, a[i]])
	}else{ ## ie indel
		(geno(vcf)$PP + geno(vcf)$NP + geno(vcf)$PB + geno(vcf)$NB)[,"TUMOUR"]
	}
}

computeMutCn <- function(vcf, bb){
	if(nrow(vcf)==0)
		return(DataFrame(MCN=numeric(),TCN=numeric(),CNF=numeric(),PMCN=numeric(), CNID=numeric()))
	altCount <- getAltCount(vcf)
	tumDepth <- getTumorDepth(vcf)
	names(altCount) <- names(tumDepth) <- NULL
	purity <- purityPloidy[ID, 'purity']
	f <- findOverlaps(vcf, bb)
	copyNum <- split(bb$copy_number[subjectHits(f)], queryHits(f))
	cloneFreq <- split(bb$clonal_frequency[subjectHits(f)], queryHits(f))
	m <- t(sapply(seq_along(altCount), function(i){
						if(!i %in% names(copyNum)) return(rep(NA,4))
						cni <- copyNum[[as.character(i)]]
						cfi <- cloneFreq[[as.character(i)]]
						totalCnTumor <- sum(cni * cfi)
						totalCnNormal <- 2
						cnStates <- unlist(sapply(seq_along(cni), 
										function(j){s <- 0:cni[j] * cfi[j]
											names(s) <- paste0(j,':',0:cni[j]) 
											return(s)}, simplify=FALSE)) * purity / (totalCnTumor * purity + (1-purity) * totalCnNormal)
						prob <- sapply(cnStates, function(s) dbinom(altCount[i],tumDepth[i],s)) + .Machine$double.eps
						prob <- prob/sum(prob)
						w <- which.max(prob)
						idx <- as.numeric(strsplit(names(prob[w]), ":")[[1]])
						names(prob) <- NULL
						c(MCN=idx[2], TCN=cni[idx[1]], CNF=cfi[idx[1]], PMCN=prob[w])
					}))
	colnames(m) <- c("MCN","TCN","CNF","PMCN")
	D <- DataFrame(m, CNID = as(f,"List"))
	colnames(D) <- gsub("[[:punct:]]+[0-9]+","", colnames(D))
	return(D)
}

addMutCn <- function(vcf, bb){
	i = header(vcf)@header$INFO
	exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=c(1,1,1,1,"."),Type=rep("Numeric",5), Description=c("Mutation copy number","Total copy number","Copy number frequency (relative to all cancer cells)", "MCN probability","BB segment ids"), row.names=c("MCN","TCN","CNF","PMCN","CNID")))
	info(vcf) <- cbind(info(vcf), computeMutCn(vcf, bb))
	return(vcf)
}

classifyMutations <- function(vcf){
	info(vcf)$MCN
}

tryExceptNull <- function(x) if(class(x)=="try-error") GRanges() else x

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

testDriver <- function(vcf) sapply(info(vcf)$VC, function(x) if(length(x) ==0) FALSE  else x %in% c('nonsense','missense','ess_splice','frameshift','inframe','cds_distrupted'))

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
	info(vcf)$MCN / info(vcf)$TCN
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



