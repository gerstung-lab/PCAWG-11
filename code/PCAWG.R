#' PCAWG downstream analysis
#' ==========
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

vcfPath <- '../santaCruz/annotated'
basePath <-  '/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/train2/santa_cruz_700_wg_output/final_output/'
dpPath <- paste0(basePath,'/2_subclones/')
mutsigDrivers <- read.table('/nfs/team78pc10/im3/Reference_data/putative_cancer_genes_MutSigCV_5000.txt')$V1
purityPloidy <- read.table(paste0(basePath,'/1_purity_ploidy/purity_ploidy.txt'), header=TRUE, row.names=1)
cnPath <- paste0(basePath,'/4_copynumber/')
bbPath <- paste0(basePath,'/4_copynumber/')

sampleInfo <- read.table("/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/train2/analysis_santa_cruz/samplesheet/output/icgc_sanger_samplesheet_single_samples.txt", header=TRUE, sep="\t", comment.char="")
sampleInfo$sample_type <- sub("/.+","",sub("/nfs/users/nfs_c/cgppipe/pancancer/automated/donors/","", sampleInfo$tumour))

sampleIds <- sub("_mutation_assignments.txt","",dir(dpPath, pattern="_mutation_assignments.txt"))
sampleIds <- intersect(sampleIds, sub("\\..+","",dir(vcfPath, pattern=".bgz$")))

dpFiles <- dir(dpPath, pattern="_subclonal_structure.txt", recursive=TRUE)

loadClusters <- function(ID){
	file <- paste0(dpPath,"/",grep(ID, dpFiles, value=TRUE))
	read.table(file, header=TRUE, sep="\t")
}

allClusters <- lapply(sampleIds, loadClusters)
names(allClusters) <- sampleIds

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

getTumorDepth <- function(vcf){
	if("cavemanVersion" %in% names(meta(header(vcf)))){
		rowSums(getTumorCounts(vcf))
	}else{
		(geno(vcf)$PR + geno(vcf)$NR)[,"TUMOUR"]
	}
}

getAltCount <- function(vcf){
	if("cavemanVersion" %in% names(meta(header(vcf)))){ ## ie subs
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
	file <- paste0(dpPath,"/",ID,"_mutation_assignments.txt")
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

c <- read.table("/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/train2/coverage/coverage_santa_cruz.txt", header=FALSE, sep="\t")
coverage <- c$V2[match(sampleIds, c$V1)]
names(coverage) <- sampleIds
rm(c)

allPosets <- sapply(sampleIds, function(ID){
			reduceToCoverRelations(applyPigeonHole(ID))
		})

#' ### Load all to VCF and annotate
#+ allVcf, cache=TRUE
allVcf <- mclapply(sampleIds, loadVcf, mc.cores=5)

#' Add mutation copy numbers
#+ allVcfMCN, cache=TRUE
allVcf <- mclapply(sampleIds, function(ID) addMutCn(allVcf[[ID]], allBB[[ID]]), mc.cores=5)

names(allVcf) <- sampleIds

#' Remove spurious clusters
for(ID in sampleIds)
	info(allVcf[[ID]])$DPC[!info(allVcf[[ID]])$DPC %in% allClusters[[ID]]$cluster[allClusters[[ID]]$proportion < 1] ] <- NA

#' Classify mutations
for(ID in sampleIds){
	vcf <- allVcf[[ID]]
	class <- rep(2,nrow(vcf))
	class[info(vcf)$DPC < max(allClusters[[ID]]$cluster[allClusters[[ID]]$proportion < 1])] <- 3
	class[info(vcf)$MCN > 1 & class==2] <- 1
	class <- factor(class, levels=1:3, labels=c("early","late","subclonal"))
	class[is.na(info(vcf)$DPC) | is.na(info(vcf)$MCN)] <- NA
	info(allVcf[[ID]])$CLS <- class
	info(header(allVcf[[ID]])) <- rbind(info(header(allVcf[[ID]])), DataFrame(Number="1",Type="String",Description="Mutation classification {early, late, subclonal}", row.names="CLS"))
}

#' Add IDs
for(ID in sampleIds){
	meta(header(allVcf[[ID]])) <- rbind(meta(header(allVcf[[ID]])), DataFrame(Value=ID, row.names="ID"))
}

allMutationTable <- sapply(allVcf, function(vcf) table(info(vcf)$CLS))
o <- order(colSums(allMutationTable))
barplot(allMutationTable[,o]+.5, col=RColorBrewer::brewer.pal(3,"Set1"), border=NA)

#' Output for dN/dS
inigo <- do.call("rbind", sapply(sampleIds, function(i) {x <- allVcf[[i]]; DataFrame(CHR=seqnames(x), POS=start(x), REF=ref(x), ALT=unlist(alt(x)), SAMPLE=i, CLS=info(x)$CLS)}))
write.table(as.data.frame(inigo), file="../santaCruz/santaCruz575SubsDp.txt", sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE)

#' Total mutation frequencies
t <- sapply(allVcf, function(x) table(isDeamination(x), info(x)$CLS), simplify="array")
nMut <- colSums(t,dims=2)
tumourType <- factor(sampleInfo$sample_type[match(sampleIds, sampleInfo$tumour_id)])
mg14:::ggPlot(nMut, tumourType, xlab="", log='y', ylab="All mutations", pch=19)

nDeam <- colSums(t["TRUE",,])
mg14:::ggPlot(nDeam, tumourType, xlab="", log='y', ylab="CpG>TpG", pch=19)

nDeamTrunk <- colSums(t["TRUE",1:2,])
mg14:::ggPlot(nDeamTrunk, tumourType, xlab="", log='y', ylab="CpG>TpG (clonal)", pch=19)

mg14:::ggPlot(realTimeAll, tumourType, xlab="", log='y', ylab="Aging signature", ylim=c(1,200), pch=19)

lBranch <- sapply(allVcf, function(x) branchLengths(x, type="deam"))
lTrunk <- sapply(lBranch, function(x) x[length(x)])

#mg14:::ggPlot(lBranch, tumourType, xlab="", log='y', ylab="Branch lengths")

mg14:::ggPlot(lTrunk, tumourType, xlab="", log='y', ylab="Trunk lengths", pch=19)

powerBranch <- sapply(sampleIds, function(ID) sapply(allClusters[[ID]]$proportion, function(pr) power(pr/purityPloidy[ID,2], round(coverage[ID]))))
avgPower <- mapply(power, purityPloidy[sampleIds,1]/purityPloidy[sampleIds,2], round(coverage), err=1e-3)
avgWeightTrunk <- sapply(allVcf, function(vcf) avgWeights(vcf[na.omit(info(vcf)$CLS!="subclonal")], type="deam"))


#' ### Age
icgc_info <- read.table("../ref/donor.all_projects.tsv", header=TRUE, sep="\t")
tcga_uuid <- read.table("../ref/tcga_uuid2barcode.txt", header=FALSE, sep="\t")
pcawg_info <- read.table("../ref/pcwag_data_release_aug_2015_v1.tsv", header=TRUE, sep="\t")
pcawg_info$tcga_id <- as.character(tcga_uuid$V4[match(pcawg_info$submitter_donor_id, tcga_uuid$V3)])
pcawg_info$age <- icgc_info$donor_age_at_diagnosis[match(ifelse(is.na(pcawg_info$tcga_id), as.character(pcawg_info$submitter_donor_id), pcawg_info$tcga_id), icgc_info$submitted_donor_id)]
rm(tcga_uuid)

age <- pcawg_info$age[match(sampleIds, pcawg_info$sanger_variant_calling_file_name_prefix)]

#' Variation explained
fit <- lm(log(nDeamTrunk) ~ log(age) + tumourType + log(avgPower) + log(avgWeightTrunk))
summary(fit)

fit <- glm(nDeamTrunk ~ log(age) + log(avgPower) + log(1/avgWeightTrunk) + tumourType, family="poisson")
summary(fit)
col <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"))
plot(nDeamTrunk[-fit$na.action], exp(predict(fit)), log='xy', col=col[tumourType[-fit$na.action]], pch=19, xlab="C>T@CpG [observed]", ylab="C>T@CpG [predicted]")
abline(0,1)
plot(age[-fit$na.action], nDeamTrunk[-fit$na.action]/exp(predict(fit, newdata=data.frame(age=1/exp(1), avgPower=avgPower, avgWeightTrunk=avgWeightTrunk, tumourType=tumourType)[-fit$na.action,])),  col=col[tumourType[-fit$na.action]], pch=19, xlab="Age", ylab="Residual C>T@CpG w/o age")
abline(0,1)

#' ### Signatures
snmf <- new.env()
load("../ref/snmf.RData", envir=snmf)
predictRealTime <- function(x, signatures=snmf$snmfFit$P[-1,]){
	snmf$snmfFit$P[1,1] / snmf$weight * snmf$nmSolve(x, signatures, maxIter=100)[1,]
}

set.seed(42)
sigList <- sapply(sampleIds, function(ID) table(info(allVcf[[ID]])$DPC, tncToPyrimidine(allVcf[[ID]])))
sigData <- sapply(sigList, colSums)
n <- 6
nmfFit <- mynmf(sigData,n)
f <- factor(substr(rownames(d),3,5))
b <- barplot(t(nmfFit$P), las=2, col=brewer.pal(8,"Set1"), border=NA,names.arg=rep("", nrow(sigData)))
mtext(at=b,rownames(sigData), side=1, las=2, cex=.7, col=brewer.pal(6,"Set1")[f])

cor(snmf$snmfFit$P[-1,], nmfFit$P[,])
newSig <- nmfFit$P
w <- which.max(cor(snmf$snmfFit$P[-1,1], nmfFit$P[,]))
newSig <- newSig[,c(w, setdiff(1:n,w))]

realTimeBranch <- sapply(sigList, function(x) predictRealTime(t(x), newSig))
realTimeAll <- predictRealTime(sigData, newSig)

fit <- lm(log(age) ~ log(realTimeAll) + tumourType + log(avgPower) + log(avgWeightTrunk))
summary(fit)

#' ### Compute graphs (posets)
toGraph <- function(edgelist, branch.length, edge.labels, node.labels=1:max(edgelist)){
	g <- graph.edgelist(edgelist)
	E(g)$weight <- branch.length
	E(g)$name <- edge.labels
	V(g)$name <- node.labels
	return(g)	
}

na.rm <- function(x) x[!is.na(x)]


allGraphs <- mclapply(sampleIds, function(ID){
			edgelist <- reduceToCoverRelations(applyPigeonHole(ID))
			branch.length <- branchLengths(allVcf[[ID]], type='deam')
			b <- as.numeric(names(branch.length))
			edgelist <- edgelist[edgelist[,2] %in% b & edgelist[,1] %in% b,,drop=FALSE]
			edgelist <- rbind(edgelist, cbind(setdiff(b, edgelist[,1]),max(b)+1))
			edgelist <- edgelist[,2:1, drop=FALSE]
			labels <- sapply(split(unlist(info(allVcf[[ID]])$DG), info(allVcf[[ID]])$DPC), function(x) paste(na.rm(x), collapse=","))
			r <- match(1:max(edgelist),sort(unique(as.numeric(edgelist)), d=T))
			branch.length <- branch.length[as.character(edgelist[,2])]
			labels <- labels[as.character(edgelist[,2])]
			edgelist <- cbind(r[edgelist[,1]], r[edgelist[,2]]) ## reverse node ids
			node.size <- branchLengths(allVcf[[ID]])
			c <- loadClusters(ID)
			alleleFrequency <- c$proportion[match(b, c$cluster)]
			node.labels <- c("Germline", paste0(b, ", n=", node.size, ", f=", round(alleleFrequency,2)))
			g <- toGraph(edgelist, branch.length, edge.labels=labels, node.labels=node.labels)
			V(g)$size <- c(1-purityPloidy[ID,"purity"],alleleFrequency)
			return(g)
		}, mc.cores=4)
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
posetDist <- function(g) {
	e <- get.edgelist(g)
	w <-  c(0,E(g)$weight)
	names(w) <- c("Germline", e[,2])
	d <- shortest.paths(g, mode='out')
	d <- d - rep(w[colnames(d)], each=length(w))/2
	diag(d) <- NA
	d
}

allDist <- sapply(allGraphs, function(g){ #all
			posetDist(graph = g)
		})

#' Map to drivers
allDrivers <- sapply(allVcf, function(x) {sapply(split(unlist(info(x)$DG), info(x)$DPC), na.rm)})
driverTable <- table(unlist(allDrivers))
driverTable <- sort(driverTable, decreasing=TRUE)

minDriver <- 5
recMutations <- names(driverTable)[which(driverTable >= minDriver)]
#allMutations <- allMutations[allMutations!="ARID1A"]
recMutationsGermline <- c("germline",recMutations)

getMutationCluster <- function(allMutations, vcf){
	m <- match(allMutations, unlist(info(vcf)$DG))
	info(vcf)$DPC[m]
}

genotypes <- mclapply(allVcf, function(vcf) sapply(recMutations, function(g) {w <- which(unlist(info(vcf)$DG)==g); if(length(w)==0) 0 else as.numeric(info(vcf)$CLS[w[1]])}), mc.cores=4)
genotypes <- simplify2array(genotypes)
genotypes <- do.call("data.frame", sapply(recMutations, function(g) factor(genotypes[g,], labels=c("none","early","late","subclonal"), levels=0:3), simplify=FALSE))

#' Combine into array
allDistMutations <- array(0, c(rep(length(recMutationsGermline),2), length(sampleIds)), dimnames=list(recMutationsGermline, recMutationsGermline, sampleIds))
for(n in sampleIds){
	c <- getMutationCluster(recMutations, allVcf[[n]])
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
			m <- getMutationCluster(recMutations, allVcf[[ID]])
			l <- branchLengths(allVcf[[ID]], type="deam")
			1/l[m]
		})
rownames(weights) <- recMutations

#' # Fit model
#+ fit, results='asis'
observed <- !is.na(allDistMutations) & !is.infinite(allDistMutations) & c(TRUE,rep(FALSE, length(recMutations))) # only distance to root for the moment..
w <- which(observed, arr.ind = TRUE)
y <- allDistMutations[observed]
x <- CoxHD:::MakeInteger(factor(w[,2], levels=seq_along(recMutationsGermline), labels=recMutationsGermline)) - CoxHD:::MakeInteger(factor(w[,1], levels=seq_along(recMutationsGermline), labels=recMutationsGermline))

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
vcf <- (d1*v1 + d0 * v0)/(d0+d1)  
sd <- sqrt(vcf)
plot(c, r, xlab="Time [C>T@CpG]", yaxt='n', pch=19, cex=sqrt(table(u)[gsub("`","",names(r))]/2), xlim=pmax(0,range(c(pmin(5000,c+sd), c-sd))), ylab="")
segments(pmax(0,c-sd),r, pmin(5000,c+sd),r)
text(pmax(1,c+sd),r, paste0(names(c), ", n=",table(u)[names(c)]), pos=4)


#' Distributions
boxplot(t(allDistMutations[1,-1,][order(c),]), las=2, horizontal=TRUE)

#' ### Pairwise precendences
allPrecMutations <- rowSums(!is.na(allDistMutations) & !is.infinite(allDistMutations), dims=2)
Matrix(allPrecMutations[recMutations, recMutations])

#' ### Intervals

distAsRange <- function(g){
	e <- get.edgelist(g)
	w <-  c(0,E(g)$weight)
	names(w) <- c("Germline", e[,2])
	d <- shortest.paths(g, mode='out')
	y <- shift(IRanges(-w[colnames(d)],0), d["Germline", ])
	names(y) <- paste0(colnames(d), ", genes=",E(g)$name[match(colnames(d), e[,2])])
	y
}

allRanges <- Reduce("c", sapply(allGraphs, distAsRange))

allRangesDrivers <- sapply(recMutations, function(d) allRanges[grep(d, names(allRanges))])

i <- 1
d <- names(driverTable[driverTable >= 5])
col <- rainbow(length(d),  v=0.8)
mm <- matrix(numeric(0), nrow=0, ncol=2)
for(n in d){
	c <- coverage(allRangesDrivers[[n]])
	if(i==1) plot(c, col=col[i], type="S", log='y')
	else lines(c, col=col[i], type="S")
	m <- round(round(max(c)/2))
	w <- which(c==m)[1]
	mm <- rbind(mm, c(w,m))
	i <- i+1
}
points(mm[,1], mm[,2], col=col, pch=(15:19)[seq_along(d) %% 5 + 1])

legend("topright", d, lty=1, col=col, text.font=3, pch=(15:19)[seq_along(d)%%5 +1] , ncol=2)


#' Power
power <- function(f,n, theta=6.3, err=1e-4) sum((log10(dbinom(0:n, n, 0:n/n) / dbinom(0:n,n,err)) > theta) * dbinom(0:n,n,f))

#' ### WGD
wgd <- sapply(sampleIds, function(ID){
			if(purityPloidy[ID,2]<3)
				return(rep(NA,2))
			vcf <- allVcf[[ID]]
			t <- table(info(vcf)$MCN, info(vcf)$TCN, isDeamination(vcf), abs(info(vcf)$CNF - purityPloidy[ID,"purity"]) < 0.01)
			c(t["2","4",2, "TRUE"], rowSums(t["2",,2,"TRUE",drop=FALSE]))
		})

#+ wdgTime, 3,3
o <- order(colMeans(wgd), na.last=NA)
i <- seq_along(o)
plot(wgd[,o[i]], rep(i,each=2), xlab="Time [C>T@CpG]", ylab="", pch=rep(c(1,19), 25))
segments(wgd[1,o[i]], i, wgd[2,o[i]], i)
legend("bottomright", c("MCN=2; CN=4", "MCN=2"), pch=c(1,19))


wgdWeight <- sapply(sampleIds, function(ID){
			if(is.na(wgd[1,ID])) return(NA)
			bb <- allBB[[ID]]
			ix <- abs(bb$clonal_frequency - purityPloidy[ID, "purity"]) < 0.01 & bb$minor_cn==2 & bb$copy_number==4
			sum(as.numeric(width(bb[ix]))) / 3e9
		})



glm(nDeamTrunk ~ offset(log(1/avgWeightTrunk*wgd[1,]/wgdWeight)), subset=wgdWeight>0, family="poisson")

plot(nDeamTrunk*avgWeightTrunk, wgd[1,]/wgdWeight); abline(0,1)
quantile((wgd[1,]/wgdWeight)/(nDeamTrunk*avgWeightTrunk), na.rm=TRUE)

fit <- glm(nDeamTrunk ~ log(age) + log(avgPower) + log(1/avgWeightTrunk) + tumourType, family="poisson")
plot(sort(wgd[1,]/wgdWeight)/exp(coef(fit)[1])); abline(0,1)

p <- predict(fit, newdata=data.frame(age=1/exp(1), avgPower=avgPower, avgWeightTrunk=1, tumourType=tumourType)[-fit$na.action,], type='response', se.fit=TRUE)
plot(age[-fit$na.action], (wgd[1,]/wgdWeight)[-fit$na.action] / p$fit, ylim=c(0,100), xlab="Age at diagnosis", ylab="Infered age during WGD"); abline(0,1)
e <- (t(sapply(wgd[1,], function(w) qpois(c(0.025,0.975), w)))/wgdWeight)[-fit$na.action,] / p$fit
segments(age[-fit$na.action], e[,1], age[-fit$na.action], e[,2])



library(lme4)
y <- c(wgd[1,],nDeamTrunk)
z <- factor(rep(1:length(nDeamTrunk),2))
x <- c(rep(1, length(nDeamTrunk)), 1/avgWeightTrunk)
u <- c(rep(c(0,1), each=length(nDeamTrunk)))
fit <- glm(y ~ z + offset(log(x)) + u, family="poisson")
summary(fit)

#' WGD test
wgdTest <- function(vcf){
	ID <- meta(header(vcf))["ID","Value"]
	bb <- allBB[[ID]]
	if(is.null(bb)) return(NA)
	w <- which(abs(bb$clonal_frequency - purityPloidy[ID, "purity"]) < 0.01 & bb$minor_cn==2 & bb$copy_number==4)
	if(length(w) == 0) return(1)
	ix <- sapply(info(vcf)$CNID, function(x) ifelse(length(x)==0, NA, x %in% w)) & isDeamination(vcf)
	t <- table(as.factor(unlist(seqnames(vcf)))[ix], factor(info(vcf)$MCN, levels=0:2)[ix])[c(1:22, "X"),c("1","2")]
	t <- t[rowSums(t)!=0,,drop=FALSE]
	if(nrow(t) <=1) return(1)
	mcn <- factor(rep(colnames(t), each=nrow(t))) 
	chr <- factor(rep(rownames(t), ncol(t)))
	fit2 <- glm(as.numeric(t) ~ chr*mcn, data=data.frame(chr, mcn), family="poisson")
	fit1 <- glm(as.numeric(t) ~ chr+mcn, data=data.frame(chr, mcn), family="poisson")
	anova(fit2, fit1, test="Chisq")
}

sapply(allVcf[1:100], wgdTest)

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
boxplot(driverTable(driverCnTime[o,][w[o],]), horizontal=TRUE, las=2, col=gisticCn$Type[o][w[o]], xlab="pi0")

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


tcgaInfo <- read.table("../ref/patient_metadata.tsv", header=TRUE, sep="\t", fill=TRUE)
icgcInfo <- read.table("../ref/donor.all_projects.tsv", header=TRUE, sep="\t", fill=TRUE)
pcawgInfo <- read.table("../ref/pcwag_data_release_aug_2015_v1.tsv", header=TRUE, sep="\t")


#' Load vcf
indelPath <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/indel/indelAnnotated"
i <- dir(indelPath, pattern="*.vcf.gz")

testIndel <- function(vcf) sapply(info(vcf)$VC, function(x) if(length(x) ==0) FALSE  else any(x %in% c('frameshift','inframe','ess_splice','SO:0001577:complex_change_in_transcript', 'SO:0001582:initiator_codon_change', 'splice_region')))

indelVcf <- mclapply(sampleIds, function(ID){
			f <- grep(ID, i, value=TRUE)
			v <- readVcf(paste0(indelPath,"/",f), genome="GRCh38")
			#v <- addDriver(v, mutsigDrivers)
			#v <- addMutCn(v, allBB[[ID]])
			#if(nrow(v) > 0) v[testIndel(v)] else v
		}, mc.cores=4)
names(indelVcf) <- sampleIds


# Assign indels to cluster
indelVcf <- sapply(sampleIds, function(ID){
			f <- grep(ID, i, value=TRUE)
			vcf <- indelVcf[[ID]]
			vcf <- addDriver(vcf, mutsigDrivers)
			vcf <- addMutCn(vcf, allBB[[ID]])
			i = header(vcf)@header$INFO
			exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="Numeric",Description="DP cluster", row.names="DPC"))
			i = header(vcf)@header$INFO
			exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="Numeric",Description="DP cluster probability", row.names="DPP"))
			# add clusters
			if(nrow(vcf)==0){
				info(vcf)$DPC <- info(vcf)$DPP <- numeric(0)
				return(vcf)
			}
			altCounts <-  getAltCount(vcf)
			
			tumourDepth <- getTumorDepth(vcf)
			clusters <- allClusters[[ID]]$proportion
			t <- info(vcf)$TCN
			dp <- sapply(seq_along(altCounts), function(i){
						if(is.na(t[i])|t[i]==0) return(c(NA,NA))
						prob <- sapply(allClusters[[ID]]$proportion, function(p) sum(dbinom(x=altCounts[i], size=max(altCounts[i],tumourDepth[i]), prob = pmin(p,1) * (1:t[i])/t[i]))) * allClusters[[ID]]$n_ssms
						prob <- prob/sum(prob)
						w <- which.max(prob)
						c(allClusters[[ID]]$cluster[w],prob[w])
					})
			info(vcf)$DPC <- as.integer(dp[1,])
			info(vcf)$DPP <- dp[2,]
			
			# classify
			class <- rep(2,nrow(vcf))
			class[info(vcf)$DPC < max(allClusters[[ID]]$cluster[allClusters[[ID]]$proportion < 1])] <- 3
			class[info(vcf)$MCN > 1 & class==2] <- 1
			class <- factor(class, levels=1:3, labels=c("early","late","subclonal"))
			class[is.na(info(vcf)$DPC) | is.na(info(vcf)$MCN)] <- NA
			info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number="1",Type="String",Description="Mutation classification {early, late, subclonal}", row.names="CLS"))
			info(vcf)$CLS <- class
			
			return(vcf)
		})


lTotal <- nDeamTrunk*avgWeightTrunk + (nDeam - nDeamTrunk) *2/purityPloidy[sampleIds,"ploidy"]
		

## Time drivers rel to CN
timeCn <- Vectorize(function(ID, cnid){
	# gains only
	bb <- allBB[[ID]]
	vcf <- allVcf[[ID]]
	r <- nDeamTrunk[ID] * avgWeightTrunk[ID] / 3e9
	if(!is.na(wgd[1,ID]))
		return(wgd[1,ID]/wgdWeight[ID])
	if(bb$copy_number[cnid]==3){
		(sum(isDeamination(vcf) & info(vcf)$cnid == cnid & info(vcf)$MCN==2)) * 3e9/width(bb)[cnid]
	}else if(bb$copy_number[cnid]==4){
		(sum(isDeamination(vcf) & info(vcf)$cnid == cnid & info(vcf)$MCN==2)) * 3e9/width(bb)[cnid]
	}else if(bb$copy_number[cnid]==2 & bb$major_cn[cnid]==1){
		0
	}else if(bb$copy_number[cnid]==2 & bb$major_cn[cnid]==2){
		(sum(isDeamination(vcf) & info(vcf)$cnid == cnid & info(vcf)$MCN==2)) * 3e9/width(bb)[cnid]
	}else return(NA)
},vectorize.args = "cnid")


timeDriverCn <- function(ID){
	do.call("c",lapply(c("sub", "indel"), function(type){
						vcf <- if(type=="sub") allVcf[[ID]] else indelVcf[[ID]]
						rgs <- GRanges(genes=CharacterList(), type=character(), sample=character())
						if(nrow(vcf)==0) return(rgs)
						w <- which(na.omit(sapply(!is.na(info(vcf)$DG),any) & info(vcf)$CLS != 'subclonal'))
						cnid <- sapply(info(vcf)$CNID[w], function(c) {
									if(length(c)==1) return(c)
									else c[which.max(allBB[[ID]][c]$clonal_frequency)]
								})
						if(length(w)>0){
							t <- timeCn(ID, cnid)
							v <- info(vcf)$MCN[w]!=1
							s <- ifelse(v, 0, t)
							e <- ifelse(v, t, nDeamTrunk[ID]*avgWeightTrunk[ID])
							cls <- as.character(info(vcf)$CLS[w])
							cls[allBB[[ID]]$copy_number[cnid]==2 & allBB[[ID]]$major_cn[cnid]==1] <- "diploid"
							n <- is.na(s) | is.na(e) | is.infinite(t)
							if(sum(!n)>0) 
								rgs <- c(rgs, GRanges(factor(cls[!n], levels=c(levels(info(vcf)$CLS), "diploid")), IRanges(pmin(s[!n],e[!n]),e[!n]), genes=info(vcf)$DG[w][!n], type=rep(type, sum(!n)), sample=rep(ID, sum(!n))))
						}
						w <- which(na.omit(sapply(!is.na(info(vcf)$DG),any) & info(vcf)$CLS == 'subclonal'))
						if(length(w)>0){
							s <- rep(nDeamTrunk[[ID]]*avgWeightTrunk[ID], length(w))
							e <- rep(nDeamTrunk[[ID]]*avgWeightTrunk[ID] + (nDeam[[ID]] - nDeamTrunk[[ID]]) *2/purityPloidy[ID,"ploidy"], length(w))
							rgs <- c(rgs, GRanges(rep('subclonal', length(w)), IRanges(s,e), genes=info(vcf)$DG[w], type=rep(type,length(w)), sample=rep(ID, length(w))))
						}
						return(rgs)
}))}

driverTimes <- GRangesList(sapply(sampleIds, timeDriverCn))

u <- unlist(driverTimes)
t <- table(unlist(u$genes), u$type, as.character(seqnames(u)))
barplot(sort(rowSums(t), d=TRUE), las=2)
names(age) <- sampleIds
par(mfrow=c(2,3))

pdf("times.pdf",3,3, pointsize=8); par(bty="n", mgp=c(2,0.5,0), mar=c(3,3,1,1))
#### Age of WGD!!!!
y <- (wgd[1,]/wgdWeight)/(nDeamTrunk*avgWeightTrunk)*age
plot(age, y, ylim=c(0,100), xlab="Age at diagnosis", ylab="Age of WGD"); abline(0,1)
e <- (t(sapply(wgd[1,], function(w) qpois(c(0.025,0.975), w)))/wgdWeight) /(nDeamTrunk*avgWeightTrunk)*age
segments(age, e[,1], age, e[,2], col='grey')
d <- density(na.omit(age)); lines(d$x, d$y*.1/max(d$y) * 100)
d <- density(na.omit(y)); lines(d$y*.1/max(d$y) * 100, d$x)
legend("topleft", legend=sum(!is.na(wgd[1,])), lty=1, bty="n", col=1)

clr <- c(col[1:3], "#888888")
#### Drivers
for(gene in c("TP53","PIK3CA","KRAS","NOTCH1","CDKN2A")){
t <- u[sapply(u$genes, function(x) gene %in% x)]
s <- age[t$sample]/lTotal[t$sample]
y <- (end(t)+start(t)) /2*s

plot(age[t$sample], y , xlim=c(0,100), ylim=c(0,100), col=clr[as.numeric(seqnames(t))], xlab="Age at diagnosis", ylab="Age at mutation", main=gene, pch=".")
abline(0,1)
e <- cbind(start(t),end(t)) *s
segments(age[t$sample], e[,1], age[t$sample], e[,2], col=clr[as.numeric(seqnames(t))])
f <- function(x) paste(names(x), x)
legend("topleft", legend=f(table(seqnames(t))), lty=1, col=clr[1:4], bty="n")
d <- density(na.omit(age[t$sample])); lines(d$x, d$y*.1/max(d$y) * 100)
w <- !is.na(e[,1])
d <- sapply(coverage(GRanges(seqnames(t)[w], IRanges(e[w,1],e[w,2])), weight=1/(e[w,2]-e[w,1]+1)), function(x) {y <- numeric(100); n <- as.numeric(x); y[1:min(length(n),100)] <- n[1:min(length(n),100)];y})
a <- apply(d, 1,cumsum)/length(t)*500
for(i in nrow(a):1) polygon(c(0,a[i,]), c(1,1:100), col=mg14::colTrans(clr[i]), border=NA, )
}
dev.off()


## inigo selection
sSub <- cbind(wMIS =	c(1.91831062196026,	0.919545485202377,	4.00188538962342),
wNON =	c(6.73416105169359,	2.38927228182595,	18.9802248220492),
wSPLICE	= c(4.73003260101133,	1.20999550398506,	18.4903235862819))


sEarly <- cbind(wMIS=	c(2.32117247053749,	1.67002387515843,	3.22620635436723),
wNON	=c(9.95274974928443,	6.25798469653712,	15.8289341338107),
wSPLICE	=c(3.76229376536249,	2.1063230165632,	6.72017267322149))

sLate	<- cbind(wMIS=c(	1.63512904723344,	1.33839800123583,	1.99764718614178),
wNON=c(	4.93163721684375,	3.65519983206903,	6.65382106476828),
wSPLICE	=c(3.11074216948118,	2.11924768949577,	4.56610942314669))

pdf("selection.pdf",3,3, pointsize=8); par(bty="n", mgp=c(2,0.5,0), mar=c(3,3,1,1))
barplot(rbind(sEarly[1,], sLate[1,],sSub[1,]), col=clr[1:3], beside=TRUE, ylab="dN/dS") -> b
legend("topleft", c("early","late","subclonal"), fil=clr[1:3], bty="n")
segments(b,rbind(sEarly[2,], sLate[2,],sSub[2,]),b,rbind(sEarly[3,], sLate[3,],sSub[3,]))
abline(h=1)
dev.off()
