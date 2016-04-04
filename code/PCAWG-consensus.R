#' ---
#' title: PCAWG-11 downstream analysis
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 4
#'     number_sections: true
#'     auto_identifiers: true
#'     table_captions: true
#' author: Moritz Gerstung
#' ---

#' # PCAWG-11 downstream analysis

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
library(Matrix)
library(CoxHD)
library(igraph)

setwd("/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/code")
source("functions.R")

basePath <-  '/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/dp/2016-03_consensus'
dpPath <- paste0(basePath,'/2_subclones/')
purityPloidy <- read.table(paste0(basePath,'/1_purity_ploidy/purity_ploidy.txt'), header=TRUE, row.names=1)
colnames(purityPloidy) <- c("purity","ploidy")
cnPath <- paste0(basePath,'/4_copynumber/')
bbPath <- paste0(basePath,'/4_copynumber/')

ref <- "/lustre/scratch112/sanger/cgppipe/PanCancerReference/genome.fa.gz" #meta(header(v))["reference",]
refLengths <- scanFaIndex(file=ref)

mutsigDrivers <- read.table('../ref/putative_cancer_genes_MutSigCV_5000.txt')$V1
census <- read.csv('../ref/Census_allThu Mar 24 17_30_12 2016.csv.txt', stringsAsFactors=FALSE)
census <- census$Gene.Symbol[grepl("Mis|F|N|M|S",census$Mutation.Types)]
cancerGenesS <- sort(union(mutsigDrivers, census))

sampleInfo <- read.table("/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/train2/analysis_santa_cruz/samplesheet/output/icgc_sanger_samplesheet_single_samples.txt", header=TRUE, sep="\t", comment.char="")
sampleInfo$sample_type <- sub("/.+","",sub("/nfs/users/nfs_c/cgppipe/pancancer/automated/donors/","", sampleInfo$tumour))

sampleIds <- sub("_mutation_assignments.txt.gz","",dir(dpPath, pattern="_mutation_assignments.txt.gz"))
sampleIds <- intersect(sampleIds, sub("\\..+","",dir(vcfPath, pattern="*.bgz")))

dpFiles <- dir(dpPath, pattern="_subclonal_structure.txt", recursive=TRUE)

allClusters <- mclapply(sampleIds, loadClusters, mc.cores=8)
names(allClusters) <- sampleIds

ploidy <- read.table("/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/icgc_pancan_full/processed_data/copynumberConsensus/consensus_001/ploidy.tsv", header=TRUE, sep="\t", row.names=1)
purity <- read.table("/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/icgc_pancan_full/processed_data/copynumberConsensus/consensus_001/purity/consensus_purities.20160319.txt", header=TRUE, sep="\t", row.names=1)
s <-  setdiff(sampleIds, rownames(purityPloidy))
purityPloidy <- rbind(purityPloidy, data.frame(purity=purity[s,1],ploidy=ploidy[s,1], row.names=s))

#' ### Battenberg
#+ allBB, cache=TRUE
allBB <- mclapply(sampleIds, loadBB, mc.cores=8) 
names(allBB) <- sampleIds

#' Add consensus CN where BB is missing
for(w in which(sapply(allBB, length)==0)){
	cnPath <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/cn/consensus_001/all"
	ID <- names(allBB)[w]
	file <- paste0(cnPath, "/",ID,"_segments.txt")
	tab <- read.table(file, header=TRUE, sep='\t')
	allBB[[w]] <- GRanges(tab$chromosome, IRanges(tab$start, tab$end), strand="*", tab[-3:-1])
}
			
#' ### Subs 
#' Load all annotated Subs files
#+ allVcf, cache=TRUE, cache.lazy=FALSE
allVcf <- mclapply(sampleIds, function(ID){
			readVcf(grep(ID, dir(vcfPath, pattern=".complete_annotation.vcf.bgz", full.names=TRUE), value=TRUE), genome="GRCh37")
		}, mc.cores=8)
names(allVcf) <- sampleIds

#' Reload a few (not run)
#+ eval=FALSE
#w <- setdiff(sampleIds, sub("_.+","",dir(bbPath, pattern="*_segments.txt.gz"))) 
w <- setdiff(sub("\\..+","",dir(vcfPath, pattern=".complete_annotation.vcf.bgz")), sampleIds)
tmp <- mclapply(w, function(ID){
			readVcf(grep(ID, dir(vcfPath, pattern=".complete_annotation.vcf.bgz", full.names=TRUE), value=TRUE), genome="GRCh37")
		}, mc.cores=8)
names(tmp) <- w
for(s in w)
	allVcf[[s]] <- tmp[[s]]
rm(tmp)
sampleIds <- names(allVcf)

#' ### Indels
#' Load indels after exoneration to get reliable VAF estimates
indelPath <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/indel/indelVafAnnotatedConsensus"
i <- dir(indelPath, pattern="*.vcf.gz$")

#' Load annotated vcfs
#+ indelVcf, cache=TRUE
indelVcf <- mclapply(sampleIds, function(ID){
			f <- grep(ID, i, value=TRUE)
			v <- readVcf(paste0(indelPath,"/",f), genome="GRCh37")
			#v <- addDriver(v, mutsigDrivers)
			#v <- addMutCn(v, allBB[[ID]])
			#if(nrow(v) > 0) v[testIndel(v)] else v
		}, mc.cores=8)
names(indelVcf) <- sampleIds

#' Fix empty vcfs
for(w in which(sapply(indelVcf, ncol)==0))
	indelVcf[[w]] <- indelVcf[[1]][0,]


#' Compute mutation copy numbers and assign indels to subclonal clusters
#+ indelAnnotated, cache=TRUE
indelAnnotated <- mclapply(sampleIds, function(ID){
			vcf <- indelVcf[[ID]]
			meta(header(vcf)) <- rbind(meta(header(vcf)), DataFrame(Value=ID, row.names="ID"))
			rownames(info(header(vcf))) <- sub("NA","NOA",rownames(info(header(vcf))) )
			colnames(info(vcf)) <- sub("NA","NOA",colnames(info(vcf)))
			
			vcf <- addDriver(vcf, cancerGenes)
			
			exptData(vcf)$header@header$INFO <- rbind(header(vcf)@header$INFO, DataFrame(Number=1,Type="Integer",Description="DP cluster", row.names="DPC"))
			exptData(vcf)$header@header$INFO <- rbind(header(vcf)@header$INFO, DataFrame(Number=1,Type="Float",Description="DP cluster probability", row.names="DPP"))
			info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number="1",Type="String",Description="Mutation classification: {clonal [early/late/NA], subclonal}", row.names="CLS"))
			# add clusters
			if(nrow(vcf)==0){
				info(vcf)$DPC <- info(vcf)$DPP <-  numeric(0)
				info(vcf)$CLS <- character(0)
				return(vcf)
			}
			meta(header(vcf)) <- rbind(meta(header(vcf)), DataFrame(Value=ID, row.names="ID"))
			
			clusters <- allClusters[[ID]]
			# Add mutation copy numbers
			vcf <-  addMutCn(vcf, allBB[[ID]], clusters)
			
			info(vcf)$DPC <- sapply(info(vcf)$CNF, function(f) {if(is.null(f) | is.na(f)) NA else {d <- abs(f-clusters$proportion); w <- which.min(d); if(d[w]>0.05) NA else clusters$cluster[w]}})
			info(vcf)$DPP <- rep(NA, nrow(vcf))
			
			# Remove spurious clusters
			info(vcf)$DPC[!info(vcf)$DPC %in% clusters$cluster[clusters$proportion < 1] ] <- NA
			
			# Classify mutations
			class <- rep(2,nrow(vcf))
			class[info(vcf)$DPC < max(clusters$cluster[clusters$proportion < 1])] <- 4
			class[info(vcf)$MCN > 1 & class==2 ] <- 1
			i <- info(vcf)
			class[ (i$MJCN == 1 | i$MNCN == 1) & i$MCN == 1] <- 3
			class <- factor(class, levels=1:4, labels=c("clonal [early]","clonal [late]","clonal [NA]", "subclonal"))
			class[is.na(info(vcf)$MCN)] <- NA
			
			info(vcf)$CLS <- class
			
			return(vcf)
		}, mc.cores=8)


#' ## Genotypes
#' Compute all genotypes, including zygousity
#+ genotypes, cache=TRUE
subGenotypes <- simplify2array(mclapply(allVcf,  getGenotype, mc.cores=8, useNA='always'))
indelGenotypes <- simplify2array(mclapply(indelAnnotated,  getGenotype, mc.cores=8, useNA='always'))

allGenotypes <- aperm(abind::abind(subs=subGenotypes,indels=indelGenotypes, along=5), c(1,5,2,3,4))

#t <- (apply(allGenotypes, c(1,3), sum) + apply(indelGenotypes, c(1,3), sum))

#' Most abundant
#+ genesFrequency, fig.width=7
t <- t(asum(allGenotypes[cancerGenesS,,,,], c(2,4,5)))
barplot(t[,order(-colSums(t))[1:100]], col=RColorBrewer::brewer.pal(5,"Set1"), legend=TRUE ,las=2)

#' Normalised by CDS length
#+ genesDensity, fig.width=7
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
entrez <- select(org.Hs.eg.db, cancerGenesS, "ENTREZID", "SYMBOL")
txids <- select(TxDb.Hsapiens.UCSC.hg19.knownGene, na.omit(entrez$ENTREZID), "TXID", "GENEID")$TXID
cdsByGene <- cdsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")[na.omit(entrez[,2])]
cdsLength <- sapply(cdsByGene, function(x) sum(width(x)))
names(cdsLength) <- entrez$SYMBOL[match(names(cdsLength), entrez$ENTREZID)]

n <- t/rep(cdsLength[match(cancerGenesS, names(cdsLength))], each=5)
barplot(n[,order(-colSums(n))[1:100]], col=RColorBrewer::brewer.pal(5,"Set1"), legend=TRUE ,las=2)

n <- n[,!is.na(n[1,])]

#' Cumulative effects
#+ genesCumulative
tt <- t(n[,order(-colSums(n))])
RColorBrewer::brewer.pal(10,"Spectral") -> s
m <- as.matrix(aggregate(tt/rep(colSums(tt), each=nrow(tt)), list(round(log2(1:nrow(tt)))), sum))
barplot(m[,-1], col=s, legend.text=2^m[,1])

set1 <- RColorBrewer::brewer.pal(9, "Set1")
r <- tt/rep(colSums(tt), each=nrow(tt))
plot(cumsum(r[,1]), col=set1[1], type='l', xlab="cancer gene rank", ylab="fraction of mutations", log='', ylim=c(0,1))
for(j in 2:5) lines(cumsum(r[,j]), col=set1[j])
legend("topleft", col=set1[1:5], lty=1, colnames(tt), bty="n")
abline(h=0.5, lty=3)

#P <- apply(t, 1, function(x) pbinom(x[2], sum(x), prob=1-1712568/(18578183+1712568), lower.tail=TRUE))
#
#data.frame(t,P=P, Q=p.adjust(P, method="BH"), sig=mg14::sig2star(p.adjust(P, method="BH")))


#' ## Burden and signatures

#' ### Age
icgc_info <- read.table("../ref/donor.all_projects.tsv", header=TRUE, sep="\t")
tcga_uuid <- read.table("../ref/tcga_uuid2barcode.txt", header=FALSE, sep="\t")
pcawg_info <- read.table("../ref/pcwag_data_release_aug_2015_v1.tsv", header=TRUE, sep="\t")
pcawg_info$tcga_id <- as.character(tcga_uuid$V4[match(pcawg_info$submitter_donor_id, tcga_uuid$V3)])
pcawg_info$age <- icgc_info$donor_age_at_diagnosis[match(ifelse(is.na(pcawg_info$tcga_id), as.character(pcawg_info$submitter_donor_id), pcawg_info$tcga_id), icgc_info$submitted_donor_id)]
rm(tcga_uuid)

age <- pcawg_info$age[match(sampleIds, pcawg_info$sanger_variant_calling_file_name_prefix)]

tumourType <- factor(sub("-.+","",pcawg_info$dcc_project_code[match(sampleIds, pcawg_info$sanger_variant_calling_file_name_prefix)]))


#' ### Total mutation frequencies, split by deaminations at CpG
#+ tabDeam, cache=TRUE
tabDeam <- simplify2array(mclapply(allVcf, function(x) table(isDeamination(x), classifyMutations(x)), mc.cores=8))

#+ tabDeamPlot, fig.width=7
nMut <- colSums(tabDeam,dims=2)
tumourType <- factor(sub("-.+","",sampleInfo$sample_type)[match(sampleIds, sampleInfo$tumour_id)])
mg14:::ggPlot(nMut, tumourType, xlab="", log='y', ylab="All mutations", pch=19)

nDeam <- colSums(tabDeam["TRUE",,])
mg14:::ggPlot(nDeam, tumourType, xlab="", log='y', ylab="CpG>TpG", pch=19)

nDeamTrunk <- colSums(tabDeam["TRUE",1:2,])
mg14:::ggPlot(nDeamTrunk, tumourType, xlab="", log='y', ylab="CpG>TpG (clonal)", pch=19)

#mg14:::ggPlot(realTimeAll, tumourType, xlab="", log='y', ylab="Aging signature", ylim=c(1,200), pch=19)

lBranch <- sapply(allVcf, function(x) branchLengths(x, type="deam"))
lTrunk <- sapply(lBranch, function(x) x[length(x)])

#mg14:::ggPlot(lBranch, tumourType, xlab="", log='y', ylab="Branch lengths")

mg14:::ggPlot(lTrunk, tumourType, xlab="", log='y', ylab="Trunk lengths", pch=19)

#powerBranch <- sapply(sampleIds, function(ID) sapply(allClusters[[ID]]$proportion, function(pr) power(pr/purityPloidy[ID,2], round(coverage[ID]))))
#avgPower <- mapply(power, purityPloidy[sampleIds,1]/purityPloidy[sampleIds,2], round(coverage), err=1e-3)
avgWeightTrunk <- unlist(mclapply(allVcf, function(vcf) avgWeights(vcf[na.omit(info(vcf)$CLS!="subclonal")], type="deam"), mc.cores=8))


#' Variation explained
fit <- lm(log(nDeamTrunk + 1) ~ log(age) + tumourType + log(avgWeightTrunk))
summary(fit)

fit <- glm(nDeamTrunk ~ log(age) + log(1/avgWeightTrunk) + tumourType, family="poisson")
summary(fit)
col <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"))
plot(nDeamTrunk[-fit$na.action], exp(predict(fit)), log='xy', col=col[tumourType[-fit$na.action]], pch=19, xlab="C>T@CpG [observed]", ylab="C>T@CpG [predicted]")
abline(0,1)
plot(age[-fit$na.action], nDeamTrunk[-fit$na.action]/exp(predict(fit, newdata=data.frame(age=1/exp(1), avgWeightTrunk=avgWeightTrunk, tumourType=tumourType)[-fit$na.action,])),  col=col[tumourType[-fit$na.action]], pch=19, xlab="Age", ylab="Residual C>T@CpG w/o age")
abline(0,1)

#' ### COSMIC Signatures
#' #### Run all 30 signatures
#+ sigDecomp30, cache=TRUE
signatures <- read.table("http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt", header=TRUE, sep="\t")
sigTable <- simplify2array(mclapply(allVcf, function(vcf) table(classifyMutations(vcf), tncToPyrimidine(vcf)), mc.cores=8))
sigTable <- aperm(sigTable, c(2,3,1))
dimnames(sigTable)[[2]] <- sampleIds
S <- as.matrix(signatures[match(dimnames(sigTable)[[1]],as.character(signatures[,3])),1:30+3])
rownames(S) <- dimnames(sigTable)[[1]]
w <- which(colSums(sigTable) != 0)
nmfFit <- nmSolve(matrix(sigTable, nrow=96)[,w], S, maxIter=10000, tol=1e-5)
M <- matrix(0,nrow=30, ncol=prod(dim(sigTable)[2:3]))
M[,w] <- nmfFit
sigDecomp30 <- array(M, dim=c(30, dim(sigTable)[2:3]))
rm(M)
dimnames(sigDecomp30) <- c(list(colnames(signatures[,1:30+3])),dimnames(sigTable)[2:3])

#+ sigDecomp30Plot, fig.width=7
barplot(t(sigDecomp30[,1,]))

sigClonal <- t(sigDecomp30[,,1]+sigDecomp30[,,2])

mg14:::ggPlot(sigClonal[,1] +1, tumourType, xlab="", log='y', ylab="Signature 1 (clonal)", pch=19, ylim=c(1,1e5), las=2)

fit <- VGAM::vglm(round(sigClonal) ~ log(age) + log(1/avgWeightTrunk) + tumourType, family="poissonff")
summary(fit)

fit <- lm(log(sigClonal[,1] +1) ~ log(age) + tumourType + log(avgWeightTrunk))
summary(fit)

fit <- glm(age ~  log(1/avgWeightTrunk) + tumourType + log(sigClonal + 1), family=gaussian(link="log"))
summary(fit)


#' #### Active signatures
#+ eval=FALSE
library(png)
p <- readPNG("~/Desktop/matrix.png")
w <- p[ seq(0,2156,l=32)[-32]+30,seq(5,3583,l=41)[-41]+42,1] < .5
rownames(w) <- paste("Signature",1:31)
colnames(w) <- c("ACCAx","ALL","LAML","BLCA","BRCA","CESC","CHSA","CLLE","COAD","GBM","LGG", "HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","LUSQx","MALY","HLYx","PEME","MELA","MMx","NACA","NBL","ESAD","ORCA","OS","OV","PACA","PARAx","PILOx","PRAD","STAD","THCA","URCAx","UCEC","UCSAx","SKCM")

#' Read signature activities
#+ sigDecomp, cache=TRUE
sigActivity <- read.table("../ref/sigActivity.txt", header=TRUE, row.names=1, sep="\t")
sigActivityHere <- as.matrix(sigActivity)[1:30,match(levels(tumourType), colnames(sigActivity))] 
sigActivityHere[,is.na(colnames(sigActivityHere))] <- rowSums(sigActivity)[1:30] > 10
colnames(sigActivityHere) <- levels(tumourType)
sigActivityHere[8,"PBCA"] <- 1

w <- colSums(sigTable) != 0

M <- matrix(0,nrow=30, ncol=prod(dim(sigTable)[2:3]))
rownames(M) <- colnames(S)


for(t in levels(tumourType)){
	v <- tumourType==t 
	E <- nmSolve(matrix(sigTable, nrow=96)[,which(w & v)], S[,sigActivityHere[,t]==1], maxIter=5000)
	M[rownames(E), which(v&w)] <- E
}

sigDecomp <- array(M, dim=c(30, dim(sigTable)[2:3]))
rm(M)

dimnames(sigDecomp) <- c(list(colnames(signatures[,1:30+3])),dimnames(sigTable)[2:3])

#+ signaturesTumourTypes, fig.width=7
mg14:::ggPlot(sigDecomp[5,,2]+1, tumourType, log='y', ylab="Signature 5")
mg14:::ggPlot(sigDecomp[1,,2]+1, tumourType, log='y', ylab="Signature 1")

#+ signaturesTime, fig.width=7
t <- t(asum(sigDecomp, 2))
barplot(t, col=RColorBrewer::brewer.pal(4,"Set1"), legend=TRUE, las=2)
barplot(t/rep(colSums(t), each=4),col=RColorBrewer::brewer.pal(4,"Set1"), legend=TRUE, las=2)



#' ### TNC probabilities
#+ tncProb, cache=TRUE
tncProb <- sapply(1:96, function(i) S[i,] * sigDecomp, simplify='array')
s <- apply(tncProb, 2:4, sum)
tncProb <- tncProb / rep(s, each=30)
rm(s)
dimnames(tncProb)[[4]] <- rownames(S)


#' ## WGD
#' ### Test for WGD
library(mixtools)
mixmdl = normalmixEM(purityPloidy[,2], k=3)
plot(mixmdl,which=2, n=50)
lines(density(purityPloidy[,2]), lty=2, lwd=2)

isWgd <- sapply(sampleIds, function(ID) !(which.max(mixmdl$posterior[rownames(purityPloidy)==ID,])!=3 | purityPloidy[ID,2] < 2))

#' ### Using C>T @ CpG
#+ wgdDeam, cache=TRUE
wgdDeam <- simplify2array(mclapply(allVcf, function(vcf){
					ID <- meta(header(vcf))["ID","Value"]
			if(which.max(mixmdl$posterior[rownames(purityPloidy)==ID,])!=3 | purityPloidy[ID,2] < 2)
				return(rep(NA,3))
			w <- isDeamination(vcf) & abs(info(vcf)$CNF - purityPloidy[ID,"purity"]) < 0.01
			t <- table(factor(info(vcf)$MCN[w], levels=1:20), factor(info(vcf)$TCN[w], levels=1:20))
			c(t["1","4"], t["2","4"], rowSums(t["2",,drop=FALSE]))
		}, mc.cores=8))
colnames(wgdDeam) <- sampleIds

#+ wdgTime, 3,3
o <- order(colMeans(wgdDeam), na.last=NA)
i <- seq_along(o)
plot(wgdDeam[2:3,o[i]], rep(i,each=2), xlab="Time [C>T@CpG]", ylab="", pch=rep(c(1,19), 25), log='x')
segments(wgdDeam[2,o[i]], i, wgdDeam[3,o[i]], i)
legend("bottomright", c("MCN=2; CN=4", "MCN=2"), pch=c(1,19))


#' #### Direct using age
a <- jitter(age)
plot(a, age * (2*wgdDeam[2,]+1)/(1+2*wgdDeam[2,] + wgdDeam[1,]), col=col[tumourType], pch=ifelse(colSums(wgdDeam[1:2,])==0,NA,19), xlab='Age at diagnosis', ylab="Age at WGD", cex=1.5)
d <- density(na.omit(a))
lines(d$x, 500*d$y)
d <- density(na.omit( ifelse(colSums(wgdDeam[1:2,])==0,NA,age * (2*wgdDeam[2,]+1)/(1+2*wgdDeam[2,] + wgdDeam[1,]))))
lines(500*d$y,d$x)
abline(0,1)
ci <- sapply(1:1000, function(foo){
			n1 <- rpois(n=length(wgdDeam[1,]),lambda=wgdDeam[1,]+1)
			n2 <- rpois(n=length(wgdDeam[2,]),lambda=wgdDeam[2,]+1)
			age * (2*n2)/(2*n2 + n1)
		})

ci <- apply(ci,1, function(x) if(all(is.na(x))) rep(NA,2) else quantile(x,c(.025, 0.975), na.rm=TRUE))

segments(a,ci[1,], a,ci[2,], col=col[tumourType], lwd=ifelse(colSums(wgdDeam[1:2,])==0,NA,1))

tWgd <- (2*wgdDeam[2,]+1)/(1+2*wgdDeam[2,] + wgdDeam[1,])
mg14:::ggPlot(tWgd, tumourType)

#' ### Using TNC
#+ wgdTnc, cache=TRUE
wgdTnc <- simplify2array(mclapply(allVcf[isWgd], function(vcf){
					ID <- meta(header(vcf))["ID","Value"]
					if(which.max(mixmdl$posterior[rownames(purityPloidy)==ID,])!=3 | purityPloidy[ID,2] < 2)
						return(matrix(NA,nrow=30, ncol=3))
					c <- classifyMutations(vcf)
					bb <- allBB[[ID]]
					w <- which(vcf %over% bb[bb$minor_cn==2 & bb$major_cn==2] & c!="subclonal")
					pre <- info(vcf)$MCN[w] == 2 &  info(vcf)$TCN[w] == 4
					post <- info(vcf)$MCN[w] == 1 &  info(vcf)$TCN[w] == 4
					tnc <- tncToPyrimidine(vcf[w])
					
					wgd <- sapply(1:30, function(i){
								c(sum(tncProb[i,ID,1,tnc[pre]]), sum(tncProb[i,ID,2,tnc[post]]))
							})
					
					ci <- sapply(1:1000, function(foo){
								n1 <- rpois(n=length(wgd[1,]),lambda=wgd[1,]+1)
								n2 <- rpois(n=length(wgd[2,]),lambda=wgd[2,]+1)
								(2*n2)/(2*n2 + n1)
							})
					
					ci <- apply(ci,1, function(x) if(all(is.na(x))) rep(NA,2) else quantile(x,c(.025, 0.975), na.rm=TRUE))
					
					cbind(hat=(2*wgd[2,]+0.5)/(1+2*wgd[2,] + wgd[1,]), CI=t(ci))
					
				}, mc.cores=8))


rowMeans((1-wgdTnc[,1,])*rep(age[isWgd],30), na.rm=TRUE)



#' ## Selection
#' ### Output for dN/dS
#+ dNdSOut, eval=FALSE
forInigo <- do.call("rbind", mclapply(allVcf, function(vcf) {DataFrame(CHR=seqnames(vcf), POS=start(vcf), REF=ref(vcf), ALT=unlist(alt(vcf)), SAMPLE=meta(header(vcf))["ID",], CLS=classifyMutations(vcf), MCN=info(vcf)$MCN, TCN=info(vcf)$TCN)}, mc.cores=8))
write.table(as.data.frame(forInigo), file="../scratch/mar16_1429samples.txt", sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE)
rm(forInigo)

#' ### From Inigo
#+ dNdS, fig.width=4
dNdS <- sapply(grep("dNdS_357drivers",dir("/nfs/users/nfs_c/cgppipe/pancancer/workspace/im3/PCAWG_Evolution_and_Heterogeneity/dNdS_vs_clonality/BARCELONA_dataset", pattern="dNdSvalues_with_CI95.txt", full.names=TRUE, recursive=TRUE), value=TRUE), function(x) as.matrix(read.table(x, sep="\t", row.names=1, header=TRUE)), simplify='array')
dimnames(dNdS)[[3]] <- sub("(.+drivers_)(.+)(.dNdScv.+)","\\2",dimnames(dNdS)[[3]])
dimnames(dNdS)[[1]] <- c(wMIS="missense",wNON="nonsense", wSPL="splice site")[dimnames(dNdS)[[1]] ]
col <- sapply(c("#444444",RColorBrewer::brewer.pal(9,"Set1")), function(c) sapply(1:3, function(f) mg14::colTrans(c,f)))
o <- c(2,4,5,6,1,3)
barplot(dNdS[,1,o], beside=TRUE, legend=TRUE, ylab="dN/dS", col=col[,c(2:5,1,6)], args.legend=list(x="topleft", bty='n')) -> b
abline(h=1)
segments(b,dNdS[,2,o],b,dNdS[,3,o])


#' # Session
#+ sessionInfo, eval=TRUE
library(devtools)
devtools::session_info()
sessionInfo()
