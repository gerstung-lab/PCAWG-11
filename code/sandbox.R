# TODO: Add comment
# 
# Author: mg14
###############################################################################


stats <- mclapply(failed, function(f){
			v <- try(readVcf(grep(f,dir(vcfPath, pattern="mnv.vcf", full.names=TRUE), value=TRUE), genome="GRCh37"))
			if(class(v)=="try-error") return("broken vcf")
			g <- grep(f,dir(vcfPath, pattern="TNC.vcf.bgz", full.names=TRUE), value=TRUE)
			if(length(g)==0) return("no TNC")
			w <- try(readVcf(g, genome="GRCh37"))
			if(class(w)=="try-error"){
				unlink(g)
				return("broken TNC, deleted")
			}
			if(nrow(v)!= nrow(w)){
				cat(nrow(v), nrow(w), g, "\n", sep="\t")
				unlink(g)
				return("broken TNC, deleted")
			}
			return("ok")
		}, mc.cores=8)


for(f in failed)
	for(ff in grep(f,dir(vcfPath, full.names=TRUE), value=TRUE))
		print(ff)

library(survival)
s <- Surv(icgc_info$donor_survival_time[match(ifelse(is.na(pcawg_info$tcga_id), as.character(pcawg_info$submitter_donor_id), pcawg_info$tcga_id), icgc_info$submitted_donor_id)],
		icgc_info$donor_vital_status[match(ifelse(is.na(pcawg_info$tcga_id), as.character(pcawg_info$submitter_donor_id), pcawg_info$tcga_id), icgc_info$submitted_donor_id)]=="deceased")
s <- s[match(sampleIds, pcawg_info$sanger_variant_calling_file_name_prefix),]

w <- !is.na(s[,1])
tt <- droplevels( tumourType[w])
nSub <- sapply(allClusters[sampleIds], nrow)-1
G <- t(apply(allGenotypes, c(1,4), sum))
c <- coxph(s[w] ~ strata(tt) + log(age[w]) + nSub[w])
summary(c)

T <- model.matrix(~tt-1)
Z <- cbind(T, G[w,])
d <- CoxHD::CoxRFX(Z, s[w], groups=c(rep("Tumour Type", ncol(T)), rep("Genes", ncol(G))), nu=0.1, which.mu=c("Genes"))


#' SNMF

snmf <- new.env()
load("../ref/snmf.RData", envir=snmf)

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

#'




s <- t(apply(sigDecomp[,,1:3], 1:2, sum))
s <- t(apply(sigDecomp30[,,1:3], 1:2, sum))

t <- factor(paste(tumourType))
p <- o <- rep(mean(age, na.rm=TRUE), length(age))
for(i in 1:100){
	l <- lm(age/o ~ s -1)
	p[-l$na.action] <- pmax(0.01,predict(l))
	m <- glm(age/p ~ t + log(1/avgWeightTrunk), family=gaussian(link="log"))
	o[-m$na.action] <- exp(predict(m))
}

tmp <- sapply(unique(as.character(t[!is.na(age)])), function(tt){
			lms <- summary(lm(age*avgWeightTrunk ~ . -1, data=as.data.frame(s)[,sigActivityHere[,tt]==1], subset=t==tt))
			#lms$r.squared
			lms$coefficients[c('Signature.1','Signature.5'),]
		}, simplify='array')

sapply(unique(as.character(t[!is.na(age)])), function(tt){
			c(rsq=summary(lm(age ~ . -1, data=as.data.frame(s)[,sigActivityHere[,tt]==1]))$r.squared,
			rsq.cn=summary(lm(age ~ . -1, data=as.data.frame(s/avgWeightTrunk)[,sigActivityHere[,tt]==1]))$r.squared)
		}, simplify='array')

m <- sapply(split(asum(s,2), t), mean)
l <- glm(age ~ log(m[t]) + log(1/avgWeightTrunk), family=gaussian(link="log"))
slope <- age
slope[-l$na.action] <- exp(predict(l))
	
summary(lm(age/slope ~ s))

f <- lm(age ~ s[,5] + tumourType)
plot(age[-f$na.action], f$fitted.values)

library(glmnet)

x <- cbind(avgWeightTrunk, s, model.matrix(~factor(paste(tumourType)) -1))
poorMansImpute <- function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}
x <- apply(x, 2, poorMansImpute)

g <- cv.glmnet(y=na.omit(age), x=x[!is.na(age),], standardize=TRUE, alpha=0.5)

#' HDP signatures
library(hdp)

t <- factor(paste(tumourType))

ppindex <- c(0, rep(1, nlevels(t)), as.numeric(t)+1, rep(as.numeric(t) + nlevels(t) + 1, each=3))
cpindex <- ppindex +1
hdp <- hdp_init(ppindex, 
		cpindex, 
		hh=rep(1, 96), 
		alphaa=rep(1, length(unique(ppindex))), 
		alphab=rep(1, length(unique(ppindex))))


#' CN TNC

timeCnTnc <- function(ID, cnids = 1:length(allBB[[ID]]), signature=5){
			i <- signature
			# gains only
			bb <- allBB[[ID]]
			vcf <- allVcf[[ID]]
			tnc <- tncToPyrimidine(vcf)
			c <- classifyMutations(vcf)
			t(sapply(cnids, function(cnid){
						if(bb$copy_number[cnid]==4 & bb$minor_cn[cnid]==2){
							w <- which(vcf %over% bb[cnid]  & c!="subclonal")
							pre <- info(vcf)$MCN[w] == 2 &  info(vcf)$TCN[w] == 4
							post <- info(vcf)$MCN[w] == 1 &  info(vcf)$TCN[w] == 4
							npre <- sum(tncProb[i,ID,1,tnc[w][pre]])
							npost <- sum(tncProb[i,ID,2,tnc[w][post]])
							return(c(npre=npre,npost=npost, t=(2*npost + .5)/(2*npost + npre+1)))
						}			
						if(bb$copy_number[cnid]==3 & bb$minor_cn[cnid]==1){
							w <- which(vcf %over% bb[cnid]  & c!="subclonal")
							pre <- info(vcf)$MCN[w] == 2
							post <- info(vcf)$MCN[w] == 1
							npre <- sum(tncProb[i,ID,1,tnc[w][pre]])
							npost <- sum(tncProb[i,ID,2,tnc[w][post]])
							return(c(npre=npre,npost=npost,(3*npost + .5)/(2*npost + 3 * npre+1)))
						}else if(bb$copy_number[cnid]==4){
							return(rep(NA,3))
						}else if(bb$copy_number[cnid]==2 & bb$major_cn[cnid]==1){
							return(rep(NA,3))
						}else if(bb$copy_number[cnid]==2 & bb$major_cn[cnid]==2){
							w <- which(vcf %over% bb[cnid]  & c!="subclonal")
							pre <- info(vcf)$MCN[w] == 2
							post <- info(vcf)$MCN[w] == 1
							npre <- sum(tncProb[i,ID,1,tnc[w][pre]])
							npost <- sum(tncProb[i,ID,2,tnc[w][post]])
							return(c(npre=npre,npost=npost,(2*npost + .5)/(2*npost + 1 * npre+1)))			
						}
						else return(rep(NA,3))
					}))
		}


timeDriverCn <- function(ID){
	do.call("c",lapply(c("sub", "indel"), function(type){
						vcf <- if(type=="sub") allVcf[[ID]] else indelAnnotated[[ID]]
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


m = read.table("/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/im3/PCAWG_Evolution_and_Heterogeneity/dNdS_vs_clonality/dNdS_allgenes_allmuts/dNdScv/Annotated_substitutions.txt", header=1, sep="\t", stringsAsFactors=F)
counts = as.matrix(t(sapply(split(m$impact, m$sampleID), function(x) table(x)[c("Synonymous","Missense")])))
counts[is.na(counts)] = 0
plot(counts)
abline(a=0,b=1,col="blue")


TP53 <- do.call("rbind",mclapply(allVcf, function(v) {w <- which(unlist(info(v)$DG)=="TP53"); DataFrame(info(v)[w,], ID=rep(meta(header(v))["ID",], length(w)))}, mc.cores=8))
TP53 <- do.call("rbind",mclapply(allVcf, function(v) {w <- which(grepl("^TP53\\|",info(v)$VD)); DataFrame(info(v)[w,], ID=rep(meta(header(v))["ID",], length(w)))}, mc.cores=8))

a <- sapply(strsplit(TP53$VD,"\\|"), `[`,5)

re <- function(x, re){
	r <- regexpr(re,x)
	regmatches(x,r)
}

for(i in 1:25+10){
	klaR::triplot(info(allVcf[[i]])$PEAR, info(allVcf[[i]])$PLAT, info(allVcf[[i]])$PSUB, label=c("Early", "Late", "Subclonal"), main=substr(names(allVcf)[[i]], 1,10), pch=16, cex=1, grid=FALSE)
	z <- klaR::tritrafo(info(allVcf[[i]])$PEAR, info(allVcf[[i]])$PLAT, info(allVcf[[i]])$PSUB)
	z <- z[!apply(is.na(z),1,any),]
	contour(kde2d(z[,1], z[,2]), add=TRUE, col='red')
}


# For Peter
vcfDriver <- mclapply(allVcf, function(vcf) vcf[unlist(info(vcf)$DG) %in% cancerGenesS], mc.cores=8)
for(vcf in vcfDriver)
	writeVcf(vcf, file=paste0("../scratch/for_pvl/driver_subs/", meta(header(vcf))["ID",], ".driver.vcf"))

vcfDriverIndel <- mclapply(indelAnnotated, function(vcf) vcf[unlist(info(vcf)$DG) %in% cancerGenesS], mc.cores=8)
for(vcf in vcfDriverIndel)
	writeVcf(vcf, file=paste0("../scratch/for_pvl/driver_indel/", meta(header(vcf))["ID",], ".driver.vcf"))


# For Phil
t <- factor(sampleInfo$sample_type[match(sampleIds, sampleInfo$tumour_id)])
s <- asum(allGenotypes[cancerGenesS,,,,], c(2:3))
g <- sapply(levels(t), function(x) asum(s[,,which(t==x)], 3)/length(which(t==x)), simplify='array')
projects <- read.table("../ref/ICGC-projects.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)
p <- projects$V2; names(p) <- projects$V1 
names(dimnames(g)) <- c("Gene","Homozygous","Project")
dimnames(g)[[3]] <- p[dimnames(g)[[3]]]
wb <- xlsx::createWorkbook()
for(h in dimnames(g)[[1]]){
	sheet <- xlsx::createSheet(wb, sheetName=h)
xlsx::addDataFrame(as.data.frame(round(t(g[h,,])*100,1)), sheet=sheet)
}
xlsx::saveWorkbook(wb, file="../scratch/pancan_driver.xlsx")


#' Poor man's dNdS
nsTable <- simplify2array(mclapply(allVcf[1:8], function(vcf) {
					subset <- grepl(paste(paste0("^",cancerGenesS,"\\|"), collapse="|"), info(vcf)$VD)
					v <- vcf[subset]
					table(factor(testDriver(v), levels=c(TRUE,FALSE), labels=c("N","S")), classifyMutations(v), tncToPyrimidine(v))}, mc.cores=8))

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
entrez <- select(org.Hs.eg.db, cancerGenesS, "ENTREZID", "SYMBOL")
txids <- select(TxDb.Hsapiens.UCSC.hg19.knownGene, na.omit(entrez$ENTREZID), "TXID", "GENEID")$TXID
cdsByGene <- cdsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")[na.omit(entrez[,2])]
cdsByTx <- cdsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "tx")[txids]

seqlevels(cdsByGene) <- sub("chr","",seqlevels(cdsByGene))
w <- sapply(cdsByGene, function(x) {s <- seqnames(x) ; length(unique(s))==1 & all(s %in% seqnames(ref))})
cds <- sapply(cdsByGene[w], function(x) unlist(scanFa(file=ref, x)))

cdsLength <- sapply(cdsByGene, function(x) sum(width(x)))
names(cdsLength) <- entrez$SYMBOL[match(names(cdsLength), entrez$ENTREZID)]

n <- t/rep(cdsLength[match(cancerGenesS, names(cdsLength))], each=5)

barplot(n[,order(-colSums(n))[1:100]], col=RColorBrewer::brewer.pal(5,"Set1"), legend=TRUE ,las=2)


#' ## Graphs and timing
#' Posets
allPosets <- sapply(sampleIds, function(ID){
			reduceToCoverRelations(applyPigeonHole(ID))
		})

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


#' Plot all posets
#+ allGraphs, fig.width=2, fig.height=2
par(mar=c(1,1,1,1))
i <- 1; for(g in allGraphs){
	plotPoset(g);
	title(main=names(allGraphs)[i])
	i<-i+1
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

genotypes <- sapply(allVcf, function(vcf) table(factor(unlist(info(vcf)$DG), levels=as.character(mutsigDrivers)), factor(info(vcf)$CLS, levels=c("early","late","subclonal"))), simplify='array')

genotypes <- mclapply(allVcf, function(vcf) sapply(mutsigDrivers, function(g) {w <- which(unlist(info(vcf)$DG)==g); if(length(w)==0) 0 else as.numeric(info(vcf)$CLS[w[1]])}), mc.cores=8)
genotypes <- simplify2array(genotypes)
genotypes <- do.call("data.frame", sapply(recMutations, function(g) factor(genotypes[g,], labels=c("none","early","late","subclonal"), levels=0:3), simplify=FALSE))

genotypes <- simplify2array(mclapply(allVcf, function(vcf) table(factor(info(vcf)$DG, levels=mutsigDrivers), factor(info(vcf)$CLS, levels=c("early","late","subclonal"))), mc.cores=8))


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


#' ??
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
						vcf <- if(type=="sub") allVcf[[ID]] else indelAnnotated[[ID]]
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
