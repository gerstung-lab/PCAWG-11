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
G <- t1(apply(allGenotypes, c(1,4), sum))
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
b <- barplot(t1(nmfFit$P), las=2, col=brewer.pal(8,"Set1"), border=NA,names.arg=rep("", nrow(sigData)))
mtext(at=b,rownames(sigData), side=1, las=2, cex=.7, col=brewer.pal(6,"Set1")[f])

cor(snmf$snmfFit$P[-1,], nmfFit$P[,])
newSig <- nmfFit$P
w <- which.max(cor(snmf$snmfFit$P[-1,1], nmfFit$P[,]))
newSig <- newSig[,c(w, setdiff(1:n,w))]

realTimeBranch <- sapply(sigList, function(x) predictRealTime(t1(x), newSig))
realTimeAll <- predictRealTime(sigData, newSig)

fit <- lm(log(age) ~ log(realTimeAll) + tumourType + log(avgPower) + log(avgWeightTrunk))
summary(fit)

#'




s <- t1(apply(sigDecomp[,,1:3], 1:2, sum))
s <- t1(apply(sigDecomp30[,,1:3], 1:2, sum))

t1 <- factor(paste(tumourType))
p <- o <- rep(mean(age, na.rm=TRUE), length(age))
for(i in 1:100){
	l <- lm(age/o ~ s -1)
	p[-l$na.action] <- pmax(0.01,predict(l))
	m <- glm(age/p ~ t1 + log(1/avgWeightTrunk), family=gaussian(link="log"))
	o[-m$na.action] <- exp(predict(m))
}

tmp <- sapply(unique(as.character(t1[!is.na(age)])), function(tt){
			lms <- summary(lm(age*avgWeightTrunk ~ . -1, data=as.data.frame(s)[,sigActivityHere[,tt]==1], subset=t1==tt))
			#lms$r.squared
			lms$coefficients[c('Signature.1','Signature.5'),]
		}, simplify='array')

sapply(unique(as.character(t1[!is.na(age)])), function(tt){
			c(rsq=summary(lm(age ~ . -1, data=as.data.frame(s)[,sigActivityHere[,tt]==1]))$r.squared,
			rsq.cn=summary(lm(age ~ . -1, data=as.data.frame(s/avgWeightTrunk)[,sigActivityHere[,tt]==1]))$r.squared)
		}, simplify='array')

m <- sapply(split(asum(s,2), t1), mean)
l <- glm(age ~ log(m[t1]) + log(1/avgWeightTrunk), family=gaussian(link="log"))
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

t1 <- factor(paste(tumourType))

ppindex <- c(0, rep(1, nlevels(t1)), as.numeric(t1)+1, rep(as.numeric(t1) + nlevels(t1) + 1, each=3))
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
			t1(sapply(cnids, function(cnid){
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
counts = as.matrix(t1(sapply(split(m$impact, m$sampleID), function(x) table(x)[c("Synonymous","Missense")])))
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
t1 <- factor(sampleInfo$sample_type[match(sampleIds, sampleInfo$tumour_id)])
s <- asum(allGenotypes[cancerGenesS,,,,], c(2:3))
g <- sapply(levels(t1), function(x) asum(s[,,which(t1==x)], 3)/length(which(t1==x)), simplify='array')
projects <- read.table("../ref/ICGC-projects.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)
p <- projects$V2; names(p) <- projects$V1 
names(dimnames(g)) <- c("Gene","Homozygous","Project")
dimnames(g)[[3]] <- p[dimnames(g)[[3]]]
wb <- xlsx::createWorkbook()
for(h in dimnames(g)[[1]]){
	sheet <- xlsx::createSheet(wb, sheetName=h)
xlsx::addDataFrame(as.data.frame(round(t1(g[h,,])*100,1)), sheet=sheet)
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

n <- t1/rep(cdsLength[match(cancerGenesS, names(cdsLength))], each=5)

barplot(n[,order(-colSums(n))[1:100]], col=RColorBrewer::brewer.pal(5,"Set1"), legend=TRUE ,las=2)





#' #### Model based
wgdWeight <- sapply(sampleIds, function(ID){
			if(is.na(wgdDeam[2,ID])) return(NA)
			bb <- allBB[[ID]]
			ix <- abs(bb$clonal_frequency - purityPloidy[ID, "purity"]) < 0.01 & bb$minor_cn==2 & bb$copy_number==4
			sum(as.numeric(width(bb[ix]))) / 3e9
		})

glm(nDeamTrunk ~ offset(log(1/avgWeightTrunk*wgdDeam[2,]/wgdWeight)), subset=wgdWeight>0, family="poisson")

plot(nDeamTrunk*avgWeightTrunk, wgdDeam[2,]/wgdWeight); abline(0,1)
quantile((wgdDeam[2,]/wgdWeight)/(nDeamTrunk*avgWeightTrunk), na.rm=TRUE)

fit <- glm(nDeam ~ log(age) + offset(log(1/avgWeightTrunk)) + tumourType -1, family="poisson")
plot(sort(wgdDeam[2,])/exp(coef(fit)[1])); abline(0,1)

p <- predict(fit, newdata=data.frame(age=1/exp(1), avgWeightTrunk=avgWeightTrunk, tumourType=tumourType)[-fit$na.action,], type='response', se.fit=TRUE)
plot(age[-fit$na.action], (wgdDeam[2,]/wgdWeight)[-fit$na.action] / p$fit, ylim=c(0,100), xlab="Age at diagnosis", ylab="Infered age during WGD"); abline(0,1)
e <- (t1(sapply(wgdDeam[2,], function(w) qpois(c(0.025,0.975), w)))/wgdWeight)[-fit$na.action,] / p$fit
segments(age[-fit$na.action], e[,1], age[-fit$na.action], e[,2])


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
boxplot(t1(allDistMutations[1,-1,][order(c),]), las=2, horizontal=TRUE)

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


summary(lm((1-tWgd)*age~ tumourType + purityPloidy[sampleIds,1],))

w <- setdiff(sampleIds, names(allVcf))
tmp <- mclapply(w, function(ID){
			try(readVcf(grep(ID, dir(vcfPath, pattern=".complete_annotation.vcf.bgz", full.names=TRUE), value=TRUE), genome="GRCh37"))
		}, mc.cores=MCCORES)
names(tmp) <- w
table(sapply(tmp, class))
allVcf <- c(allVcf, tmp)
allVcf <- allVcf[sampleIds]
save(allVcf, file="allVcf.RData")

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


timeAcn <- function(vcf, bb=allBB[[meta(header(vcf))['ID',1]]], subset=info(vcf)$CLS != "subclonal"){
	v <- vcf[which(subset)]
	cnid <- sapply(info(v)$CNID, paste, collapse=",")
	n <- table(cnid)
	tcn <- info(v)$TCN
	if(is.null(tcn))
		tcn <- (info(v)$MJCN + info(v)$MNCN)
	cnw <- info(v)$MCN / tcn
	res <- aggregate(cnw, list(cnid), function(x,...) 1/mean(x,...), na.rm=TRUE)
	colnames(res) <- c("cnid", "avg.ploidy")
	m <- match(seq_along(bb), as.numeric(as.character(res$cnid)))
	T0 <- (bb$major_cn!=0) + (bb$minor_cn!=0)
	res <- data.frame(res[m,], MJCN=bb$major_cn, MNCN=bb$minor_cn, time=(bb$copy_number-res$avg.ploidy[m])/(bb$copy_number-T0), n=n[m])
	return(res)
}


bbTimeDeam <- mclapply(allVcf, function(vcf) timeAcn(vcf, subset=isDeamination(vcf)), mc.cores=MCCORES)
td <- sapply(bbTimeDeam, function(x) weighted.mean(x[x$MJCN > 1,"time"],x[x$MJCN > 1,"n"], na.rm=TRUE))
plot(tWgd, td)


#' All sig, all clonal
sigTableClonal <- asum(sigTable[,,1:3],3)
w <- colSums(sigTableClonal) != 0

M <- matrix(0,nrow=30, ncol=ncol(sigTableClonal))
rownames(M) <- colnames(S)


for(t1 in levels(tumourType)){
	v <- tumourType==t1 
	E <- nmSolve(matrix(sigTableClonal, nrow=96)[,which(w & v)], S[,sigActivityHere[,t1]==1], maxIter=10000, tol=1e-4)
	M[rownames(E), which(v&w)] <- E
}

sigDecompClonal <- array(M, dim=c(30, ncol(sigTableClonal)))
rm(M)

tncProbClonal <- sapply(1:96, function(i) S[i,] * sigDecompClonal, simplify='array')
s <- apply(tncProbClonal, 2:3, sum)
tncProbClonal <- tncProbClonal / rep(s, each=30)
rm(s)
dimnames(tncProbClonal)[[1]] <- colnames(S)
dimnames(tncProbClonal)[[3]] <- rownames(S)
dimnames(tncProbClonal)[[2]] <- colnames(sigTableClonal)


#' PCAWG signatures
library(xlsx)
PCAWG_signatures <- read.table("../ref/PCAWG_signature_patterns.txt", header=TRUE, sep="\t")
PCAWG_S <- as.matrix(PCAWG_signatures[,-2:-1])
rownames(PCAWG_S) <- paste0(substr(PCAWG_signatures$Mutation.Subtype,1,1), "[", PCAWG_signatures$Mutation.Type, "]",substr(PCAWG_signatures$Mutation.Subtype,3,3))
PCAWG_S <- PCAWG_S[match(rownames(sigTable),rownames(PCAWG_S)),]
PCAWG_activity <- read.table("../ref/PCAWG_signatures_in_samples.txt", header=TRUE, sep="\t")
PCAWG_A <- as.matrix(PCAWG_activity[,-2:-1])
rownames(PCAWG_A) <- PCAWG_activity[,2]


PCAWG_sigDecomp <- array(0, dim=c(dim(sigTable)[-1], dim(PCAWG_S)[2]), dimnames=c(dimnames(sigTable)[-1], list(colnames(PCAWG_S))))

for(sample in rownames(PCAWG_sigDecomp)){
	if(class(try({
				w <- PCAWG_A[sample,]!=0
				E <- nmSolve(sigTable[,sample,], PCAWG_S[,w], maxIter=10000, tol=1e-4)
				PCAWG_sigDecomp[sample,,w] <- t(E)
			})) == "try-error") PCAWG_sigDecomp[sample,,] <- NA
}

t1 <- asum(PCAWG_sigDecomp, 2)
m <- match(rownames(t1), rownames(PCAWG_A))
plot(t1[,2]+0.1, PCAWG_A[m,2], log='xy')


plot(t1[,1]+0.1, PCAWG_A[m,1], log='xy')

sapply(1:ncol(t1), function(i)  cor(t1[,i]+0.1, PCAWG_A[m,i], use='c')) 
par(mfrow=c(7,8), mar=c(2,2,2,1), mgp=c(2,.5,0))
for(i in 1:ncol(t1)) 
	plot(round(t1[,i]), PCAWG_A[m,i], log='xy', main=colnames(t1)[i])

par(mfrow=c(1,1), mar=c(3,3,1,1))
scatterpie(x=t1[,1], y=PCAWG_A[m,1], prop=PCAWG_sigDecomp[,,1], xlim=c(0,1e3), ylim=c(0,1e3), radius=10, col=set1[c(1,2,4,3)], xlab="Moritz", ylab="Ludmil")

plot(rowSums(t1), rowSums(PCAWG_A[m,]))

PCAWG_sigDecompClonal <- t1(asum(PCAWG_sigDecomp[,1:3,],2)) 


#' PCAWG TNC prob
PCAWG_tncProbClonal <- sapply(1:96, function(i) PCAWG_S[i,] * PCAWG_sigDecompClonal, simplify='array')
s <- apply(PCAWG_tncProbClonal, 2:3, sum)
PCAWG_tncProbClonal <- PCAWG_tncProbClonal / rep(s, each=nrow(PCAWG_tncProbClonal))
rm(s)
dimnames(PCAWG_tncProbClonal)[[1]] <- colnames(PCAWG_S)
dimnames(PCAWG_tncProbClonal)[[3]] <- rownames(PCAWG_S)
dimnames(PCAWG_tncProbClonal)[[2]] <- colnames(sigTableClonal)

#' Merging signatures
m <- sub("Signature.(PCAWG.)*","",sub("[a-z]$","",colnames(PCAWG_S)))
s <- sub(".+[0-9]+(?=[a-z]*$)","", colnames(PCAWG_S), perl=TRUE)
n <- sapply(unique(m), function(mm) paste0(mm, paste(s[m==mm], collapse=",")))
PCAWG_sigDecomp_merged <- sapply(unique(m), function(mm) asum(PCAWG_sigDecomp[,,m==mm, drop=FALSE],3), simplify='array')
dimnames(PCAWG_sigDecomp_merged)[[3]] <- n

#' PCAWG TNC prob (merged)
PCAWG_tncProbClonal_merged <- sapply(unique(m), function(mm) asum(PCAWG_tncProbClonal[m==mm,,, drop=FALSE],1), simplify='array')
dimnames(PCAWG_tncProbClonal_merged)[[3]] <- sub("Signature[.PCAWG]+","S",dimnames(PCAWG_tncProbClonal_merged)[[3]])
PCAWG_tncProbClonal_merged <- aperm(PCAWG_tncProbClonal_merged, c(3,1,2))

#' Variances

i <- intersect(w, which(newTumourTypes=="Kidney-CCRCC"))[1]

#sigTable[,sample,], PCAWG_S[,w]

V <- solve(poisI(PCAWG_S[,PCAWG_A[i,]!=0], PCAWG_sigDecomp[i,1,PCAWG_A[i,]!=0],sigTable[,i,1]))
U <- poisdl(PCAWG_S[,PCAWG_A[i,]!=0], PCAWG_sigDecomp[i,1,PCAWG_A[i,]!=0],sigTable[,i,1])

#' Timing using ACN
timeAcnTnc <- function(vcf, bb=allBB[[meta(header(vcf))['ID',1]]], subset=info(vcf)$CLS != "subclonal", select = c("all","pq"), tncProb){
	select <- match.arg(select)
	ID <- meta(header(vcf))['ID',1]
	v <- vcf[which(subset)]
	cnid <- as.character(as.list(info(v)$CNID))#sapply(info(v)$CNID, paste, collapse=",")
	if(select=='all'){
		seg <- cnid
		bbSeg <- as.character(seq_along(bb))
	}else{
		o <- findOverlaps(bb, pq, select='first')
		bbSeg <- paste(names(pq)[o], bb$major_cn,bb$minor_cn, round(bb$clonal_frequency*2,1)/2, sep=":")
		seg <- bbSeg[as.numeric(cnid)]
	}
	tnc <- tncToPyrimidine(v)
	n <- table(seg)
	tcn <- info(v)$TCN
	if(is.null(tcn))
		tcn <- (info(v)$MJCN + info(v)$MNCN)
	cnw <- info(v)$MCN / tcn
	cnw30 <- cnw * t1(tncProb[,ID,tnc])
	res <- aggregate(cnw30, list(seg), function(x,...) sum(x,...), na.rm=TRUE)
	res0 <- aggregate(t1(tncProb[,ID,])[tnc,], list(seg), function(x,...) sum(x,...), na.rm=TRUE)
	acn <- data.frame(seg=res0[,1],res0[,-1]/res[,-1])
	w <- !duplicated(bbSeg)
	m <- match(bbSeg[w], as.character(acn[,1]))
	TCN0 <- (bb$major_cn[w]!=0) + (bb$minor_cn[w]!=0)
	time <- (bb$copy_number[w]-acn[m,-1])/(bb$copy_number[w]-TCN0)
	time[time > 1] <- 1
	acn <- data.frame(cnid=acn[m,1], avg.ploidy=acn[m,-1], MJCN=bb$major_cn[w], MNCN=bb$minor_cn[w], time=time, n=n[acn[m,1]])[!is.na(m),]
	return(acn)
}

#' All signatures
bbTimeTnc <- mclapply(allVcf, function(vcf) timeAcnTnc(vcf,tncProb=PCAWG_tncProbClonal_merged), mc.cores=MCCORES)
bbTimeTncAll <- do.call("rbind", bbTimeTnc)
bbTimeTncAll$sampleId <- factor(sub("\\..+","",rownames(bbTimeTncAll)))

tt <- unlist(sapply(names(bbTimeTnc), function(x) rep(x, nrow(bbTimeTnc[[x]]))))
o <- rank(aggregate(bbTimeTncAll$time.Signature.1, list(tumourType[tt]), median, na.rm=TRUE)[,2])
boxplot(bbTimeTncAll$time.Signature.1 ~ tumourType[tt], subset=bbTimeTncAll$n>10 & bbTimeTncAll$time.Signature.1<=1, at=o, ylim=c(0,1), col=col)

#' Merged signatures


plot(bbTimeTnc[[1]]$time.Signature.3 ~ bbTimeTnc[[1]]$time.Signature.5)

par(las=2)
names(age) <- sampleIds
names(isWgd) <- sampleIds
w <- bbTimeTncAll$n>100 & bbTimeTncAll$time.Signature.1<=1
o <- rank(aggregate(((1-bbTimeTncAll$time.Signature.1) * age[tt])[w], list(tumourType[tt][w]), median, na.rm=TRUE)[,2])
par(mfrow=c(3,1))
boxplot((1-bbTimeTncAll$time.Signature.1) * age[tt] ~ tumourType[tt], subset=w, at=o, ylim=c(0,75), col=rep(col,2)[o], ylab="Time lag", main="all samples")
boxplot((1-bbTimeTncAll$time.Signature.1) * age[tt] ~ tumourType[tt], subset=w & !isWgd[tt] , at=o, ylim=c(0,75), col=rep(col,2)[o], ylab="Time lag", main="no WGD")
boxplot((1-bbTimeTncAll$time.Signature.1) * age[tt] ~ tumourType[tt], subset=w & isWgd[tt] , at=o, ylim=c(0,75), col=rep(col,2)[o], ylab="Time lag", main="with WGD")

par(mfrow=c(1,1))
o <- rank(aggregate(bbTimeTncAll$time.Signature.1[w], list(tumourType[tt][w]), median, na.rm=TRUE)[,2])
boxplot(bbTimeTncAll$time.Signature.1 ~  tumourType[tt], subset=w, at=o, ylim=c(0,1), col=rep(col,2)[o], ylab="Relative time (sig 1)", las=2)


pdf("sigTime_merged.pdf",25,25,pointsize=8)
pairs(bbTimeTncAll[bbTimeTncAll$n > 100, grep("time", colnames(bbTimeTncAll))], xlim=c(0,1), ylim=c(0,1), col="#00000022", pch=16, upper.panel=NULL,lower.panel=function(x,y,...){
			points(x,y,...)
			z <- cbind(x,y) %*% R(pi/4)
			zx <- z[,1]
			zy <- z[,2]
			t <- seq(0,sqrt(2),l=100)
			try({
						#contour(kde2d(a[w], y[w], n=200), col=mg14::colTrans(rep(col,2)[i]), main=sub("time.","",x), ylim=c(-20,0), xlim=c(-20,0), xlab="Time before diagnosis (age x sig 1)", ylab="Time before diagnosis (age x sig x)")
						p <- predict(loess(zy ~ zx, data=bbTimeTncAll, span=0.5), newdata=data.frame(zx=t))
						f <- (cbind(t,p) %*% R(-pi/4))
						abline(0,1, col='blue')
						lines(f[,1], f[,2], col='#FF0000', lwd=2)
				}, silent=TRUE)
		}, text.panel=function(x=0.5,y=0.5, txt, ...){
			text(x,y,sub('time.','',sub('Signature(.PCAWG.)*','S.',txt)), ...)
		}, cex.labels=3
)
dev.off()

s <- t1(PCAWG_sigDecompClonal)
s[s<=1] <- NA
sr <- s[,-1]/(s[,1]+1)
par(mar=c(3,3,1,3), mgp=c(2,.5,0), bty="L", las=0)
plot(NA, NA, xlab="Relative time (sig 1)", xlim=c(0,1), ylim=c(1e-1,1e3), log='y', ylab="# Signature x mutations/signature 1 mutation")
i <- 0
for(n in colnames(sr)[1:20])try({
	x <- bbTimeTncAll$time.Signature.PCAWG.1
	y <- bbTimeTncAll[[paste0("time.",n)]]
	w <- bbTimeTncAll$n > 100 & !is.na(x) &! is.na(y)
	z <- (cbind(x[w],y[w]) %*% R(pi/4))
	#plot(x[w],y[w])
	fit <- 	cobs(z[,1],z[,2], w=1/sqrt(bbTimeTncAll$n[w]), degree=1, lambda=0.1, nknots=3, pointwise = matrix(c(0,0,0,sqrt(2),0,0), ncol=3, byrow=TRUE))
	p <- predict(fit, newdata=data.frame(x=seq(0,sqrt(2),0.01)))
	sp <- p %*% R(-pi/4)
	#lines(zp, col='red')
	ds <- diff(sp[,2]) / diff(sp[,1])
	m <- median(sr[,n], na.rm=TRUE)
	print(length(fit$knots))
	lines(sp[-1,1] + diff(sp[,1]), ds * m, col=c(set1[1],"black",set1[2])[1+ (ds[1]==ds[99]) + 2*(ds[99] > ds[1])])
	if(abs(log10(ds[1]/ds[99]))>0.2)
		mtext(sub("Signature.(PCAWG.)*","S",n), side=4, at=ds[99] * m, las=2)
	#text(i,labels=sub("Signature.(PCAWG.)*","S",n),  y=ds[1] * m, pos=4, cex=0.66)
	#i <- i+0.01
})
	

#' Getting chromosome arms
library(biovizBase)
data(hg19IdeogramCyto, package = "biovizBase")
cytoCol = getOption("biovizBase")$cytobandColor
cyto <- hg19IdeogramCyto
cyto <- sort(cyto)
seqlevels(cyto) <- sub("chr","", seqlevels(cyto))
chrCyto <- split(cyto, seqnames(cyto))
detach(package:biovizBase)

pq <- paste0(seqnames(unlist(chrCyto)), sub("[0-9].+","", unlist(chrCyto)$name))
mcols(chrCyto)$pq <- pq
pq <- unlist(GRangesList(sapply(split(unlist(chrCyto), pq), reduce)))

allRanges <- GRangesList(sapply(names(bbTimeTnc), function(ID) granges(allBB[[ID]])[bbTimeTnc[[ID]]$cnid]))

chr <- sapply(names(bbTimeTnc), function(ID) as.character(seqnames(allBB[[ID]])[bbTimeTnc[[ID]]$cnid]))
boxplot(bbTimeTncAll$time.Signature.1 ~ unlist(chr))

sigTimeLoess <- sapply(1:30, function(i){
			t <- try(predict(loess(as.formula(paste0("time.Signature.",i," ~ time.Signature.5")), data=bbTimeTncAll, subset=bbTimeTncAll$n>10 & bbTimeTncAll$time.Signature.1<=1 , span=0.5), newdata=data.frame(time.Signature.5=seq(0,1,l=1000))))
			if(class(t)!='try-error')
				t
			else(rep(NA, 1000))
		})

x <- seq(0,1,l=1000)
col <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"))
plot(x, sigTimeLoess[,1], type='l', col=col[1], xlab="Time (sig 5)", ylab="Time (sig x)")

for(i in 2:30)
	lines(x, sigTimeLoess[,i], col=rep(col,2)[i], lty=ifelse(i<=15, 1,2))
legend("topleft", paste0("Signature ", 1:30), col=rep(col,2)[1:30],  bty='n', cex=.7, ncol=2, lty=rep(c(1,2), each=15))
abline(0,1)



R <- function(theta) matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol=2)

#' 2d clustering
par(mfrow=c(5,6), mar=c(3,3,3,1), mgp=c(2,.5,0))
n <- table(bbTimeTncAll$sampleId)
for(i in (1:30)[-1]){
	x <- paste0("time.Signature.",i)
	w <- ! (is.na(bbTimeTncAll[[x]]) | is.na(bbTimeTncAll$time.Signature.1) | is.infinite(bbTimeTncAll[[x]])) & bbTimeTncAll$n > 10 #& bbTimeTncAll$MJCN <= 2 & bbTimeTncAll$MNCN <= 2 
	t1 <- seq(0,1,l=100)
	try({
				contour(kde2d(bbTimeTncAll$time.Signature.1[w], bbTimeTncAll[[x]][w] ), col=mg14::colTrans(rep(col,2)[i]),  ylim=c(0,1), xlim=c(0,1),  main=sub("time.","",x), xlab="Time (sig 1)", ylab="Time (sig x)")
				#p <- predict(loess(as.formula(paste0(x, "~ time.Signature.1")), data=bbTimeTncAll, subset=w , span=0.5), newdata=data.frame(time.Signature.1=t), se=TRUE)
				a <- bbTimeTncAll$time.Signature.1[w]
				b <- bbTimeTncAll[[x]][w]
				z <- cbind(a,b) %*% R(pi/4)
				c <- cobs(z[,1], z[,2],  w=1/n[bbTimeTncAll$sampleId[w]], degree=2, lambda=0.01) #pointwise = matrix(c(0,0,0,sqrt(2),0,0), ncol=3)
				p <- as.data.frame(predict(c, z=seq(0,sqrt(2),l=100)))
				#polygon( c(t, rev(t)), c(p$fit+ 1.96*p$se.fit, rev(p$fit - 1.96*p$se.fit)), col=paste0(rep(col,2)[i], "88"), border=NA)
				#f <- p$fit
				f <- (as.matrix(p[,1:2]) %*% R(-pi/4))[,2]
				lines(t1, f, col=rep(col,2)[i], lwd=2)
				lines(t1[-1], diff(f)/diff(t1)/5, col=rep(col,2)[i], lwd=2)
				abline(h=0)
			}, silent=FALSE)
}

#' 2d clustering - real time
par(mfrow=c(5,6), mar=c(3,3,3,1), mgp=c(2,.5,0))
for(i in (1:30)[-1]){
	x <- paste0("time.Signature.",i)
	a <- (bbTimeTncAll$time.Signature.1-1) * age[bbTimeTncAll$sampleId]
	b <- (bbTimeTncAll[[x]]-1) * age[bbTimeTncAll$sampleId]
	
	w <- ! (is.na(a) | is.na(b) | is.infinite(b)) & bbTimeTncAll$n > 100 & !is.na(a)
	t1 <- seq(-50,0,l=100)
	try({
				contour(kde2d(a[w], b[w]), col=mg14::colTrans(rep(col,2)[i]), main=sub("time.","",x), ylim=c(-50,0), xlim=c(-50,0), xlab="Time before diagnosis (age x sig 1)", ylab="Delta")
				c <- cobs(a[w], b[w], constraint="increase", pointwise = matrix(c(0,0,0,-50,0,-50), ncol=3), w=1/n[bbTimeTncAll$sampleId[w]], degree=2, lambda=0.5, nknots=10)
				p <- as.data.frame(predict(c, z=t1, interval="confidence"))
				#polygon( c(t, rev(t)), c(p$cb.up, rev(p$cb.lo)), col=paste0(rep(col,2)[i], "88"), border=NA)
				lines(t1, p$fit, col=rep(col,2)[i], lwd=2)
				lines(t1[-1], diff(p$fit)/diff(t1)*10 - 50, col=rep(col,2)[i], lwd=2)
				
				abline(0,1)
			})
}

#' 2d clustering - rlative to diagnosis
par(mfrow=c(5,6), mar=c(3,3,3,1), mgp=c(2,.5,0))
for(i in (1:30)[-1]){
	x <- paste0("time.Signature.",i)
	a <- I(-(1-bbTimeTncAll$time.Signature.1) * age[bbTimeTncAll$sampleId])
	y <- log10((1-bbTimeTncAll[[x]] )/(1-bbTimeTncAll$time.Signature.1))
	w <- ! (is.na(bbTimeTncAll[[x]]) | is.na(bbTimeTncAll$time.Signature.1) | is.infinite(bbTimeTncAll[[x]])) & bbTimeTncAll$n > 50 & !is.na(a) & !is.na(y) &!is.infinite(y)
	t1 <- seq(-50,0,l=100)
	try({
				contour(kde2d(a[w], y[w], n=100), col=mg14::colTrans(rep(col,2)[i]), main=sub("time.","",x), ylim=c(-1,1), xlim=c(-50,0), xlab="Time before diagnosis (age x sig 1)", ylab="log10FC")
				p <- predict(loess(as.formula(paste0("y ~ a")), data=bbTimeTncAll, subset=w , span=0.5, weights=bbTimeTncAll$n), newdata=data.frame(a=t1), se=TRUE)
				ci <- c(p$fit+ 1.96*p$se.fit, rev(p$fit - 1.96*p$se.fit))
				polygon( c(t1, rev(t1))[!is.na(ci)], na.omit(ci), col=paste0(rep(col,2)[i], "88"), border=NA)
				lines(t1, p$fit, col=rep(col,2)[i], lwd=2)
				abline(h=0)
			})
}



#' 2d clustering - rlative to diagnosis
par(mfrow=c(5,6), mar=c(3,3,3,1), mgp=c(2,.5,0))
for(i in (1:30)[-1]){
	x <- grep("time", colnames(bbTimeTncAll), value=TRUE)[i]#paste0("time.Signature.PCAWG.",i)
	a <- I(-(1-bbTimeTncAll$time.Signature.PCAWG.1) * age[bbTimeTncAll$sampleId])
	y <- I(-(1-bbTimeTncAll[[x]]) * age[bbTimeTncAll$sampleId])
	z <- cbind(a,y) %*% R(pi/4)
	zx <- z[,1]
	zy <- z[,2]
	w <- ! (is.na(bbTimeTncAll[[x]]) | is.na(bbTimeTncAll$time.Signature.PCAWG.1) | is.infinite(bbTimeTncAll[[x]])) & bbTimeTncAll$n > 50 & !is.na(a) & !is.na(y) &!is.infinite(y) & bbTimeTncAll$MJCN ==2 & bbTimeTncAll$MNCN==1 #&! (a==0 | y==0)
	t1 <- seq(-50*sqrt(2),0,l=100)
	try({
				contour(kde2d(a[w], y[w], n=200), col=mg14::colTrans(rep(col,2)[i]), main=sub("time.","",x), ylim=c(-20,0), xlim=c(-20,0), xlab="Time before diagnosis (age x sig 1)", ylab="Time before diagnosis (age x sig x)")
				p <- predict(loess(zy ~ zx, data=bbTimeTncAll, subset=w , span=0.5, weights=bbTimeTncAll$n), newdata=data.frame(zx=t1), se=TRUE)
				ci <- c(p$fit + 1.96*p$se.fit, rev(p$fit - 1.96*p$se.fit))
				polygon( cbind(c(t1, rev(t1))[!is.na(ci)], na.omit(ci)) %*% R(-pi/4), col=paste0(rep(col,2)[i], "88"), border=NA)
				lines(cbind(t1, p$fit)%*% R(-pi/4) , col=rep(col,2)[i], lwd=2)
				abline(0,1)
			})
}


bbTimeTncAvg <- as.data.frame(t1(sapply(bbTimeTnc, function(x){
			colSums(x[2:64] * x$n) / sum(x$n) 
		})))
bbTimeTncAvg$sampleId <- sampleIds

par(mfrow=c(5,6), mar=c(3,3,3,1), mgp=c(2,.5,0))
for(i in (1:30)[-1]){
	x <- paste0("time.Signature.",i)
	a <- I(bbTimeTncAvg$time.Signature.1 * age[bbTimeTncAvg$sampleId])
	w <- ! (is.na(bbTimeTncAvg[[x]]) | is.na(bbTimeTncAvg$time.Signature.1) | is.infinite(bbTimeTncAvg[[x]]))  & !is.na(a)
	t1 <- seq(0,100,l=100)
	try({
				contour(kde2d(bbTimeTncAvg$time.Signature.1[w] * age[bbTimeTncAvg$sampleId[w]], (bbTimeTncAvg$time.Signature.1[w] - bbTimeTncAvg[[x]][w] )* age[bbTimeTncAvg$sampleId[w]]), col=mg14::colTrans(rep(col,2)[i]), main=sub("time.","",x), ylim=c(-20,20), xlab="Time (age x sig 1)", ylab="Delta")
				p <- predict(loess(as.formula(paste0("(time.Signature.1 - ",x, ") * age[bbTimeTncAvg$sampleId] ~ a")), data=bbTimeTncAvg, subset=w , span=0.5), newdata=data.frame(a=t1), se=TRUE)
				ci <- c(p$fit+ 1.96*p$se.fit, rev(p$fit - 1.96*p$se.fit))
				polygon( c(t1, rev(t1))[!is.na(ci)], na.omit(ci), col=paste0(rep(col,2)[i], "88"), border=NA)
				lines(t1, p$fit, col=rep(col,2)[i], lwd=2)
				abline(h=0)
			})
}


#' PQ-timing

bbTimeTncPq <- mclapply(allVcf, function(vcf) timeAcnTnc(vcf, select='pq'), mc.cores=MCCORES)
bbTimeTncPqAll <- do.call("rbind", bbTimeTncPq)


allSigAcnCobs <- sapply(2:30, function(i){
			a <- mclapply(bbTimeTncPq, function(b) try({
									x <- b$time.Signature.1
									y <- b[,paste0("time.Signature.",i)]
									c <- cobs(x,y, constraint="increase", pointwise = matrix(c(0,0,0,1,0,1), ncol=3), w=b$n[!is.na(x) & !is.na(y)], degree=1)
									predict(c, z=seq(0,1,l=100))[,2]
								}, silent=TRUE), mc.cores=8)
			c <- sapply(a, class) != "try-error"
			b <- sapply(a[c], I, simplify='array')
		})


par(mfrow=c(5,6), mar=c(3,3,3,1), mgp=c(2,.5,0), las=1)
for(j in seq_along(allSigAcnCobs)) try({
				b <- allSigAcnCobs[[j]]
				c <- apply(b[,], 2, diff)/0.01
				q <- apply(c, 1, quantile, c(0.25, 0.5, 0.75))
				plot(seq(0,1,l=99), c[,1], type='s', ylim=c(0,2), col='lightgrey', main=paste("Signature", j+1), ylab="Relative activity", xlab="Rel time (sig 1)")
				for(i in 2:ncol(c)) lines(seq(0,1,l=99), c[,i], col='lightgrey', type='s')
				#lines(seq(0,1,l=99), rowMeans(c[,-798]))
				x <- seq(0,1,l=99)
				polygon(c(x, rev(x)), c(q[1,], rev(q[3,])), col=paste0(rep(col,2)[j], "88"), border=NA)			
				lines(x, q[2,], lty=2, col=rep(col,2)[j], lwd=2)
				#lines(seq(0,1,l=99), q[1,], lty=3, col=col[j])
				#lines(seq(0,1,l=99), q[3,], lty=3, col=col[j])
				abline(h=1, lty=3)
			})

par(mfrow=c(1,1))
ID <- sampleIds[1]
plot(bbTimeTncPq[[ID]][,'time.Signature.1'], bbTimeTncPq[[ID]][,'time.Signature.3'], xlab="Relative time (sig 1)", ylab="Relative time (sig 3)", main=ID)
b <- allSigAcnCobs[[2]]
abline(0,1,lty=3)
lines(seq(0,1,l=100), b[,ID], lty=2)
lines(x, diff(b[,ID])*20)
legend("topleft", pch=c(1, NA,NA), lty=c(NA, 2,3), legend=c("Segments", "Fit", "Activity"))



#' Early-late timing
par(mfrow=c(1,3), cex=1)
s <- aperm(PCAWG_sigDecomp_merged, c(3,1,2)) #sigDecomp
s[s==0]<-NA
or <- function(x) x[2,1]/x[1,1]*x[1,2]/x[2,2]
1/sapply(2:nrow(s), function(i) apply(s[c(1,i),,1:2]+10, 2, or)) -> oddsRatioEarlyLate
r <- order(-apply(oddsRatioEarlyLate,2, median, na.rm=TRUE))
colnames(oddsRatioEarlyLate) <- rownames(s)[-1]#paste0("sig", 2:30)
boxplot(oddsRatioEarlyLate[,r], log='x', ylim=c(1e-2,1e2), horizontal=TRUE, ylim=c(0,30), at=1:ncol(oddsRatioEarlyLate), col=col, lty=1, staplewex=0, pch=16, cex=.5, las=1, xlab="OR late/early v sig 1", notch=TRUE)
abline(v=1, lty=3)

#' Clonal-subclonal timing
s <- abind::abind(asum(s[,,1:3], 3), s[,,4], along=3)
1/sapply(2:nrow(s), function(i) apply(s[c(1,i),,1:2]+10, 2, or)) -> oddsRatioClonalSubclonal
colnames(oddsRatioClonalSubclonal) <- rownames(s)[-1]#paste0("sig", 2:30)
boxplot(oddsRatioClonalSubclonal[,r], log='x', ylim=c(1e-2,1e2), horizontal=TRUE, ylim=c(0,30), at=1:ncol(oddsRatioClonalSubclonal), col=col, lty=1, staplewex=0, pch=16, cex=.5, las=1, xlab="OR subclonal/clonal v sig 1", notch=TRUE)
abline(v=1, lty=3)

#' Early/late v Clonal/subclonal
e <- apply(oddsRatioEarlyLate,2, median, na.rm=TRUE)
f <- apply(oddsRatioClonalSubclonal,2, median, na.rm=TRUE)
par(bty='n')
plot(NA,NA, xlim=c(1,2.3), ylim=range(c(e,f), na.rm=TRUE), xlab="", xaxt='n', ylab="Median OR v sig 1", log='y')
mtext(side=1, at=1:2, c("late/early", "subclonal/clonal"), las=1)
segments(1,e[r],2,f[r], col=col)
w <- which(!is.na(f))
f <- sort(log(f[w]), d=TRUE)
f <- mindist(f, mindist=0.2)
text(2, exp(f), names(f), pos=4)
abline(h=1, lty=3)


#' All comparisons E/L
s <- aperm(PCAWG_sigDecomp_merged, c(3,1,2)) #sigDecomp
s[s==0]<-NA
or <- function(x) x[2,1]/x[1,1]*x[1,2]/x[2,2]
sapply(1:nrow(s), function(j) sapply(1:nrow(s), function(i) apply(s[c(j,i),,1:2]+100, 2, or)), simplify='array') -> orel
morel <- apply(orel, 2:3, median, na.rm=TRUE)
o <- order(morel[1,])
image(morel[o,o], breaks=10^c(-2,seq(-1,1,l=99),2), col=colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(100), x=1:36, y=1:36, xaxt='n', yaxt='n', asp=1, xlab='Signature', ylab='Signature')
axis(side=1, at=seq_along(o), labels=o, las=2)
axis(side=2, at=seq_along(o), labels=o, las=2)

#' All comparisons C/S
s <- abind::abind(asum(s[,,1:3], 3), s[,,4], along=3)
sapply(1:nrow(s), function(j) sapply(1:nrow(s), function(i) apply(s[c(j,i),,1:2]+100, 2, or)), simplify='array') -> orel
morcs <- apply(orel, 2:3, median, na.rm=TRUE)
#o <- order(morcs[1,])
image(morcs[o,o], breaks=10^c(-2,seq(-1,1,l=99),2), col=colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(100), x=1:36, y=1:36, xaxt='n', yaxt='n', asp=1, xlab='Signature', ylab='Signature')
axis(side=1, at=seq_along(o), labels=o, las=2)
axis(side=2, at=seq_along(o), labels=o, las=2)

#' Time v absolute numbers
t1 <- bbTimeTnc[[1]]
w <- t1$MJCN <= 2
v <- allVcf[[1]]
b <- allBB[[1]][t1$cnid[w]]
v <- v[ as.character(as.list(info(v)$CNID)) %in% as.character(t1$cnid[w])]
u <- table(as.character(as.list(info(v)$CNID)), info(v)$MCN)
l <- width(b)
y <- u[as.character(t1$cnid[w]),2]
x <- t1$time.Signature.5[w]
plot(x, y/l, col=as.numeric(seqnames(b)))
segments(x, y/l + sqrt(y)/l, x, y/l - sqrt(y)/l)

timeTotalTnc <- function(vcf, bb=allBB[[meta(header(vcf))['ID',1]]], subset=info(vcf)$CLS != "subclonal", select = c("all","pq")){
	select <- match.arg(select)
	ID <- meta(header(vcf))['ID',1]
	v <- vcf[which(subset)]
	cnid <- as.character(as.list(info(v)$CNID))#sapply(info(v)$CNID, paste, collapse=",")
	if(select=='all'){
		seg <- cnid
		bbSeg <- as.character(seq_along(bb))
	}else{
		o <- findOverlaps(bb, pq, select='first')
		bbSeg <- paste(names(pq)[o], bb$major_cn,bb$minor_cn, round(bb$clonal_frequency*2,1)/2, sep=":")
		seg <- bbSeg[as.numeric(cnid)]
	}
	tnc <- tncToPyrimidine(v)
	n <- table(seg)
	tcn <- info(v)$TCN
	if(is.null(tcn))
		tcn <- (info(v)$MJCN + info(v)$MNCN)
	tab <- table(seg, info(v)$MJCN, info(v)$MCN, tnc)
	
	sig <- apply(tab, 1:3,  `%*%`, t1(tncProbClonal[,1,]))
	sig <- asum(sig[,,,2], 3)
	
	

	cnw <- info(v)$MCN / tcn
	w <- seg %in% bbSeg
	cnw30 <- colSums(cnw[w] * t1(tncProbClonal[,ID,tnc[w]]))
	
	l <- width(bb)
	
	m <- match(dimnames(sig)$seg,bbSeg)

	relW <- l[m]/sum(l[m], na.rm=TRUE)
	
	res <- t1(sig) / relW / rep(cnw30, each=ncol(sig))
	res <- res[order(as.numeric(rownames(res))),]
}

wgdTncClonal <- simplify2array(mclapply(allVcf[isWgd], function(vcf){
					ID <- meta(header(vcf))["ID","Value"]
					if(which.max(mixmdl$posterior[rownames(purityPloidy)==ID,])!=3 | purityPloidy[ID,2] < 2)
						return(matrix(NA,nrow=30, ncol=3))
					c <- classifyMutations(vcf, reclassify='all')
					bb <- allBB[[ID]]
					w <- which(vcf %over% bb[bb$minor_cn==2 & bb$major_cn==2] & c!="subclonal")
					pre <- info(vcf)$MCN[w] == 2 &  info(vcf)$MJCN[w] == 2 &  info(vcf)$MNCN[w] == 2
					post <- info(vcf)$MCN[w] == 1 &  info(vcf)$MJCN[w] == 2 &  info(vcf)$MNCN[w] == 2
					tnc <- tncToPyrimidine(vcf[w])
					
					wgd <- sapply(1:30, function(i){
								c(sum(tncProbClonal[i,ID,tnc[post]]), sum(tncProbClonal[i,ID,tnc[pre]]))
							})
					ci <- sapply(1:1000, function(foo){
								n1 <- rpois(n=length(wgd[1,]),lambda=wgd[1,]+1)
								n2 <- rpois(n=length(wgd[2,]),lambda=wgd[2,]+1)
								(2*n2)/(2*n2 + n1)
							})
					
					ci <- apply(ci,1, function(x) if(all(is.na(x))) rep(NA,2) else quantile(x,c(.025, 0.975), na.rm=TRUE))
					
					cbind(hat=(2*wgd[2,]+0.5)/(1+2*wgd[2,] + wgd[1,]), CI=t1(ci))
					
				}, mc.cores=MCCORES))

relSigActWgd <- t1(sigDecompClonal[,isWgd])/colSums(sigDecompClonal[,isWgd])

#' 2D raw - full TNC
par(mfrow=c(5,6), mar=c(3,3,3,1), mgp=c(2,.5,0))
for(i in (1:30)[-1]){
	x <- wgdTnc[1,1,]
	y <- wgdTnc[i,1,]
	R <- function(theta) matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol=2)
	z <- cbind(x,y) %*% R(pi/4)
	zx <- z[,1]
	zy <- z[,2]
	w <- which((x != .5 | y != .5 ) & !(is.na(x) | is.na(y)) & relSigActWgd[,1] > 0.01 & relSigActWgd[,i] > 0.01) 
	t1 <- seq(0,sqrt(2),l=100)
	try({
				contour(kde2d(x[w], y[w], n=200), col=mg14::colTrans(rep(col,2)[i]), main=paste("Sig",i), ylim=c(0,1), xlim=c(0,1), xlab="Time before diagnosis ( sig 1)", ylab="Time before diagnosis ( sig x)")
				points(x[w], y[w], col=rep(col,2)[i], pch=16, cex=1)
				p <- predict(loess(zy ~ zx,  subset=w , span=1.5, weights=(wgdTnc[1,3,] - wgdTnc[1,2,] + wgdTnc[i,3,]-wgdTnc[i,2,])), newdata=data.frame(zx=t1), se=TRUE)
				ci <- c(p$fit + 1.96*p$se.fit, rev(p$fit - 1.96*p$se.fit))
				polygon( cbind(c(t1, rev(t1))[!is.na(ci)], na.omit(ci)) %*% R(-pi/4), col=paste0(rep(col,2)[i], "88"), border=NA)
				lines(cbind(t1, p$fit)%*% R(-pi/4) , col=rep(col,2)[i], lwd=2)
				abline(0,1)
			})
}

#' 2D raw - clonal TNC
par(mfrow=c(5,6), mar=c(3,3,3,1), mgp=c(2,.5,0))
for(i in (1:30)[-1]){
	x <- wgdTncClonal[1,1,]
	y <- wgdTncClonal[i,1,]
	R <- function(theta) matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol=2)
	z <- cbind(x,y) %*% R(pi/4)
	zx <- z[,1]
	zy <- z[,2]
	w <- which((x != .5 | y != .5 ) & !(is.na(x) | is.na(y)) & relSigActWgd[,1] > 0.01 & relSigActWgd[,i] > 0.01) 
	t1 <- seq(0,sqrt(2),l=100)
	try({
				contour(kde2d(x[w], y[w], n=200), col=mg14::colTrans(rep(col,2)[i]), main=paste("Sig",i), ylim=c(0,1), xlim=c(0,1), xlab="Time before diagnosis ( sig 1)", ylab="Time before diagnosis ( sig x)")
				points(x[w], y[w], col=rep(col,2)[i], pch=16, cex=1)
				p <- predict(loess(zy ~ zx,  subset=w , span=1.5, weights=(wgdTnc[1,3,] - wgdTnc[1,2,] + wgdTnc[i,3,]-wgdTnc[i,2,])), newdata=data.frame(zx=t1), se=TRUE)
				ci <- c(p$fit + 1.96*p$se.fit, rev(p$fit - 1.96*p$se.fit))
				polygon( cbind(c(t1, rev(t1))[!is.na(ci)], na.omit(ci)) %*% R(-pi/4), col=paste0(rep(col,2)[i], "88"), border=NA)
				lines(cbind(t1, p$fit)%*% R(-pi/4) , col=rep(col,2)[i], lwd=2)
				abline(0,1)
			})
}

#' 2D WGD
par(mfrow=c(5,6), mar=c(3,3,3,1), mgp=c(2,.5,0))
for(i in (1:30)[-1]){
	x <- I(-(1-wgdTnc[1,1,]) * age[isWgd])
	y <- I(-(1-wgdTnc[i,1,]) * age[isWgd])
	R <- function(theta) matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol=2)
	z <- cbind(x,y) %*% R(pi/4)
	zx <- z[,1]
	zy <- z[,2]
	w <- (x != .5 | y != .5 ) & !(is.na(x) | is.na(y)) 
	t1 <- seq(-50*sqrt(2),0,l=100)
	try({
				contour(kde2d(x[w], y[w], n=200), col=mg14::colTrans(rep(col,2)[i]), main=paste("Sig",i), ylim=c(-50,0), xlim=c(-50,0), xlab="Time before diagnosis (age x sig 1)", ylab="Time before diagnosis (age x sig x)")
				points(x[w], y[w], col=rep(col,2)[i], pch=16, cex=1)
				p <- predict(loess(zy ~ zx,  subset=w , span=1.5, weights=(wgdTnc[1,3,] - wgdTnc[1,2,] + wgdTnc[i,3,]-wgdTnc[i,2,])), newdata=data.frame(zx=t1), se=TRUE)
				ci <- c(p$fit + 1.96*p$se.fit, rev(p$fit - 1.96*p$se.fit))
				polygon( cbind(c(t1, rev(t1))[!is.na(ci)], na.omit(ci)) %*% R(-pi/4), col=paste0(rep(col,2)[i], "88"), border=NA)
				lines(cbind(t1, p$fit)%*% R(-pi/4) , col=rep(col,2)[i], lwd=2)
				abline(0,1)
			})
}


#' 2D WGD clonal TNC
par(mfrow=c(5,6), mar=c(3,3,3,1), mgp=c(2,.5,0))
for(i in (1:30)[-1]){
	x <- I(-(1-wgdTncClonal[1,1,]) * age[isWgd])
	y <- I(-(1-wgdTncClonal[i,1,]) * age[isWgd])
	R <- function(theta) matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol=2)
	z <- cbind(x,y) %*% R(pi/4)
	zx <- z[,1]
	zy <- z[,2]
	w <- which( (x != .5 | y != .5 ) & !(is.na(x) | is.na(y)) & relSigActWgd[,1] > 0.01 & relSigActWgd[,i] > 0.01) 
	t1 <- seq(-50*sqrt(2),0,l=100)
	try({
				contour(kde2d(x[w], y[w], n=200), col=mg14::colTrans(rep(col,2)[i]), main=paste("Sig",i), ylim=c(-50,0), xlim=c(-50,0), xlab="Time before diagnosis (age x sig 1)", ylab="Time before diagnosis (age x sig x)")
				points(x[w], y[w], col=rep(col,2)[i], pch=16, cex=1)
				p <- predict(loess(zy ~ zx, subset=w , span=0.5, weights=(wgdTncClonal[1,3,] - wgdTncClonal[1,2,] + wgdTncClonal[i,3,]-wgdTncClonal[i,2,])), newdata=data.frame(zx=t1), se=TRUE)
				ci <- c(p$fit + 1.96*p$se.fit, rev(p$fit - 1.96*p$se.fit))
				polygon( cbind(c(t1, rev(t1))[!is.na(ci)], na.omit(ci)) %*% R(-pi/4), col=paste0(rep(col,2)[i], "88"), border=NA)
				lines(cbind(t1, p$fit)%*% R(-pi/4) , col=rep(col,2)[i], lwd=2)
				abline(0,1)
			})
}


#' Sig activity WGD
par(mfrow=c(5,6), mar=c(3,3,3,1), mgp=c(2,.5,0))
for(i in (1:30)[-1]){
	x <- I(-(1-wgdTncClonal[1,1,]) * age[isWgd])
	y <- I(-(1-wgdTncClonal[i,1,]) * age[isWgd])
	R <- function(theta) matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol=2)
	z <- cbind(x,y) %*% R(pi/4)
	zx <- z[,1]
	zy <- z[,2]
	w <- (x != .5 | y != .5 ) & !(is.na(x) | is.na(y)) 
	t1 <- seq(-50*sqrt(2),0,l=100)
	try({
				p <- predict(loess(zy ~ zx, subset=w , span=0.5, weights=(wgdTncClonal[1,3,] - wgdTncClonal[1,2,] + wgdTncClonal[i,3,]-wgdTncClonal[i,2,])), newdata=data.frame(zx=t1), se=TRUE)
				ciup <- p$fit + 1.96*p$se.fit 
				cilo <- p$fit - 1.96*p$se.fit
				xy <- cbind(t1, p$fit)%*% R(-pi/4)
				plot(xy[-1,1]+diff(xy[,1]), diff(xy[,2]) , col=rep(col,2)[i], lwd=2, type='l',main=paste("Sig",i), xlim=c(-50,0), xlab="Time before diagnosis (age x sig 1)", ylab="Time before diagnosis (age x sig x)")
				ciupxy <- cbind(t1, ciup) %*% R(-pi/4)
				ciloxy <- cbind(t1, cilo) %*% R(-pi/4)
				ciupxy2 <- cbind(ciup)
				polygon( c(ciupxy[-1,1]+diff(ciupxy[,1]), rev(ciloxy[-1,1]+diff(ciloxy[,1]))), c(diff(ciupxy[,2]), rev(diff(ciloxy[,2]))), col=paste0(rep(col,2)[i], "88"), border=NA)
				abline(0,1)
			})
}


#' Aggregate all time to WGD
bbTimeTncWgd <- do.call("rbind",lapply(bbTimeTnc, function(x){
			w <- x$MJCN==2 & x$MNCN==2
			sapply(x[w,32:64], weighted.mean, w=x[w,"n"])
		}))

bbTimeTncPq <- simplify2array(mclapply(names(allRanges), function(ID){
			x <- allRanges[[ID]]
			y <- bbTimeTnc[[ID]]
			f <- findOverlaps(pq, x)
			t1(sapply(1:length(pq), function(i){
						w <- subjectHits(f)[which(queryHits(f) == i)]
						colSums(y[w,2:64] * y$n[w]) / sum(y$n[w])
					}))
		}, mc.cores=8))

x <- asum(bbTimeTncPq[,33:63,],3, na.rm=TRUE) / asum(!is.na(bbTimeTncPq[,33:63,]) & !is.nan(bbTimeTncPq[,33:63,]),3)

d <- apply(wgdDeam[1:2,], 2, function(x) {w <- x/sum(x); a <- 1/sum(w/sum(w) * c(1:2)/4); (4-a)/(4-2)})

f <- t1(asum(sigDecomp[,,1:3],3))

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
t1 <- table(unlist(u$genes), u$type, as.character(seqnames(u)))
barplot(sort(rowSums(t1), d=TRUE), las=2)
names(age) <- sampleIds
par(mfrow=c(2,3))

pdf("times.pdf",3,3, pointsize=8); par(bty="n", mgp=c(2,0.5,0), mar=c(3,3,1,1))
#### Age of WGD!!!!
y <- (wgd[1,]/wgdWeight)/(nDeamTrunk*avgWeightTrunk)*age
plot(age, y, ylim=c(0,100), xlab="Age at diagnosis", ylab="Age of WGD"); abline(0,1)
e <- (t1(sapply(wgd[1,], function(w) qpois(c(0.025,0.975), w)))/wgdWeight) /(nDeamTrunk*avgWeightTrunk)*age
segments(age, e[,1], age, e[,2], col='grey')
d <- density(na.omit(age)); lines(d$x, d$y*.1/max(d$y) * 100)
d <- density(na.omit(y)); lines(d$y*.1/max(d$y) * 100, d$x)
legend("topleft", legend=sum(!is.na(wgd[1,])), lty=1, bty="n", col=1)

clr <- c(col[1:3], "#888888")
#### Drivers
for(gene in c("TP53","PIK3CA","KRAS","NOTCH1","CDKN2A")){
	t1 <- u[sapply(u$genes, function(x) gene %in% x)]
	f <- age[t1$sample]/lTotal[t1$sample]
	y <- (end(t1)+start(t1)) /2*f
	
	plot(age[t1$sample], y , xlim=c(0,100), ylim=c(0,100), col=clr[as.numeric(seqnames(t1))], xlab="Age at diagnosis", ylab="Age at mutation", main=gene, pch=".")
	abline(0,1)
	e <- cbind(start(t1),end(t1)) *f
	segments(age[t1$sample], e[,1], age[t1$sample], e[,2], col=clr[as.numeric(seqnames(t1))])
	f <- function(x) paste(names(x), x)
	legend("topleft", legend=f(table(seqnames(t1))), lty=1, col=clr[1:4], bty="n")
	d <- density(na.omit(age[t1$sample])); lines(d$x, d$y*.1/max(d$y) * 100)
	w <- !is.na(e[,1])
	d <- sapply(coverage(GRanges(seqnames(t1)[w], IRanges(e[w,1],e[w,2])), weight=1/(e[w,2]-e[w,1]+1)), function(x) {y <- numeric(100); n <- as.numeric(x); y[1:min(length(n),100)] <- n[1:min(length(n),100)];y})
	x <- apply(d, 1,cumsum)/length(t1)*500
	for(i in nrow(x):1) polygon(c(0,x[i,]), c(1,1:100), col=mg14::colTrans(clr[i]), border=NA, )
}
dev.off()


tSubclones <- sapply(sampleIds, function(ID){
			x <- lBranch[[ID]]
			clone <- which.max(allClusters[[ID]]$n_ssms)
			if(clone==1) return(c(NA,NA))
			s <- x[1:(clone-1)]
			c <- x[clone:length(x)]
			c(max=sum(s), mean=mean(s))/sum(x)
		})




f <- function(x, t=.5, f=2) if(x>t) f else 1

F <- function(x, t=.5, f=2) ifelse(x<t, x, t+f*(x-t))

time <- function(x, t=.5, f=2) F(x,t,f) / (F(1,t,f))

f <- t1(sapply(1:1000, function(foo){
	t0 <- 0.9
	a <- 0.1
	t <- runif(1)
	mu=100
	x12 <- rpois(1,mu*t) 
	x11 <- rpois(1,mu*(1-t)) + rpois(1, mu) 
	t1 <- 3*x12/(x11+2*x12)
	x22 <- rpois(1,mu*F(t, t0, a))
	x21 <- rpois(1,mu*(F(1,t0,a)-F(t,t0,a))) + rpois(1,mu*F(1,t0,a)) 
	t2 <- 3*x22/(x21 + 2*x22)
	return(c(t1,t2))
}))

plot(f[,1], f[,2], log='')

c <- cobs(f[,1],f[,2], constraint="increase", pointwise = matrix(c(0,0,0,1,0,1), ncol=3),  degree=1, nknots=10)
p <- predict(c)
lines(p, col='red')
diff(c$coef)/diff(c$knots)


#' Poisson regression v sigActivity
ID <- sampleIds[1]
f <- glm(sigTableClonal[,ID] ~ S[,sigActivityHere[,tumourType[ID]]==1]-1, family=poisson)

dir.create("../scratch/sig")
for(f in colnames(sigTable)) write.table(sigTable[,f,], file=paste0("../scratch/sig/", f,".txt"), quote=FALSE, sep="\t")


w <- runif(10, 0, 1000)
Y <- matrix(rpois(96*10, S[,1:10] %*% w), ncol=10)
n <- nn.poisglm(Y[,1,drop=FALSE], S[,1:10])

g <- glm(Y[,1] ~ S[,1:10]-1, family=poisson(link=identity))


ab <- sigDecomp[,,1:2]
cd <- asum(sigDecomp[,,1:2],1)
sigOddsEarlyLate <- ab[,,1]/ab[,,2] * cd[,2]/cd[,1]

ab[ab<=10] <- NA
pseudo <- 1
sigOddsEarlyLateSig1 <- t1(t1((ab[,,1]+pseudo)/(ab[,,2]+pseudo)) * (ab[1,,2]+pseudo)/(ab[1,,1]+pseudo))
rownames(sigOddsEarlyLateSig1) <- rownames(sigDecomp)

mg14:::ggPlot(t1(sigOddsEarlyLateSig1), rep(rownames(sigOddsEarlyLateSig1), each=ncol(sigOddsEarlyLateSig1)), log='y', las=2)


plot(sigDecomp[1,7,], sigDecomp[5,7,])
segments(sigDecomp[1,7,], sigDecomp[5,7,] - sqrt(sigDecompVar[5,7,]), sigDecomp[1,7,] , sigDecomp[5,7,] + sqrt(sigDecompVar[5,7,]))
segments(sigDecomp[1,7,] - sqrt(sigDecompVar[1,7,]), sigDecomp[5,7,], sigDecomp[1,7,] + sqrt(sigDecompVar[1,7,]), sigDecomp[5,7,] )



breakClass <- function(bb){
	before <- c(0,diff(bb$copy_number))
	after <- c(diff(bb$copy_number),0)
	table(early=factor(before %% 2 == 0 & before != 0, levels=c(TRUE, FALSE)), chr=factor(as.character(seqnames(bb)), levels=c(1:22, "X","Y")))
}

breakClasses <- sapply(allBB, breakClass, simplify='array')

#' cancerTiming
library(cancerTiming)
timeSpell <- function(vcf, bb=allBB[[meta(header(vcf))['ID',1]]], subset=info(vcf)$CLS != "subclonal", select = c("all","pq")){
	select <- match.arg(select)
	ID <- meta(header(vcf))['ID',1]
	v <- vcf[which(subset)]
	cnid <- as.character(as.list(info(v)$CNID))#sapply(info(v)$CNID, paste, collapse=",")
	if(select=='all'){
		seg <- cnid
		bbSeg <- as.character(seq_along(bb))
	}else{
		o <- findOverlaps(bb, pq, select='first')
		bbSeg <- paste(names(pq)[o], bb$major_cn,bb$minor_cn, round(bb$clonal_frequency*2,1)/2, sep=":")
		seg <- bbSeg[as.numeric(cnid)]
	}
	x <- getAltCount(v)
	n <- rowSums(getTumourCounts(v))
	no
	sapply(unique(seq), function(s){
				ix <- seq==s
				xx <- x[ix]
				nn <- n[ix]
			})
	
}

#' Time and fit signatures simulatenously

sigTime <-  function(vcf, bb=allBB[[meta(header(vcf))['ID',1]]], subset=info(vcf)$CLS != "subclonal", select = c("all","pq")){
	select <- match.arg(select)
	ID <- meta(header(vcf))['ID',1]
	v <- vcf[which(subset)]
	cnid <- as.character(as.list(info(v)$CNID))#sapply(info(v)$CNID, paste, collapse=",")
	if(select=='all'){
		seg <- cnid
		bbSeg <- as.character(seq_along(bb))
	}else{
		o <- findOverlaps(bb, pq, select='first')
		bbSeg <- paste(names(pq)[o], bb$major_cn,bb$minor_cn, round(bb$clonal_frequency*2,1)/2, sep=":")
		seg <- bbSeg[as.numeric(cnid)]
	}
	tnc <- tncToPyrimidine(v)
	n <- table(seg)
	tcn <- info(v)$TCN
	if(is.null(tcn))
		tcn <- (info(v)$MJCN + info(v)$MNCN)
	tab <- table(seg , MCN=factor(pmin(info(v)$MCN,3), levels=1:2), tnc)	
	
	s <- intersect(bbSeg[bb$copy_number <=4 & bb$major_cn <=2], rownames(tab))
	
	O <- (cbind(bb$copy_number, (bb$minor_cn>0) + (bb$major_cn>0) ) * width(bb)/1e6)[as.numeric(s),]

	T <- matrix(tab[s,,] / rep(O, 96), ncol=96)
	
	
	#D <- outer(O, S[,sigActivityHere[,tumourType[ID]]==1])
	#D <- matrix(D, ncol=3)
	D <- S[,sigActivityHere[,tumourType[ID]]==1]
	
	E <- nmSolve(t1(T), D, maxIter=10000, tol=1e-4)
	
	E <- array(E, dim=c(dim(E)[1], dim(tab[s,,])[-3]))
}

#' MMR
mmrGenes <- c("MSH2","MSH3","MSH6","MLH1","PMS2", "EXO1","PCNA","RPA1","RPA2","RPA3","LIG1","POLD1","POLD2","POLD3","POLD4")
mmrSignatures <- grep(paste0("\\.",paste(c(6, 15, 20, 26), collapse="|")), dimnames(PCAWG_sigDecomp)[[3]], value=TRUE)

mmrGenotypes <- allGenotypes[intersect(mmrGenes, rownames(allGenotypes)),,,,] 
g <- which(asum(mmrGenotypes, c(1:4))!=0)
round(asum(PCAWG_sigDecomp[g,,mmrSignatures], 3))/round(PCAWG_sigDecomp[g,,1]+1)

s <- which(asum(mmrGenotypes[,,,,g],c(1:2,4))["subclonal",]!=0)

round(asum(PCAWG_sigDecomp[g[s],,mmrSignatures], 3))/round(PCAWG_sigDecomp[g[s],,1]+1)

w <- which(asum(PCAWG_sigDecomp[,,mmrSignatures],2:3)>1)

table(MMR_mut=seq_along(sampleIds)%in%g, MMR_sig=seq_along(sampleIds)%in%w)

#' POLE
e <- which(asum(allGenotypes["POLE",,,,], 1:3) > 0)
PCAWG_sigDecomp_merged[e,,"S10"]
which(asum(PCAWG_sigDecomp_merged[,,"S10"],2)>1)


w <- which(asum(PCAWG_sigDecomp_merged[,,"S10"],2)>1)
PCAWG_sigDecomp_merged[w,,"S10"]
round(PCAWG_sigDecomp_merged[w,,"S10"])/round(PCAWG_sigDecomp_merged[w,,1]+1)

table(POLE_mut=seq_along(sampleIds)%in%e, POLE_sig=seq_along(sampleIds)%in%w)


#' BRCA
b <- which(asum(allGenotypes[c("BRCA1","BRCA2"),,,"TRUE",], 1:3) > 0)
PCAWG_sigDecomp_merged[b,,"S3"]

round(PCAWG_sigDecomp_merged[b,,"S3"])/round(PCAWG_sigDecomp_merged[b,,1]+1)

w <- which(asum(PCAWG_sigDecomp_merged[,,"S3"],2)>1)
table(BRCA_mut=seq_along(sampleIds)%in%b, BRCA_sig=seq_along(sampleIds)%in%w)


t1(asum(allGenotypes[c('BRCA1','BRCA2'),,,,b],c(1:3)))

g <- asum(allGenotypes,2:4)
cbind(asum(g[,seq_along(sampleIds)%in%w],2),asum(g[,!seq_along(sampleIds)%in%w],2))

table(newTumourTypes,TP53=g['TP53',]>0, BRCA_sig=seq_along(sampleIds)%in%w)


t1 <- table(newTumourTypes,TP53=g['TP53',]>0, BRCA_sig=seq_along(sampleIds)%in%w)

s <- asum(PCAWG_sigDecomp_merged[,1:3,],2)
f <- lm(log(age) ~ log(s+1) + newTumourTypes + log(purityPloidy[sampleIds,2]))
summary(f)

g <- lm(age ~ s + newTumourTypes + purityPloidy[sampleIds,2])
summary(g)

plot(f)

summary(lm(s[,5] ~ age + newTumourTypes))

summary(lm(rowSums(s) ~ age + newTumourTypes))


#' WGD - new tumour Types
col <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"))
col <- rep(col, each=2)[1:nlevels(newTumourTypes)]
names(col) <- levels(newTumourTypes)

mg14:::ggPlot(tWgd, newTumourTypes, las=2)#, col=unlist(sapply(levels(newTumourTypes)[order(aggregate(tWgd, list(newTumourTypes), median, na.rm=TRUE)[,2], na.last=TRUE)], function(t) rep(col[t], table(newTumourTypes)[t]))), pch=19, ylab="Relative time", las=2)
mg14:::ggPlot(age*(1-tWgd), newTumourTypes, col=unlist(sapply(levels(newTumourTypes)[order(aggregate(age*(1-tWgd), list(newTumourTypes), median, na.rm=TRUE)[,2], na.last=TRUE)], function(t) rep(col[t], table(newTumourTypes)[t]))), pch=19, ylab="Time lag", las=2)

par(mfrow=c(6,7), cex=0.5, mar=c(3,3,3,1), mgp=c(2,0.5,0))
for(t1 in levels(newTumourTypes)) {plot(age[newTumourTypes==t1], s[newTumourTypes==t1,1], ylim=c(0,5000), xlim=c(0,100), ylab='#mut', xlab='age', main=t1, pch=19); points(age[newTumourTypes==t1], s[newTumourTypes==t1,5], pch=19, col='red')}

chr <- unlist(sapply(names(bbTimeTnc), function(n) {
			s <- as.character(seqnames(allBB[[n]]))
			as.character(s[as.numeric(as.character(bbTimeTnc[[n]]$cnid))])
		}
))

chr <- factor(chr, levels=c(1:22, "X","Y"))

par(mfrow=c(6,7), cex=0.5, mar=c(3,3,3,1), mgp=c(2,0.5,0))
for(t1 in levels(newTumourTypes)) {
	w <- newTumourTypes[bbTimeTncAll$sampleId]==t1 &! isWgd[bbTimeTncAll$sampleId] & bbTimeTncAll$n > 100 
	v <- newTumourTypes[bbTimeTncAll$sampleId]==t1 & isWgd[bbTimeTncAll$sampleId] & bbTimeTncAll$n > 100 
	try({
				boxplot(bbTimeTncAll$time.Signature.PCAWG.1[w] ~ chr[w], main=t1, at=1:24, xlim=c(0,26), ylim=c(0,1))
				boxplot(bbTimeTncAll$time.Signature.PCAWG.1[v],  at=25, border='red', add=TRUE, width=1)
			})
}


apply(sapply(1:100, function(foo){
t <- rmultinom(1, 40000, S[,6]) + rmultinom(1, 10000, S[,20]) + rmultinom(1, 10000, S[,26]) +  rmultinom(1, 2000, S[,1])  +  rmultinom(1, 1000, S[,17]) +  rmultinom(1, 2000, S[,5])
nmSolve(as.matrix(t, ncol=1), S[,c(1,5,6,17,20,26)])
}),1,quantile)



w <- which(newTumourTypes=="Breast-AdenoCa")
f <- ssnmf(sigTableClonal[,w]+.001, S[,c(1,5)], s=5, maxIter=1000)

L <- t1(sapply(p, function(pp) dbinom(altCount[w],tumDepth[w],pp)))

pi <- rep(1, nrow(L))
for(i in 1:100){
	post <- L * piState[cnStates[1:k,"state"]] * cnStates[1:k,"piCn"] * pi
	post <- t1(post)/colSums(post)
	pi <- colMeans(post) / piState[cnStates[1:k,"state"]] / cnStates[1:k,"piCn"]
	print(pi)
}



hist(altCount[w]/tumDepth[w], n=50, col='black', probability=TRUE)
f <- rowSums(sapply(1:6, function(i) dbinom(0:50, 50, p[i])*piState[cnStates[i,"state"]] * cnStates[i,"piCn"] * pi[i]))
lines(seq(0,1,length=length(f)), f*50, col='red')


tAdjust <- function(t, t0, r){
	x <- ifelse(t < t0,	t, (t0 + (t-t0)*r))
	x/(t0+(1-t0)*r)
}

tAdjustInv <- function(t, t0, r){
	t1 <- t0 / (t0 + (1-t0) * r)
	x <- ifelse(t < t1,	(t0 + (1-t0) * r) * t, ((t0 + (1-t0) * r) * t - t0 + t0*r)/r )
}

that <- seq(0,1,l=100)
tadj <- sapply(1:1000, function(foo){
			r <- 5 #rgamma(1, 5,1)
			t0 <- rbeta(1,5,1)
			tAdjustInv(that, t0, r)
		})

qadj <- apply(tadj, 1, quantile, c(0.025,.5,.975))
plot(that, qadj[2,], type='l', lwd=2)
lines(that, qadj[1,], type='l')
lines(that, qadj[3,], type='l')



wgdFit <- function(ID){
	v <- allVcf[[ID]]
	w <- which(info(v)$MJCN==2 & info(v)$MNCN==2& sapply(info(v)$CNID, length)==1 & isDeamination(v))
	majcni <- 2
	mincni <- 2
	purity <- purityPloidy[ID,1]
	cfi <- purityPloidy[ID,1]
	effCnTumor <- sum((majcni + mincni)*cfi)
	cnNormal <- 2
	effCnNormal <- as.numeric(cnNormal) * (1-purity)
	
	cnStates <- matrix(0, nrow=10000, ncol=5)
	colnames(cnStates) <- c("state","m","f","n.m.s","pi.m.s")
	
	if(any(is.na(majcni))) next
	
	multFlag <- rep(FALSE, length(cfi))
	
	if(length(cfi)>1){ # multiple (subclonal) CN states, if so add clonal option (ie. mixture of both states), subclonal states only change by 1..delta(CN)
		majcni <- c(rep(1, 2) + c(0, diff(majcni)), max(majcni))
		mincni <- c(rep(1, 2) + c(0, diff(mincni)), max(mincni))
		cfi <- c(cfi, purity)
		multFlag <- c(multFlag, TRUE)
	}
	
	clusters <- allClusters[[ID]]
	
	a <- sapply(clusters$proportion, function(p) all(abs(p-cfi) > 0.05)) # subclone(s) not coinciding with CN change
	if(any(a)){ # assume subclones have derived from most abundant CN state
		majcni <- c(majcni, rep(majcni[which.max(cfi)]>0, sum(a))+0)
		mincni <- c(mincni, rep(mincni[which.max(cfi)]>0, sum(a))+0)
		cfi <- c(cfi, clusters$proportion[a])
		multFlag <- c(multFlag, rep(FALSE, sum(a)))
	}
	totcni <- majcni+mincni
	pi.s <- sapply(cfi, function(p) ifelse(min(abs(clusters$proportion - p)) < 0.05, clusters$n_ssms[which.min(abs(clusters$proportion - p))], 1))
	pi.s <- pi.s/sum(pi.s)
	if(all(totcni==0)) next
	
	k <- 0
	for( j in seq_along(majcni)){
		if(majcni[j]==0) {
			f <- m <- pi.m.s <- n.m.s <- 0 # allele frequency
			l <- 1
		}else{
			l <- 1:max(majcni[j], mincni[j]) # mincni>majcni can occur if minor allele changes subclonally
			m <- l
			n.m.s <- rep(1, length(l)) #multiplicity of cn states
			if(!multFlag[j]){ # single subclone, or no subclonal cn
				f <- l * cfi[j] / (effCnTumor + effCnNormal)
				if(mincni[j] > 0)
					n.m.s[1:min(majcni[j], mincni[j])] <- 2
				pi.m.s <- n.m.s/sum(n.m.s)
			}else{ # coexisting cn subclones, use mixture
				f <- (cfi[1] * pmin(l, majcni[1]) + cfi[2] * pmin(l,majcni[2])) / (effCnTumor + effCnNormal) # Variant on major allele
				if(mincni[j] > 0){ # Variant on minor allele
					ll <- 1:mincni[j]
					f <- c(f,  (cfi[1] * pmin(ll, mincni[1]) + cfi[2] * pmin(ll,mincni[2])) / (effCnTumor + effCnNormal))
					d <- !duplicated(f) # remove duplicates
					n.m.s <- table(f)[as.character(f[d])] # table duplicates
					f <- f[d]
					#piCnState[1:mincni[j]] <- 2 
					m <- c(l,ll)[d]
					l <- seq_along(m)
				}
				pi.m.s <- n.m.s/sum(n.m.s)
			}
		}
		cnStates[k + l,"state"]=j
		cnStates[k + l,"m"]=m
		cnStates[k + l,"f"]=f
		cnStates[k + l,"pi.m.s"]=pi.m.s
		cnStates[k + l,"n.m.s"]=n.m.s
		k <- k + length(l)
	}
	hh <-w
	altCount <- getAltCount(v)
	tumDepth <- getTumorDepth(v)
	
	L <- matrix(sapply(pmin(cnStates[1:k,"f"],1), function(pp) dbinom(altCount[hh],tumDepth[hh],pp) + .Machine$double.eps), ncol=k)
	P.m.sX <- cnStates[1:k,"pi.m.s"]
	for(em.it in 1:100){
		P.xsm <- L * rep(pi.s[cnStates[1:k,"state"]] * P.m.sX, each=nrow(L)) # P(X,s,m)
		P.sm.x <- P.xsm/rowSums(P.xsm) # P(s,m|Xi)
		P.sm.X <- colMeans(P.sm.x) # P(s,m|X) / piState[cnStates[1:k,"state"]] / cnStates[1:k,"pi.m.s"]
		P.s.X <- sapply(split(P.sm.X, cnStates[1:k,"state"]), sum)
		P.m.sX <- P.sm.X / P.s.X[cnStates[1:k,"state"]]
	}
	
#			boot <- sapply(1:100, function(foo) {Lb <- L[sample(1:nrow(L), replace=TRUE),]
#						P.m.sX <- cnStates[1:k,"pi.m.s"]
#						for(em.it in 1:100){
#							P.xsm <- Lb * rep(pi.s[cnStates[1:k,"state"]] * P.m.sX, each=nrow(L)) # P(X,s,m)
#							P.sm.x <- P.xsm/rowSums(P.xsm) # P(s,m|Xi)
#							P.sm.X <- colMeans(P.sm.x) # P(s,m|X) / piState[cnStates[1:k,"state"]] / cnStates[1:k,"pi.m.s"]
#							P.s.X <- sapply(split(P.sm.X, cnStates[1:k,"state"]), sum)
#							P.m.sX <- P.sm.X / P.s.X[cnStates[1:k,"state"]]
#						}
#						return(P.m.sX)})
#			
#			
#			S <- L * rep(pi.s[cnStates[1:k,"state"]], each=nrow(L)) / rowSums(L * rep(pi.s[cnStates[1:k,"state"]] * P.m.sX, each=nrow(L))) # Sufficient statistic
#			S <- sapply(1:j, function(jj) rowSums(S[,cnStates[1:k,"state"]==jj, drop=FALSE]) - S[,cnStates[1:k,"state"]==jj, drop=FALSE][,1])
#			V <- solve(t(S)%*% S) # Cramer-Rao bound on the variance of pi.m.s
#			a <- P.m.sX*(1-P.m.sX)/diag(V) * P.m.sX - P.m.sX
#			b <- a * (1-P.m.sX)/P.m.sX
	
	
	
	P.sm.x[apply(is.na(P.sm.x)|is.nan(P.sm.x),1,any),] <- NA
	
	pwr <- sapply(1:k, function(j) power(cnStates[j,"f"], round(mean(tumDepth[w]))))
	
	opp <- 3-cnStates[1:k,"m"]
	
	T <- P.sm.X/opp/pwr
	
	T <- T/sum(T)
	
	names(T) <- make.names( c("LateClonal","EarlyClonal", rep("Subclonal", length(T)-2)), unique=TRUE)
	
	# correct for acceleration
	max.a <- 10
	

	
	ctn <- sapply(1:100, function(foo){
				a <- runif(1,1,max.a)
				ta <- runif(1,.5,1)
				ctn <- c(correct.t(P.sm.X, ta, a), rep(a, length(T)-2))
			})
	T.ct <- T/ctn
	T.ct <- t(T.ct)/ colSums(T.ct)
	colnames(T.ct) <- names(T)
	return(list(T=T, T.ct=T.ct))
}

correct.t <- function(pi, ta, a){
	pi <- pi[1:2]/sum(pi[1:2])
	t <- 2*pi[2]/(pi[1]+2*pi[2])
	if(t < ta)
		crct <- c(1+(a-1)*(ta-t)/(1-t),1)
	else
		crct <- c(a,1+(a-1)*(t-ta)/t)
	return(crct)
}

t.accel <- function(pi, ta, a){
	tmin <- (pi[1]*(a-1)*ta + 2*pi[2]*a)/(pi[1]*a + 2*pi[2]*a)
	tmax <- (2*pi[2] * (ta + a*(1-ta)))/(pi[1]+ 2*pi[2])
	if(ta < tmin) tmin else tmax
}

wgd.fits <- sapply(sampleIds[isWgd], function(ID) {try(wgdFit(ID))})

tWgdCr <- sapply(wgd.fits, function(x){
			if(class(x)=="try-error") return(matrix(NA, ncol=6, nrow=3))
			else{
				T.ct <- x$T.ct
				T <- matrix(0, ncol=6, nrow=100)
				T[,1:ncol(T.ct)] <- T.ct
				apply(T, 2, quantile, c(0.025,.5,0.975))
			}
		}, simplify='array') 


x <-colSums(tWgdCr[2,c(1,3:6),]) * age[isWgd] + .5
o <- sapply(split(x, newTumourTypes[isWgd]), median, na.rm=TRUE)
l <- sub(".NA","",names(o)[order(o)])
t1 <- table(factor(newTumourTypes[isWgd], levels=l))
w <- t1 >= 5 & names(t1)!="" & ! is.na(o[order(o)])
y <- factor(newTumourTypes[isWgd], levels=l[w])
par(mar=c(8,5,3,1), mfrow=c(2,1), bty='L')
boxplot(x ~ y, las=2, log='', ylab="Approx. years before diagnosis", main="WGD", staplewex=0, lty=1, pch=NA, ylim=c(0,40))
boxplot(colSums(tWgdCr[2,3:6,]) * age[isWgd] ~ y, las=2, log='',ylab="Approx. years before diagnosis", main="Subclonal diversification", staplewex=0, lty=1, pch=NA, ylim=c(0,10))


pdf("WGD.pdf", 5,3, pointsize=8)
par(mar=c(9,3,3,1), mgp=c(2,.5,0), bty='L')
boxplot(x ~ y, las=2, log='', ylab="Approx. years before diagnosis", main="WGD", staplewex=0, lty=1, pch=NA, ylim=c(0,40))
dev.off()


wgdRealTime <- sapply(1:3, function(i){
			x <-colSums(tWgdCr[i,c(1,3:6),]) * age[isWgd] + .5
			o <- sapply(split(x, newTumourTypes[isWgd]), median, na.rm=TRUE)
			l <- sub(".NA","",names(o)[order(o)])
			t1 <- table(factor(newTumourTypes[isWgd], levels=l))
			w <- t1 >= 5 & names(t1)!="" & ! is.na(o[order(o)])
			y <- factor(newTumourTypes[isWgd], levels=l[w])
			sapply(split(x, y), median, na.rm=TRUE)
		})

subclRealTime <- sapply(1:3, function(i){
			x <-colSums(tWgdCr[i,c(3:6),]) * age[isWgd] + .5
			o <- sapply(split(x, newTumourTypes[isWgd]), median, na.rm=TRUE)
			l <- sub(".NA","",names(o)[order(o)])
			t1 <- table(factor(newTumourTypes[isWgd], levels=l))
			w <- t1 >= 5 & names(t1)!="" & ! is.na(o[order(o)])
			y <- factor(newTumourTypes[isWgd], levels=l[w])
			sapply(split(x, y), median, na.rm=TRUE)
		})

bottleIncidence <- function(n=100000, n1=20, n2=100){
	t1 <- rexp(n, 1e-2/n1)
	t2 <- rexp(n, 1e-2/n2)
	return(t1+t2)
}


#Overview pie charts
pdf("timing-pie.pdf", 10,20, pointsize=8)
t1 <- cut(bbTimeTncAll$time.Signature.1, seq(0,1,0.1))
par(mfrow=c(nlevels(newTumourTypes)+1,nlevels(chr)+1), cex=0.5, mar=c(0,0,0,0), mgp=c(2,0.5,0), bty="n")
for(l in levels(newTumourTypes)) {
	plot(NA,NA, xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1)); text(0.5,0.5, l)
	w <- newTumourTypes[bbTimeTncAll$sampleId]==l #&! isWgd[bbTimeTncAll$sampleId] & bbTimeTncAll$n > 10
	for(c in levels(chr)){
		v <- chr==c
		#v <- newTumourTypes[bbTimeTncAll$sampleId]==l & isWgd[bbTimeTncAll$sampleId] & bbTimeTncAll$n > 100 
		try({
					pie(table(t1[which(w & v)])+.1, radius = 2*sqrt(sum(v&w, na.rm=TRUE)/sum(newTumourTypes[bbTimeTncAll$sampleId]==l , na.rm=TRUE)), edges=50, labels="", col=rainbow(10), border=NA)
					#boxplot(bbTimeTncAll$time.Signature.PCAWG.1[v],  at=25, border='red', add=TRUE, width=1)
				})
	}
}
dev.off()

Aclonal <- PCAWG_sigDecomp_merged[,1:2,]
Asubclonal <- PCAWG_sigDecomp_merged[,c(1,4),] + PCAWG_sigDecomp_merged[,c(2,4),] + PCAWG_sigDecomp_merged[,c(3,4),] 
Asubclonal[,2,] <- Asubclonal[,2,]/3

A <- abind::abind(Aclonal, Asubclonal, along=4)
zero.na <- function(x){y <- x; y[y==0] <- NA; return(y)}
B <- zero.na(A)

levels(newTumourTypes)[levels(newTumourTypes) == ""] <- "Other"
newTumourTypes <- droplevels(newTumourTypes)

layout(matrix(seq(1, (2+nlevels(newTumourTypes)) * (1+ dim(A)[3]) ), nrow = 2+nlevels(newTumourTypes), byrow=TRUE), width=c(5, rep(1, dim(A)[3])), height=c(2, rep(1, nlevels(newTumourTypes))))
par( bty="n", mar=rep(0,4)+.1)
plot.new()
for(j in 1:dim(A)[3]){
	plot(NA,NA, xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
	text(.5,.5,sub("Signature.(PCAWG.)*","",dimnames(A)[[3]][j]), srt=45)
}
for(i in c(levels(newTumourTypes), "Pan-Cancer")){
	plot(NA,NA, xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0,2), ylim=c(0,1))
	text(2,.5,i, pos=2)
	for(j in 1:dim(A)[3]){
		w <- which(newTumourTypes == i)
		if(i == "Pan-Cancer") w <- TRUE
		pseudo.count <- 1
		a <- (B[,1,j,][w,,drop=FALSE] + pseudo.count) / asum(B[,1,,][w,,,drop=FALSE] + pseudo.count, 2, na.rm=TRUE)
		b <- (B[,2,j,][w,,drop=FALSE] + pseudo.count) / asum(B[,2,,][w,,,drop=FALSE] + pseudo.count, 2, na.rm=TRUE)
		if(all(is.na(B[w, 1, j,] )) | all(na.omit(B[w, 1, j,] )==0 | sum(!is.na(B[w,1,j,1])) <= 2)){
			plot.new()
			next
		}
		plot(NA,NA, xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0,2), ylim=c(0.2,5), log='y')
		par(xpd=TRUE)
		abline(h=c(0.5,1,2), col='grey')
		par(xpd=NA)
		z <- boxplot(b/a, plot=FALSE)
		z$stats <- apply(b/a, 2, quantile, c(0.05,0.25,0.5,0.75,0.95), na.rm=TRUE)
		t1 <- 1.05
		t2 <- 1.5
		c <- c("#990033","#FF6666","#CCCCCC","#3399FF","#006699")#RColorBrewer::brewer.pal(5,"RdBu")
		c[3] <- "grey"
		f <- sapply(1:2, function(i) {if(all(is.na(z$conf[,i]))) NA else if(all(z$conf[,i] > t2)) 1 else if(all(z$conf[,i] > t1)) 2 else if(all(z$conf[,i] < 1/t2)) 5 else if(all(z$conf[,i] < 1/t1)) 4 else 3})
		f <- c[f]
		c <- ifelse(f == "white","grey","black")
		bxp(z, pch=NA,staplewex=NA, lty=1, xaxt="n", yaxt="n", add=TRUE, xlab="", ylab="", boxfill=f, at=0.5 + 0:1, boxcol=f, medlwd=1,medcol="white", whiskcol=f, whisklwd=2, boxwex=0.2 + 0.75 * colSums(A[w, 1, j,,drop=FALSE]!=0, na.rm=TRUE)/colSums(!is.na(A[w, 1, j,,drop=FALSE])), boxlty=1)
	}
}

dev.copy(device=pdf, file="Signature-overview-2.pdf", width=8, height=7, pointsize=8)
dev.off()

# median v width
x <- sapply(levels(newTumourTypes), function(i){ 
			sapply(1:dim(A)[3], function(j){
			w <- which(newTumourTypes == i)
			if(length(w)==0) return(rep(NA,3))
			a <- (A[,1,j,][w,,drop=FALSE] + 10) / asum(A[,1,,][w,,,drop=FALSE] + 10, 2)
			b <- (A[,2,j,][w,,drop=FALSE] + 10) / asum(A[,2,,][w,,,drop=FALSE] + 10, 2)
			r <- b/a
			r[r==1] <- NA
			med <- apply(r, 2, median, na.rm=TRUE)
			pre <- colSums(A[w, 1, j,1,drop=FALSE]!=0, na.rm=TRUE)/colSums(!is.na(A[w, 1, j,1,drop=FALSE]))
			return(c(med, pre))
		}, simplify='array')}, simplify='array')


sapply(levels(newTumourTypes), function(i){ 
			w <- which(newTumourTypes == i)
			if(length(w) == 0) return(rep(NA, dim(B)[3]))
			colSums(!is.na(B[,2,,2][w,,drop=FALSE]))
		})


pi0 <- list()
tumourTypes <- character(0)
p <- "../data/20161029_Timing_BB_DPclust/"
for(f in dir(p, pattern="*.txt")){
	t <- strsplit(f, "_")[[1]][2]
	d <- read.table(paste0(p,f), header=TRUE, quote="\"", sep="\t") 
	pi0[[t]] <- cbind(d, tumourType=rep(t, nrow(d)))
	tumourTypes <- c(tumourTypes, t)
}

tumourTypes <- factor(tumourTypes)

pi0_agg <- Reduce("rbind",pi0)
pi0_agg$chr <- factor(pi0_agg$chr, levels=1:22)

#deviation coding
globalLogitModel <- lm(car::logit(pi0) ~ chr:tumourType + tumourType -1, contrasts=list(chr=contr.treatment(22,contrasts = FALSE, base=22)), weights = 1/(car::logit(pi0_agg$uCI) - car::logit(pi0_agg$lCI) + 0.5), data=pi0_agg, subset=pi0_agg$type=="SingleGain" & ! pi0_agg$Sample %in% names(which(isWgd)))

library(lme4)
lmmFit <- lmer(car::logit(pi0) ~ tumourType - 1 + (tumourType -1| chr), data=pi0_agg, subset=pi0_agg$type=="SingleGain" & pi0_agg$tumourType %in% levels(pi0_agg$tumourType)[1:8])


s <- summary(globalLogitModel)
p <- s$coefficients[-(1:nlevels(pi0_agg$tumourType)), 4]
q <- p.adjust(p, "BH")

layout(matrix(seq(1, (2+nlevels(tumourTypes)) * (1+ 22) ), nrow = 2+nlevels(tumourTypes), byrow=TRUE), width=c(8, rep(1, 22)), height=c(2, rep(1, nlevels(tumourTypes))))
par( bty="n", mar=rep(0,4)+.1, lend=3, ljoin=3, mgp=c(2,0.5,0))
plot.new()
for(j in 1:22){
	plot(NA,NA, xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
	text(.5,.5,j, srt=45)
}
for(i in c(levels(tumourTypes), "Pan-Cancer")){
	plot(NA,NA, xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0,2), ylim=c(0,1))
	text(2,.5,i, pos=2)
	x <- pi0[[i]]
	x <- x[x$type=="SingleGain" & ! x$Sample %in% names(which(isWgd)),]
	if(i == "Pan-Cancer") x <- pi0_agg
#	e <- try({
#				CHR <- factor(x$chr, levels=1:22)
#				l <- lm(car::logit(pi0) ~ CHR, data=x)
#				s <- summary(l)
#			})
	e <- try({
				lmmFit <- lmer(car::logit(pi0) ~ (1 | chr), data=x,  1/(car::logit(x$uCI) - car::logit(x$lCI) + 0.5) )
			})
	for(j in 1:22){
		plot(NA,NA, yaxt="n", xaxt="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
		if(j==1)
			axis(side=2, at=c(0,1), labels = i == levels(tumourTypes[1]), las=2)
		w <- x$chr==j 
		if(sum(w)==0) next
		par(xpd=TRUE)
		abline(h=c(0.5), col='grey',lty=3)
		par(xpd=NA)
		z <- boxplot(x$pi0[w],plot=FALSE)
		z$stats[,] <- quantile(x$pi0[w], c(0.05,0.25,0.5,0.75,0.95), na.rm=TRUE)
		t1 <- 1.05
		t2 <- 1.5
		l <- 30
		c <- colorRampPalette(rev(c("#990033","#FF6666","#CCCCCC","#3399FF","#006699")))(l)#RColorBrewer::brewer.pal(5,"RdBu")
		#c[3] <- "grey"
		#f <- sapply(1:2, function(i) {if(all(is.na(z$conf[,i]))) NA else if(all(z$conf[,i] > t2)) 1 else if(all(z$conf[,i] > t1)) 2 else if(all(z$conf[,i] < 1/t2)) 5 else if(all(z$conf[,i] < 1/t1)) 4 else 3})
		#f <- c[f]
		v <- paste0("chr",j,":tumourType",i)
		#f <- c[3]
#		try({
#					if(q[v] < 0.1)
#						if(s$coefficients[v,1] > 1) f <- c[5] else f <- c[1]
#					else if(p[v] < 0.05)
#						if(s$coefficients[v,1] > -1) f <- c[4] else f <- c[2]
#				}, silent=TRUE)
#			
		if(class(e) != "try-error"){
			r <- coef(lmmFit)$chr
			rr <- cut(r$`(Intercept)`, seq(-3,3, l=l+1)/log(2))#c(-10,-2,-1,1,2,10)*log(2))
			if(paste(j) %in% rownames(r) & sum(w) >= 5){
				f <- c[rr[match(j, rownames(r))]]
				#if(r[paste(j),1] < -2) f <- c[1] else if(r[paste(j),1] > 2) f <- c[5]
			}
		}
#		if(class(e) != "try-error"){
#			if(v %in% rownames(s$coefficients))
#				f <- ifelse(s$coefficients[v,4] < 0.001, "black","grey")
#			else
#				f <- "grey"
#		} else f <- "grey"
		#if(f != c[3]) rect(0,0,1,1)
		bxp(z, pch=NA,staplewex=NA, lty=1, xaxt="n", yaxt="n", add=TRUE, xlab="", ylab="", boxfill=f, at=0.5, boxcol=f, medlwd=1,medcol="white", whiskcol=f, whisklwd=2, boxwex=0.2 + 8*0.75 * sum(w)/nrow(x), boxlty=1)	
	}
}	

dev.copy(device=pdf, file="Timing-overview-lmm_mean.pdf", width=4, height=6, pointsize=8)
dev.off()

x <- PCAWG_sigDecomp_merged[1,,]
x <- x[,colSums(x)!=0]
par(mar=c(4,4,1,1), mgp=c(3,.5,0), bty="n", las=2)
barplot(t(x), legend=TRUE, col=RColorBrewer::brewer.pal(6, "Set1"), ylab="Number of mutations", args.legend=list(x="topleft", title="Signature", bty="n"))
dev.copy(device=pdf, file=paste0(sampleIds[1],"sig.pdf"), width=3, height=3, pointsize=8)
dev.off()

dir.create("2016-12-02")
for(o in ls()) save(list=o, file=paste0("2016-12-02/",o,".RData"))


par(mfrow=c(4,4))
i <- j <- 1; while(i <= 16){
	v <- allVcf[[which(isWgd)[j]]]
	w <- which(info(v)$MJCN==2 & info(v)$MNCN==2 & sapply(info(v)$CNID, length)==1)
	m <- info(v)$MCN[w]
	t <- table(m, as.numeric(info(v)$CNID[w]))
	b <- allBB[[meta(header(v))["ID",]]]
	n <- paste0(seqnames(b), ":", round(start(b)/1e6), "-",round(end(b)/1e6),"M")
	colnames(t) <- n[as.numeric(colnames(t))]
	s <- colSums(t) > 100
	if(any(s)){
		plot(2*t[2,s]/(2*t[2,s]+t[1,s]), ylab="time", ylim=c(0,1), xaxt="n", xlab="", main=substr(meta(header(v))["ID",],1,6))
		axis(side=1, at=1:sum(s), labels=colnames(t)[s], las=2)
		i <- i+1
	}
	j <- j + 1
}

t1 <- table(unlist(sapply(allBB[isWgd], function(b){
			paste(b$major_cn, b$minor_cn, sep=":")
		})))

t2 <- table(unlist(sapply(allBB[!isWgd], function(b){
					paste(b$major_cn, b$minor_cn, sep=":")
				})))


allVcfNew <- vector(mode='list', length=length(sampleIds))
names(allVcfNew) <- sampleIds
for(s in sampleIds){
	e <- new.env()
	load(paste0("../dp/2016_04_02_vanloo_wedge_consensusSNV_battenbergCNA/5_annotated_vcf/",s,".merged.annocounts.nolowsup_conf.somatic.snv_mnv.complete_annotation.vcf.RData"), envir=e)
	allVcfNew[[s]] <- e$vcf
}

kld <- function(x,y) {
	x <- x/sum(x)
	y <- y/sum(y)
	sum(x * log(x/y))
}

confusion <- function(x,y){
	t <- table(x,y)
	sum(diag(t))/sum(t)
}

c <- sapply()

drivers <- read.table("../ref/pcawg_whitelist_coding_drivers_v1_sep302016.txt", header=TRUE, sep="\t")

getDriverVaf <- function(drivers)


d <- intersect(levels(drivers$gene),cancerGenesS)
r <- DNAStringSet(drivers$ref)
a <- DNAStringSet(drivers$alt)

finalDrivers <- VRanges(seqnames = drivers$chr, ranges=IRanges(drivers$pos, width =  width(r)), ref=r, alt=a, sampleNames  = drivers$sample_id)
mcols(finalDrivers) <- drivers[,-c(1,3,4,5,6)]
finalDrivers$MCN <- finalDrivers$MNCN <- finalDrivers$MJCN <- as.numeric(rep(NA,nrow(drivers)))
finalDrivers$CLS <- factor(NA, levels = levels(info(allVcfNew[[1]])$CLS))

for(s in levels(sampleNames(finalDrivers))){
	if(! s %in%  sampleIds) next
	for(x in c(allVcfNew[[s]], indelAnnotated[[s]])){
		w <- sampleNames(finalDrivers) == s
		f <- findOverlaps(finalDrivers[w] , x, select='first')
		mcols(finalDrivers)[w[!is.na(f)],c("MCN","MNCN","MJCN","CLS")] <- info(x)[f[!is.na(f)], c("MCN","MNCN","MJCN","CLS")]
	}
}


getDriverGenotype <- function(vcf, drivers, reclassify='missing', ...){
	id <- meta(header(vcf))["ID",1]
	cls <- classifyMutations(vcf = vcf, reclassify=reclassify)
	t <- info(vcf)$TCN
	if(is.null(t))
		t <- info(vcf)$MNCN + info(vcf)$MJCN
	hom <- factor(info(vcf)$MCN==t, levels=c(TRUE,FALSE))
	table(gene=factor(unlist(info(vcf)$DG), levels=as.character(cancerGenes)), class=cls, homozygous=hom, ...)
}


#' ## Broad 500
dpPath <- '/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/broad500/Subclonal_Structure'
dpFiles <- dir(dpPath, pattern="subclonal_structure.txt", recursive=TRUE)
simPurityPloidy <- read.table("/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/broad500/pp_table.txt", header=TRUE, sep='\t')
rownames(simPurityPloidy) <- simPurityPloidy$sample
simPurityPloidy <- simPurityPloidy[,2:3]

prec <- simplify2array(mclapply(dir("../broad500", pattern="annotated_vcf", full.names=TRUE)[-3], function(d){
			x <- sapply(dir(d, pattern=".vcf.RData", full.names=TRUE), function(f){
						e <- new.env()
						load(f, envir=e)
						ID <- meta(header(e$vcf))["ID",1]
						truth <- read.table(paste0("../broad500/Mut_Assign/",ID,".mutation_assignments.txt"), header=TRUE, sep="\t")
						truth <- GRanges(truth$chr, IRanges(truth$pos, width=1), cluster=truth$cluster)
						seqlevels(truth) <-paste( 1:22)
						truth <- sort(truth)
						purity <- simPurityPloidy[ID,"purity"]						
						clusters <- loadClusters(ID)
						clusters$proportion <- clusters$ccf * purity
						x <- mean(info(e$vcf)$CNF == clusters$proportion[truth$cluster+1])
						y <- mean((info(e$vcf)$CNF==purity) == (clusters$proportion[truth$cluster+1]==purity))
						return(c(prec.all=x,prec.clonal=y))
					})
			colnames(x) <- sub("(.+)(Sim_.+)(.no_rea.+)","\\2",colnames(x))
			return(x)
		}, mc.cores=6))
dimnames(prec)[[3]] <- sub("annotated_vcf_*", "",dir("../broad500", pattern="annotated_vcf"))
#dimnames(prec)[[3]][1] <- "betabin_xmin3"
#dimnames(prec)[[3]][4] <- "betabin_xmin0"

prec[[3]] <- prec[[3]][,match(colnames(prec[[2]]),colnames(prec[[3]]))]
prec <- simplify2array(prec)

boxplot(prec[2,,])

load("../broad500/annotated_vcf_rho0_xmin3_local/Sim_500_209.no_real_info.complete_annotation.bb_granges.RData")
bb_local <- bb
load("../broad500/annotated_vcf_rho0_xmin3_global/Sim_500_209.no_real_info.complete_annotation.bb_granges.RData")
bb_global <- bb

p <- sapply(c("bb_local","bb_global"), function(b){
			unlist(sapply(1:length(bb_local), function(i){
			b <- get(b)
			t <- b$timing_param[[i]]
			try(return(t[,"power.s"] * t[,"power.m.s"]))
			try(return(t[,"power"]))
		}))})

load("../broad500/annotated_vcf_rho0_xmin3_global/Sim_500_209.no_real_info.complete_annotation.vcf.RData")

simClust <- sapply(rownames(simPurityPloidy), function(i){
			cl <- loadClusters(i)
			cl$proportion <- cl$ccf * simPurityPloidy[i,'purity']
			return(cl)
		}, simplify=FALSE)


precClip <- sapply(dir("../broad500/clip/CliP", pattern="mutation_assignment.txt.gz", full.names=TRUE), function(f){
			ass <- read.table(gzfile(f), header=TRUE, sep="\t")
			ar <- GRanges(ass$chr, IRanges(ass$pos, width=1))
			cl <- read.table(gzfile(sub("mutation_assignment","subclonal_structure",f)), header=TRUE, sep="\t")
			cc <- factor(ass$cluster, labels=cl$proportion)
			cf <- as.numeric(as.character(cc))
			ID <- sub("(.+)(Sim_.+)(_mut.+)","\\2", f)
			truth <- read.table(paste0("../broad500/Mut_Assign/",ID,".mutation_assignments.txt"), header=TRUE, sep="\t")
			truth <- GRanges(truth$chr, IRanges(truth$pos, width=1), cluster=truth$cluster)
			seqlevels(truth) <-paste( 1:22)
			truth <- sort(truth)
			o <- findOverlaps(ar, truth, select='first')
			purity <- simPurityPloidy[ID,"purity"]						
			clusters <- loadClusters(ID)
			clusters$proportion <- clusters$ccf * purity
			x <- mean(abs(cf- clusters$proportion[truth$cluster[o]+1]) < 0.01)
			y <- mean((cf==max(cf)) == (clusters$proportion[truth$cluster[o]+1]==purity))
			return(c(prec.all=x,prec.clonal=y))
		})
colnames(x) <- sub("(.+)(Sim_.+)(.no_rea.+)","\\2",colnames(x))





ID <- grep("0176cf1d", sampleIds, value=TRUE)
par(mar=c(3,3,1,1), bty="L", mgp=c(2,.5,0))
bbplot(finalBB[[ID]], ylim=c(-0.25,5))
legend("topright", col=c("black",RColorBrewer::brewer.pal(4,"Set1")[c(1,2)]), lty=1, lwd=4, c("Total","Major","Minor"))
dev.copy2pdf(file=paste0(ID,".CN.pdf"), width=16, height=4)

applyPigeonHole
function(ID){
	c <- loadClusters(ID)
	p <- purityPloidy[ID,"purity"]
	mcf <- c$proportion#*p
	l <- sapply(1:length(mcf), function(i) mcf[i] > pmax(0,1-mcf))
	w <- which(l & upper.tri(l), arr.ind=TRUE)
	cbind(c$cluster[w[,1]], c$cluster[w[,2]])
}


applyPigeonHole <- function(clusters){
	mcf <- clusters$proportion
	l <- sapply(1:length(mcf), function(i) mcf[i] > pmax(0,1-mcf))
	a <- l & upper.tri(l)
	w <- which(a, arr.ind=TRUE)
	cbind(clusters[w[,1],1], clusters[w[,2],1])
}



#' ### CN-LOH cases
w <- which(finalPloidy < 2.6 & aep > .5)
w

#ID <- "20d1b88b-3ff6-4201-a748-6a993500c652"
par(mfrow=c(10,2),mar=c(3,3,2,1)+.1, cex=.5, mgp=c(2,.5,0))
for(ID in names(w)){
	t <- bbToTime(finalBB[[ID]])
	bbplot(finalBB[[ID]])
	title(main=ID)
	n <- countQueryHits(findOverlaps(finalBB[[ID]], finalSnv[[ID]]))
	i <- t[,"type"]=="cnloh" & n >= 50
	if(sum(i, na.rm=TRUE) > 0)
		plot(t[i,"time"], ylim=c(0,1), ylab="time", xlab="segment")
	else plot.new()
}

dev.copy2pdf(file="CNLOH.pdf", width=8, height=10, pointsize=12)


w <- c(which(isWgd)[1:13], grep("c9ad6b1|ec3998|1dbdbb", names(finalBB)))
par(mfrow=c(4,4),mar=c(3,3,2,1)+.1, cex=.5, mgp=c(2,.5,0), xpd=FALSE, las=2)
for(ID in w){
	n <- countQueryHits(findOverlaps(finalBB[[ID]], finalSnv[[ID]]))
	t <- bbToTime(finalBB[[ID]])
	l <- width(finalBB[[ID]])
	i <- n>100 & l > 10e6
	plot(NA,NA, xlab='Position', ylab="Time", ylim=c(0,1), xlim=c(0,chrOffset["MT"]))
	abline(v = chrOffset[1:24], lty=3)
	if(length(which(i)>0))
		for(j in which(i)){
			s <- start(finalBB[[ID]])[j]
			e <- end(finalBB[[ID]])[j]
			x <- chrOffset[as.character(seqnames(finalBB[[ID]])[j])]
			segments(s+x,t[j,"time"], e+x, t[j,"time"], pch=19,col=as.numeric(factor(t[,"type"]))[j], ylim=c(0,1), lwd=2)
		}
	title(main=names(finalBB)[ID], line=0)
	if(ID==1) legend("topright",col=1:3, lwd=2, legend=levels(factor(t[,"type"])))
}



.wgdTime <- function(vcf, bb, clusters, purity, n.boot=200){
	wgd <- info(vcf)$MajCN==2 & info(vcf)$MinCN==2& sapply(info(vcf)$CNID, length)==1
	cnloh <- info(vcf)$MajCN==2 & info(vcf)$MinCN==0& sapply(info(vcf)$CNID, length)==1
	single <- info(vcf)$MajCN==2 & info(vcf)$MinCN==1& sapply(info(vcf)$CNID, length)==1
	w <- which(wgd|single|cnloh)
	v <- vcf[w]
	seqnames(rowRanges(v)[which(wgd[w])]) <- factor(rep(4, sum(wgd, na.rm=TRUE)), levels=seqlevels(v))
	seqnames(rowRanges(v)[which(cnloh[w])]) <- factor(rep(2, sum(cnloh, na.rm=TRUE)), levels=seqlevels(v))
	seqnames(rowRanges(v)[which(single[w])]) <- factor(rep(3, sum(single, na.rm=TRUE)), levels=seqlevels(v))

	b <- c(GRanges(factor(4, levels=2:4), IRanges(1,max(end(v))), copy_number=4, major_cn=2, minor_cn=2, clonal_frequency=purity),
			GRanges(factor(2, levels=2:4), IRanges(1,max(end(v))), copy_number=2, major_cn=2, minor_cn=0, clonal_frequency=purity),
			GRanges(factor(3, levels= 2:4), IRanges(1,max(end(v))), copy_number=3, major_cn=2, minor_cn=1, clonal_frequency=purity))
			
	l <- computeMutCn(v, b, clusters, purity, isWgd=TRUE, n.boot=n.boot)
	b$timing_param <- l$P
	t <- bbToTime(b)
	data.frame(seqnames=NA, start=NA, end=NA, n=c(sum(wgd, na.rm=TRUE),sum(cnloh, na.rm=TRUE),sum(single, na.rm=TRUE) ), t, sample=meta(header(vcf))['ID',1])
}


allSegments <- mclapply(names(finalBB), function(ID){try({
			#if(!isWgd[ID]){
				n <- countQueryHits(findOverlaps(finalBB[[ID]], finalSnv[[ID]]))
				t <- bbToTime(finalBB[[ID]])
				w <- !is.na(t[,2])
				data.frame(as.data.frame(granges(finalBB[[ID]]))[w,1:3], n=n[w], t[w,], sample=rep(ID, sum(w)) )
			#}else{
			#	.wgdTime(finalVcfSnv[[ID]], finalBB[[ID]], finalClusters[[ID]], purityPloidy[ID,1])
			#}
					})
		}, mc.cores=6)

finalSegments <- do.call("rbind", allSegments[sapply(allSegments, class)!='try-error'])
w <- finalSegments$sample %in% names(which(isWgd))
finalSegments$seqnames[w] <- finalSegments$start[w] <- finalSegments$end[w] <- NA

weighted.sd <- function(x, w){ m <- sum(x*w)/sum(w); sqrt(sum(w * (x-m)^2/sum(w)) * length(x)/length(x[-1]))}
finalSegmentsAggregated <- do.call("rbind",lapply(split(finalSegments, paste(finalSegments$sample, finalSegments$seqnames, finalSegments$type)), function(x) data.frame(chr=x$seqnames[1], start=min(x$start), end=max(x$end), width=sum(x$end-x$start), n=sum(x$n), time=sum(x$time * x$n)/sum(x$n), sd.time=weighted.sd(x$time, x$n), type=x$type[1], sample=x$sample[1])))

write.table(finalSegments, file="2017-01-27-segments-time.txt", col.names=TRUE, sep='\t', quote=FALSE, row.names=FALSE)

write.table(finalSegmentsAggregated, file="2017-02-01-segmentsAggregated-time.txt", col.names=TRUE, sep='\t', quote=FALSE, row.names=FALSE)

#' WGD cases
w <- finalSegmentsAggregated[is.na(finalSegmentsAggregated$chr),c("type","time","n","sample","sd.time")]
n <- names(which(isWgd))
s <- split(w[,c("time","sample")], w$type)
x <- sapply(s, function(x) x$time[match(n,x$sample)])
plot(as.data.frame(x))

averageCNLOH <- function(bb){
	sum(width(bb) * (bb$major_cn == 2 & bb$minor_cn == 0) * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
}



averageStar <- sapply(finalBB, function(bb) sum(width(bb) * bb$star * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE))
aStar <- cut(averageStar, 0:3)

ac <- sapply(finalBB, averageCNLOH)
finalHom <- sapply(finalBB, averageHom)
names(finalHom) <- names(finalBB)

plot(finalHom, finalPloidy, col=.classWgd( finalPloidy, finalHom)+1, xlim=c(0,1))

x <- seq(0,1,0.01)
lines(x, pmax(2, 3-2*x))

write.table(file="HomPloidy,txt",data.frame(homozygousity=finalHom, ploidy=finalPloidy), col.names=TRUE, quote=FALSE, sep="\t")

d <- rowMeans(apply(t[,2:4],1, function(x){
					b <- timeToBeta(x)
					dbeta(seq(0,1,0.01), b[1], b[2])
				}), na.rm=TRUE)

#' ## WGD-possible
w <- which(finalPloidy < 2.8 & finalHom > .3 & .classWgd( finalPloidy, finalHom))


fractionGenomeAmplified <- function(bb){
	x <- seq(0.01,0.99,0.01)
	m <- matrix(rep(x, length(bb)), nrow=length(bb), byrow=TRUE)
	time.lo <- matrix(rep(bb$time.lo, each=length(x)), ncol=length(x), byrow=TRUE)
	time.up <- matrix(rep(bb$time.up, each=length(x)), ncol=length(x), byrow=TRUE)
	i <- m >= time.lo & m <= time.up
	w <- as.numeric(width(bb))
	colSums(i * w, na.rm=TRUE)/sum(w*bb$clonal_frequency)
}


mclapply(names(finalBB), function(ID){
	#print(ID)
	pdf(paste0("../final/annotated_006/timing/",ID,".pdf"), 12,6.5)
	bb <- finalBB[[ID]]
	t <- bbToTime(bb)	
	par(mfrow=c(2,1), mar=c(2,3,2,1), mgp=c(2,.5,0), bty="L", cex=1, las=2)
	plotBB(bb, ylim=c(0,8))
	title(main=paste0(ID,", ", donor2type[sample2donor[ID]], ", ploidy=",round(finalPloidy[ID],2), ", hom=",round(finalHom[ID],2), if(.classWgd( finalPloidy[ID], finalHom[ID])) ", WGD" else ""), font.main=1, line=0)
	par(mar=c(3,3,1,1))
	plotTiming(bb, t)
	dev.off()
}, mc.cores=6)

d <- data.frame(as.data.frame(granges(finalBB[[grep("^ec39",names(finalBB))]])), bbToTime(finalBB[[grep("^ec39",names(finalBB))]]))

ID <- grep('^ec39', names(finalSnv), value=TRUE)
v <- finalSnv[[ID]]
v <- v[seqnames(v)=="15"]
b <- finalBB[[ID]]
b <- b[seqnames(b)=="15"]

boot <- sapply(1:10, function(foo){
			set.seed(foo)
			vv <- v[sort(sample(1:nrow(v), replace=TRUE))]
			b$timing_param <- NULL
			c <- computeMutCn(vv,b,finalClusters[[ID]], purity=purityPloidy[ID,1], isWgd=TRUE)
			d <- b
			d$timing_param <- c$P
			bbToTime(d)$time
		})
boxplot(t(bo))

b0 <- GRanges(seqnames=15, IRanges(15800000, 32090715), total_cn=2, major_cn=2, minor_cn=0, clonal_frequency=0.34)

boot.merged <- sapply(1:50, function(foo){
			set.seed(foo)
			vv <- v[sort(sample(1:nrow(v), replace=TRUE))]
			c <- computeMutCn(vv,b0,finalClusters[[ID]], purity=purityPloidy[ID,1], isWgd=TRUE)
			d <- b0
			d$timing_param <- c$P
			bbToTime(d)$time
		})

cnStates <- b$timing_param[[2]]
whichStates <- 1:nrow(cnStates)
w <- v %over% b[3]
x <- getAltCount(v[w])
n <- getTumorDepth(v[w])
L <- matrix(sapply(pmin(cnStates[whichStates,"f"],1), function(pp) dtrbetabinom(x,n,pp, rho=rho, xmin=pmin(x,xmin)) + .Machine$double.eps), ncol=length(whichStates))
L <- rbind(L, rep(1,ncol(L)))

pi <- cnStates[,"P.m.sX"]
N <- as.numeric(L %*% pi)
H <- t(L/N) %*% (L/N) 
solve(H)

bi <- sapply(1:100, function(foo){
			LL <- L[sample(1:nrow(L), replace=TRUE),]
			p <- pi
			for(i in 1:100){
				P <- LL*rep(p / cnStates[,"power.m.s"], each=nrow(LL))
				P <- P/rowSums(P)
				p <- colMeans(P)
			}
			return(p)
		})


#' WG output
for(ID in names(finalSnv)){
	f <- factor(round(info(finalSnv[[ID]])$CNF,3))
	t <- data.frame(chr=as.character(seqnames(finalSnv[[ID]])), pos=start(finalSnv[[ID]]), cluster=as.numeric(f), clonal=info(finalSnv[[ID]])$CLS)
	write.table(file=paste0("../final/annotated_006/output/2_subclones/",ID,"_mutation_assignments.txt"), t, quote=FALSE, sep="\t")
	t <- as.data.frame(table(f))
	write.table(file=paste0("../final/annotated_006/output/2_subclones/",ID,"_subclonal_structure.txt"), data.frame(cluster=1:nrow(t), n_ssms=t$Freq, proportion=t$f),quote=FALSE, sep="\t")
}
write.table(file="../final/annotated_006/output/1_purity_ploidy/purity_ploidy.txt" ,data.frame(purity=finalPurity, ploidy=finalPloidy))


inputClusters <- sapply(dir(dpPath, pattern="subclonal_structure.txt", full.names=TRUE), read.table, header=TRUE, sep="\t", simplify=FALSE)
names(inputClusters) <- sub("_sub.+","", gsub(".+/","",names(inputClusters)))

cl <- sapply(inputClusters, removeSuperclones, simplify=FALSE)
ml <- sapply(cl, mergeClusters, deltaFreq=0.1, simplify=FALSE)
p <- sapply(cl, function(x) x$proportion[nrow(x)])
r <- sapply(ml, function(x) x$proportion[nrow(x)])
d <- sapply(cl, function(x) if(nrow(x)>1) (x$proportion[nrow(x)] - x$proportion[nrow(x)-1])/x$proportion[nrow(x)] else NA)
a <- sapply(cl, function(x) if(nrow(x)>1) (x$n_ssms[nrow(x)] - x$n_ssms[nrow(x)-1])/sum(x$n_ssms) else NA)

s <- apply(purityPloidy[,-2:-1],1,sd, na.rm=TRUE)

q <- purityPloidy[names(cl),2] < 2.05 & purityPloidy[names(cl),2] > 1.9

plot(purityPloidy[names(cl),1], p, col=q + 1, pch= (s[names(cl)] > 0.1)+1)

w <- abs(purityPloidy[names(cl),1]-p)>0.1




fracGenomeWgdComp <- t(sapply(finalBB[wgdPoss], function(bb) {
					fgw <- try(fractionGenomeWgdCompatible(bb)); 
					if(class(fgw)!='try-error') fgw
					else rep(NA,6)}))
rownames(fracGenomeWgdComp) <- names(finalBB)[wgdPoss]

wgdPossStar <- factor(rep(2,sum(wgdPoss)), levels=0:3, labels=c("unlikely","uninformative","likely","very likely"))
wgdPossStar[fracGenomeWgdComp[,"avg.ci"]>0.75 | fracGenomeWgdComp[,"nt.total"]/chrOffset["MT"] < 0.33 ] <- "uninformative"
wgdPossStar[fracGenomeWgdComp[,"nt.wgd"]/fracGenomeWgdComp[,"nt.total"] < 0.66] <- "unlikely"
wgdPossStar[wgdPossStar=="likely" & fracGenomeWgdComp[,"nt.wgd"]/fracGenomeWgdComp[,"nt.total"] > 0.8 & fracGenomeWgdComp[,"sd.wgd"] < 0.1 &  fracGenomeWgdComp[,"nt.total"]/chrOffset["MT"] > 0.5] <- "very likely"
names(wgdPossStar) <-  names(finalBB)[wgdPoss]
prop.table(table(wgdPossStar))

pdf("WGDposs2-timing.pdf", 12,6.5) 
j <- 1
for(ID in names(finalBB)[wgdPoss][m < 0.5]){
	if(j%%100 == 0) print(j); j <- j+1
#	for(f in dir("../final/annotated_007/snv_mnv", pattern=paste0(ID,".+bb_granges.RData"), full.names=TRUE)) load(f)
#	t <- bbToTime(bb)
	par(mfrow=c(2,1), mar=c(3,3,2,1), mgp=c(2,.5,0), bty="L", cex=1, las=2)
	plotBB(finalBB[[ID]], ylim=c(0,8))
	title(main=paste0(ID,", ", donor2type[sample2donor[ID]], ", ploidy=",round(finalPloidy[ID],2), ", hom=",round(finalHom[ID],2)), font.main=1, line=0)
	par(mar=c(3,3,2,1))
	plotTiming(finalBB[[ID]])
	abline(h=fracGenomeWgdComp[ID,"time.wgd"], lty=3)
	title(main=paste0("Timeable=", round(fracGenomeWgdComp[ID,2]/chrOffset["MT"]*100), "%, WGD=",round(fracGenomeWgdComp[ID,1]/fracGenomeWgdComp[ID,2]*100), "%, sd.WGD=",round(fracGenomeWgdComp[ID,'sd.wgd'],2), "%, avg.CI=",round(fracGenomeWgdComp[ID,'avg.ci'],2),", verdict=",wgdStar[ID]), font.main=1, line=0)
}
dev.off()

greylist <- "0ef92ff8−829f−425a−91d8−c594b6e22a2b
1021b60d−f7b2−43b0−b2cc−f282d619d533
14b8bbf2−310b−459b−b52d−a7ef510ce1cf
3a6bce45−0431−49d8−82df−b1d9a738e5a6
3db3b7b1−da1d−4b9c−a92a−c60fecf4328c 
46d35b82−e1b7−4d35−be5e−3a70fd47e421 
6dbac113−d4cf−4db5−97c9−50fa400bb47e
8454fe53−869d−41c8−b0c8−a7929d00eec3  
b35d9a68−29f4−49ab−b83e−b5151679e3af
cd9efdef−a7fb−49e5−9515−63606ae8bbfc
e93b0979−65ef−4883−9b6e−39eb17966e66
f94c4f69−8119−4eaf−97c1−5106890c14d4"
greylist <- sub(" +","",gsub("−","-",strsplit(greylist,"\n")[[1]]))

table(wgdStar[greylist])

ID <- names(finalSnv)[1]

ig <- c('0448206f-3ade-4087-b1a9-4fb2d14e1367',
'062e96d4-c623-11e3-bf01-24c6515278c0',
'0bfd1068-3fcf-a95b-e050-11ac0c4860c3',
'1bea3a72-3b73-4072-a6bb-96a90119d3ac',
'28e81540-4744-4865-b627-c7c9d8a3c2b8',
'35553150-e4ef-4539-b220-259f2d634bd7',
'3b41cb48-c623-11e3-bf01-24c6515278c0',
'3bacc189-01b8-46cc-a442-f393c0f428c6',
'43dadc68-c623-11e3-bf01-24c6515278c0',
'44083f54-0953-48e3-a704-11ad0988ad2e',
'59e2d6d1-debd-4796-ab0c-6a5673a990fc',
'5d6ad982-bb01-4233-b8fa-d129460eec79',
'63762458-902a-4329-a823-703b54cb5f9d',
'68509ede-3dcf-4a6e-9af0-4a9bb4dfa567',
'6d936ef9-b5df-44d3-831f-528bf8ddc131',
'7ccb9a4d-6f48-41c2-a630-27fde8c67d60',
'83a1b304-2ec1-44ae-a9c5-8ad3a2a46a1f',
'8dd14f0e-8601-4aa1-864c-3c49e768cdd1',
'9321341c-c622-11e3-bf01-24c6515278c0',
'9536f736-63bc-4099-bd54-740f5910f4a8',
'9aac83e4-c622-11e3-bf01-24c6515278c0',
'b47aa163-eec9-4225-940b-4373e78152e2',
'b7a7d93b-38a7-4fc3-a433-3bb0a8cb7c42',
'b994762c-c622-11e3-bf01-24c6515278c0',
'c285c2fa-24b4-47a1-874d-86e74b002b05',
'c95a2b1b-726c-4608-9fff-d57b6f1aa75a',
'd05ea63c-86a3-463a-a790-2edaa74b4da7',
'e5593865-5f8e-4a4c-b36f-73fbe64d66da',
'f064f762-c622-11e3-bf01-24c6515278c0',
'f393bb0a-9b20-a0e5-e040-11ac0d48454e',
'fc447d4f-2532-c8ea-e040-11ac0c48469f',
'fe63d42b-d471-45b6-9bdf-1a3b55465d37',
'ffe4bb51-e98a-41a7-a4e1-c3970386889c',
'5ab6a1d3-76f8-45d4-a430-d9831daa9ec4')

pdf("WGD-45poss-timing.pdf", 12,6.5) 
j <- 1
for(ID in intersect(ig, names(finalBB))){
	if(j%%100 == 0) print(j); j <- j+1
#	for(f in dir("../final/annotated_007/snv_mnv", pattern=paste0(ID,".+bb_granges.RData"), full.names=TRUE)) load(f)
#	t <- bbToTime(bb)
	par(mfrow=c(2,1), mar=c(3,3,2,1), mgp=c(2,.5,0), bty="L", cex=1, las=2)
	plotBB(finalBB[[ID]], ylim=c(0,8))
	title(main=paste0(ID,", ", donor2type[sample2donor[ID]], ", ploidy=",round(finalPloidy[ID],2), ", hom=",round(finalHom[ID],2)), font.main=1, line=0)
	par(mar=c(3,3,2,1))
	plotTiming(finalBB[[ID]])
	abline(h=fracGenomeWgdComp[ID,"time.wgd"], lty=3)
	title(main=paste0("Timeable=", round(fracGenomeWgdComp[ID,2]/chrOffset["MT"]*100), "%, WGD=",round(fracGenomeWgdComp[ID,1]/fracGenomeWgdComp[ID,2]*100), "%, sd.WGD=",round(fracGenomeWgdComp[ID,'sd.wgd'],2), "%, avg.CI=",round(fracGenomeWgdComp[ID,'avg.ci'],2),", verdict=",wgdStar[ID]), font.main=1, line=0)
}

dev.off()

wgdTimeDeam <- do.call("rbind",mclapply(names(finalBB)[isWgd], function(ID){
					bb <- finalBB[[ID]]
					t <- .wgdTime(finalSnv[[ID]][isDeamination(finalSnv[[ID]])], bb[which((bb$time.lo + bb$time)/2 <= fracGenomeWgdComp[ID,'time.wgd'] & (bb$time.up + bb$time)/2 >= fracGenomeWgdComp[ID,'time.wgd'])], cluster=finalClusters[[ID]], purity=finalPurity[ID], n.boot=10)
					t$time.up <- (5 + t$n * t$time.up)/(5 + t$n)
					t$time.lo <- (0 + t$n * t$time.lo)/(5 + t$n)
					return(t)
}, mc.cores=6))

s <- as.data.frame(split(wgdTimeDeam$time, wgdTimeDeam$type))
plot(s)

library(ape)

bb <- finalBB[[1]]
trees <- sapply(1:10, function(foo){
			w <- !is.na(bb$time)
			t <- pmin(1, pmax(0,rnorm(bb$time[w], sd = (bb$time.up[w]- bb$time.lo[w])/4)))
			h <- hclust(dist(t))
			p <- as.phylo(h)
		}, simplify=FALSE)


approxCi <- function(bb, pseudo.count = 5, min.alpha = 1e-2, min.beta=min.alpha){
	time.ci <- sapply(c(0.025, 0.975), function(q) qbeta(q, pmax(min.alpha, bb$n.snv_mnv/4 * bb$time - pseudo.count), pmax(min.beta, bb$n.snv_mnv/4 * (1-bb$time) - pseudo.count)))
	return(time.ci)
}

b$timing_param <- NULL
L <- computeMutCn(finalSnv[[ID]], b, clusters=finalClusters[[ID]], purity=finalPurity[ID], xmin=3, gender="female", isWgd=TRUE, n.boot=10)
b$timing_param <- L$P
t <- bbToTime(b)
mcols(b)[names(t)] <- DataFrame(t)

par(mfrow=c(2,2))
plotBB(bb, ylim=c(0,8))
plotTiming(bb)
plot(bb$time, b$time)
plotTiming(b)


timing_output <- do.call("rbind", sapply(names(finalBB), function(ID){
					bb <- finalBB[[ID]]
					d <- as.data.frame(bb[, c("total_cn","major_cn","minor_cn","clonal_frequency","time", "time.lo","time.up","time.star","n.snv_mnv")])
					d$sample <- ID
					return(d)
				}, simplify=FALSE))

timing_output <- timing_output[!is.na(timing_output$time),]

write.table(timing_output, file=paste0("TimingAllSegmentsStar-",Sys.Date(), ".txt"), row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

ssnmf <- function(D, S, s=0, maxIter=100, minE = 0, whichS = 1:ncol(D)){
	n <- nrow(D)
	o <- ncol(S)
	m <- ncol(D)
	P <- cbind(S, if(s>0) matrix(runif(n*s,0,1), ncol=s) else NULL)
	P <- P/rep(colSums(P), each=nrow(P))
	colnames(P) <- c(colnames(S), if(s >0) paste0("N.",1:s) else NULL)
	E <- matrix(runif((s+o)*m, 0,1), ncol=m)
	E <- E * (t(P)%*% (D / (P %*% E))) / rep(colSums(P), m)
	
	iter <- 1
	while(iter < maxIter){
		P <- P * ((D / (P %*% E)) %*% t(E)) / rep(rowSums(E), each=n)
		if(o>0)	P[,1:o] <- S
		P <- P/rep(colSums(P), each=nrow(P))
		E <- E * (t(P)%*% (D / (P %*% E))) / rep(colSums(P), m)
		E[E/rep(colSums(E), each=nrow(E)) < minE] <- 0
		if(o>0) E[-(1:o),setdiff(1:ncol(E), whichS)] <- 0
		E <- E * rep(colSums(D)/colSums(E), each=nrow(E))
		iter <- iter + 1
	}
	list(E=E,S=P)
}

r <- ssnmf(S, S=matrix(0, ncol=0, nrow=96), s=5, maxIter=1000)
par(mfrow=c(5,1))
for(i in 1:5)  barplot(r$S[,i], ylim=c(0,0.1))

par(mfrow=c(1,1))
plot(S[,1], r$S %*% r$E[,1])



S <- sigTable[,donor2type[sample2donor[names(finalSnv)]]=="Breast-AdenoCA",]
S <- asum(S,3)
write.table(S, file="96-BrCA.txt", row.names=FALSE, col.names=FALSE)

write.table(t(asum(sigTable,3)), file="96-All.txt", row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(data.frame(sample=names(finalSnv), type=donor2type[sample2donor[names(finalSnv)]]), file="TumourTypes.txt", quote=FALSE)

n <- sapply(bb$timing_param, function(x) {t <- try(nrow(x)); if(is.null(t)) NA else t})
nMutCN <- n[sapply(info(vcf)$CNID,function(x) {t <- try(x[1]);if(class(t)=='try-error') NA else t})]

par(mfrow=c(2,2))
plot(seq(0,1,l=nrow(vcf)+1)[-1], sort(info(vcf)$pMutCNTail, na.last=TRUE),log='xy', ylim=c(1e-4,1), xlab="expected",ylab="pMutCNTail")
abline(0,1, col='red')

qqnorm(qnorm(info(vcf)$pMutCNTail))
abline(0,1, col='red')


plot(info(vcf)$pMutCNTail, xlab="Position", ylab="pMutCNTail" )
plot(qnorm(info(vcf)$pMutCNTail), xlab="Position", ylab="pMutCNTail (normal trnsf'd" )

resNorm95 <- sapply(finalSnv, function(vcf) quantile(qnorm(info(vcf)$pMutCNTail), c(0.025, 0.975), na.rm=TRUE))

frac2SD <- sapply(finalSnv, function(vcf) mean(abs(qnorm(info(vcf)$pMutCNTail)) < 2 , na.rm=TRUE))



swt <- sapply(finalSnv[1:10], function(vcf) shapiro.test(qnorm(info(vcf)$pMutCNTail)))


par(mfrow=c(1,1))
boxplot(resNorm95[,2] -resNorm95[,1] ~ donor2type[sample2donor[names(finalSnv)]], las=2, ylab="Width of 95% of SN-Res")

par(mfrow=c(1,1))
boxplot(frac2SD ~ donor2type[sample2donor[names(finalSnv)]], las=2, ylab="Fraction of data inside theoretical 95% CI")
abline(h=0.95, lty=3)




frac2SD


pdf("0009b464-b376-4fbc-8a56-da538269a02f.timing.pdf", 4,1.25, pointsize=8)
par(mar=c(3,3,0.5,0.5), mgp=c(2,0.25,0), bty="L", las=2, tcl=-0.25)
plotTiming(finalBB[[1]])
dev.off()



r <- fracGenomeWgdComp[,"nt.wgd"] / fracGenomeWgdComp[, "nt.total"]
w <- which(fracGenomeWgdComp[,])



d <- data.frame(fracGenomeWgdComp, WGD=isWgd, WGD.star=wgdStar) 
colnames(d) <- c("nt.coamp", "nt.amp", "time.coamp", "sd.coamp", "avg.ci", "sd.all", "WGD", "WGD.star")
write.table(d, file="Coamp-timing.txt", col.names=TRUE, sep="\t", row.names=TRUE)

which(d$nt.coamp/chrOffset["MT"] > 0.05 & d$time.coamp < 0.02 & d$avg.ci < 0.3 & !isWgd)

save(list=setdiff(ls(), c("finalIndel","finalSnv")), file=paste0(Sys.Date(),"-temp.Rdata"))


### Consensus QC

dir <- "/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/sd11/icgc_pancan_full/consensus_clustering/201703_consensus_clustering/"
q5.c <- simplify2array(mclapply(dir(dir, pattern="output_(cicc|sc3|wm)", full.names=TRUE), function(d){
			f <- dir(d, pattern="*_assignment.RData", full.names=TRUE)
			sapply(f, function(x){
								e <- new.env()
								load(x, env=e)
								mean(e$MCN$D$pMutCNTail < 0.025 | e$MCN$D$pMutCNTail > 0.975, na.rm=TRUE)
		})}, mc.cores=3))
#save(q5, file="q5.RDataa")
load("q5.RData")
q5.c <- sapply(q5.c, function(x) {n <- sub("_assignment.RData","",sub("/.+/","",names(x))); names(x) <- n; x})
n <- sapply(q5.c, names)
u <- Reduce("intersect", n)

q5.c <- sapply(q5.c, function(x) x[u])
colnames(q5.c) <- c("cicc","sc3","wm")

boxplot(1-q5.c)

w <- which(q5.c[,] > 0.1, arr.ind=TRUE)
w <- q5.c[unique(w[,1]),]

r <- apply(apply(data.frame(mg=q5.m, q5.c[m,]), 1, rank), 1, table)

barplot(r, legend=TRUE)

#' Quick check
load("../final/annotated_010/snv_mnv/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20160830.somatic.snv_mnv.complete_annotation.vcf.RData")
mean(abs(0.5- info(vcf)$pMutCNTail) > 0.475 , na.rm=TRUE)

q5.wm <- sapply(finalSnvWm, function(vcf) mean(abs(0.5- info(vcf)$pMutCNTail) > 0.475 , na.rm=TRUE))
names(q5.wm) <- gsub(".+/","", names(q5.wm))
boxplot(1-q5.wm ~ donor2type[sample2donor[names(q5.wm)]], las=2, ylab="Fraction of data inside theoretical 95% CI")
abline(h=0.95, lty=3)

paste(paste0(sub("/.+/","",names(which(q5.wm > 0.1))),"*.pdf"), collapse=",")

q1.wm <- sapply(finalSnvWm, function(vcf) mean(abs(0.5- info(vcf)$pMutCNTail) > 0.495 , na.rm=TRUE))
names(q1.wm) <- gsub(".+/","", names(q1.wm))
boxplot(1-q1.wm ~ donor2type[sample2donor[names(q1.wm)]], las=2, ylab="Fraction of data inside theoretical 99% CI")
abline(h=0.99, lty=3)

paste(paste0(sub("/.+/","",names(which(q1.wm > 0.03 & q5.wm < 0.1))),"*.pdf"), collapse=",")

n <- sapply(finalSnv, nrow)
paste(paste0(sub("/.+/","",names(which(n > 1e5 & finalPurity[names(n)] > 0.5))),"*.pdf"), collapse=",")


u <- unique(finalBB[[4]])
u$timing_param <- NULL
u$clonal_frequency <- max(u$clonal_frequency, na.rm=TRUE)
L <- computeMutCn(finalSnv[[4]], u, clusters=finalClusters[[4]], purity=finalPurity[4], xmin=3, gender="female", isWgd=FALSE, n.boot=0)
p <- cbind(t(sapply(L$D$pAllSubclones, function(x) if(length(x)!=0) x else rep(NA, 2))), 1-L$D$pSub)
klaR::triplot(p)
colSums(p, na.rm=TRUE)
table(L$D$CNF)
finalClusters[[4]]

e <-  apply(p, 1, function(x) {t <- try(mg14::entropy(x)); if(class(t)=="try-error") NA else t/log(2)})
boxplot(e ~ L$D$CNF)

# Example with missing clusters
ID <- grep("^0fbd", names(finalSnv), value=TRUE)
cl <- read.table("/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/icgc_pancan_full/consensus_clustering/201703_consensus_clustering/consensus_clusters_sc3/0fbd94b1-bb34-4620-841b-861a0b5e0c12_subclonal_structure.txt", header=TRUE, sep="\t")
cl <- cl[order(cl$proportion),]
u <- unique(finalBB[[ID]])
u$timing_param <- NULL
u$clonal_frequency <- max(u$clonal_frequency, na.rm=TRUE)
L <- computeMutCn(finalSnv[[ID]], u, clusters=cl, purity=finalPurity[ID], xmin=3, gender="female", isWgd=isWgd[ID], n.boot=0, rho=0)
p <- cbind(t(sapply(L$D$pAllSubclones, function(x) if(length(x)!=0) x else rep(NA, 2))), 1-L$D$pSub)
klaR::triplot(p)
colSums(p, na.rm=TRUE)
table(L$D$CNF)

postFreq <- p %*% cl$proportion
q <- cumsum(cl$n_ssms/sum(cl$n_ssms))
cp <- cut(postFreq, c(0, quantile(postFreq, q, na.rm=TRUE)))
table(cp)

klaR::triplot(p, col=cp)

sapply(levels(cp), function(l) mean(p[cp==l,which(levels(cp)==l)], na.rm=TRUE))

w <- apply(p,1, function(x) {w <- which.max(x); if(length(w)==0) NA else w})
sapply(1:3, function(i) colMeans(p[w==i,i, drop=FALSE], na.rm=TRUE))

fClonal <- sapply(finalSnv, function(vcf) {t <- table(info(vcf)$CLS); sum(t[1:3]/sum(t))})
fClonalExpected <- sapply(finalClusters, function(x) x$n_ssms[which.max(x$proportion)]/sum(x$n_ssms))


# Stefan's examples
s <- gsub(" ","",strsplit("02706819-bcab-4c49-a569-a4a8c60db1c0 
1eb37b28-fac2-477a-88b3-e04291a07926 
27f87d1e-2c32-4beb-9677-62f7a286673d 
532259b8-c622-11e3-bf01-24c6515278c0 
558239c7-a160-4228-8fdf-a0a1d2f62133 
55c75a2a-f3d2-4469-9d23-604cf539d548 
5cbd429f-ffab-41ad-8016-422f1c922e99 
75ad15b9-8f9c-40c1-9ca6-1e8454fbd310","\n")[[1]])

for(ss in s){
	print(ss)
	print(read.table(paste0("/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/icgc_pancan_full/consensus_clustering/201703_consensus_clustering/consensus_clusters_wm/",ss,"_subclonal_structure.txt"), header=TRUE, sep="\t"))	
}


w <- which(donor2type[sample2donor[names(finalSnv)]]=="CNS-GBM")
which(sapply(finalBB[w], function(bb) any(seqnames(bb)=="20" & bb$clonal_frequency < max(bb$clonal_frequency, na.rm=TRUE) & width(bb) > 20e6 & bb$major_cn>1)))

i <- 1
x <- 0:40
plot(x, dbetabinom(x, 40, 0.25, rho=0), type='s')
rho <-  c(0, 0.01,0.05,0.1)
for(r in rho[-1]){
	i <- i+1
	lines(x, dbetabinom(x, 40, 0.25, rho=r), type='s',col=i)
}
legend("topright", legend=rho, lty=1, col=seq_along(rho))


#' Fix drivers
for(i in seq_along(finalSnv)){
	v <- addFinalDriver(finalSnv[[i]], finalDrivers)
	info(finalSnv[[i]])$DG <- info(v)$DG
	v <- addFinalDriver(finalIndel[[i]], finalDrivers)
	info(finalIndel[[i]])$DG <- info(v)$DG
	if(i%%10) print(i)
}

qSubclone <- sapply(timeSubclones, function(x) apply(x[,], 2, quantile, c(0,0.25,0.5,0.75,1), na.rm=TRUE), simplify='array')
qWgd <- sapply(timeWgd, function(x) apply(x[,"hat",], 2, quantile, c(0,0.25,0.5,0.75,1), na.rm=TRUE), simplify='array')

finalDriver <- asum(finalGenotypesP,2)
finalDriverType <- sapply(levels(donor2type), function(l) asum(finalDriver[,,donor2type[sample2donor[dimnames(finalDriver)[[3]]]]==l],3), simplify='array')

save(finalDriverType, file="finalDriverType.RData")

pdf("tumourtype.pdf", 20,1)
par(mar=c(0,0,0,0)+.1)
t <- table(droplevels(donor2type[sample2donor[names(finalPurity)]]))
t <- t[order(tolower(names(t)))]
image(x=c(0,cumsum(t)), y=c(0,1), z=matrix(seq_along(t),ncol=1), col=tissueColors[names(t)], xlab="", ylab="")
dev.off()

germline <- readVcf("../germline/deb9fbb6-656b-41ce-8299-554efc2379bd_het.vcf", genome="GRCh37")
g3 <- germline[seqnames(germline)==3]

sapply(0:9, function(i){
			ir <- IRanges(start=1 + i*1e7, end=1e7+i*1e7)
			b <- deepSNV::bam2R("../bam/tumour/deb9fbb6-656b-41ce-8299-554efc2379bd.bam", chr=3, start=start(ir), stop=end(ir))
			r <- ranges(g3)
			w <- width(sub(".+_","",names(r))) == 3
			r <- r[r %over% ir & w]
			bb <- b[start(r),]
			bb <- bb[,1:5] + bb[,12:16]
			bbb <- t(apply(bb, 1, function(x) sort(x, decreasing=TRUE)[1:2]))
			p <- pbinom(bbb[,2],rowSums(bbb), 0.5)
			l <- Rle(p < 1e-20)
		})

w <- donor2type[sample2donor[names(finalSnv)]] %in% c("ColoRect-AdenoCA", "Stomach-AdenoCA")

#w <- which(asum(finalGenotypes[grep("VHL", rownames(finalGenotypes)),,,,donor2type[sample2donor[names(finalSnv)]]=="Kidney-RCC",drop=FALSE], c(1,4))>=1, arr.ind=TRUE)
g <- "VHL"

g <- "Z95704.4"
t <- donor2type[sample2donor[names(finalSnv)]]%in% levels(donor2type) #=="Kidney-RCC"
a <- asum(finalGenotypes[grep(g, rownames(finalGenotypes)),,,,,drop=FALSE], c(1,4:4))
w <- which(a>=1 & rep(t, each=prod(dim(a))/dim(a)[length(dim(a))]), arr.ind=TRUE)
rownames(w) <- names(finalSnv)[w[,3]]
r <- GRangesList(apply(w, 1, function(x){
			if(x[1]==1)
				v <- finalSnv[[x[3]]]
			else
				v <- finalIndel[[x[3]]]
			vv <- v[grep(g,info(v)$DG)]
			r <- rowRanges(vv)
			mcols(r) <- DataFrame(alt.count=getAltCount(vv), dp=getTumorDepth(vv), pSub=info(vv)$pSub, row.names=NULL, status=info(vv)$CLS, MajCN=info(vv)$MajCN, MinCN=info(vv)$MinCN, pMutCNTail=info(vv)$pMutCNTail)
			
			r
		}))
d <- as.data.frame(r)
d


par(mfrow=c(3,1), mar=c(3,3,3,1), cex=.5, las=2)
boxplot(qnorm(t(finalGenotypesQ[rownames(g)[1:250],1,])), main="Subs", xaxt="n", ylab="Residual")
abline(h=c(-2,2), lty=3)
boxplot(qnorm(t(finalGenotypesQ[rownames(g)[1:250],2,])), main="Indel", ylab="Residual")
abline(h=c(-2,2), lty=3)


qEnt <- t(sapply(finalSnv, function(vcf){
			p <- as.matrix(info(vcf)[,c("pSingle","pGain","pSub")])
			e <- -rowSums(log(p^p))
			quantile(e, seq(0,1,0.1),na.rm=TRUE)
		}))

par(mfrow=c(1,1))
boxplot(qEnt[,"50%"]/log(2) ~ donor2type[sample2donor[names(finalSnv)]], las=2, ylab="Fraction of data inside theoretical 95% CI")
abline(h=0.95, lty=3)



fracGenomeWgdCompGray <- t(sapply(finalBBGray, function(bb) {
					fgw <- try(fractionGenomeWgdCompatible(bb)); 
					if(class(fgw)!='try-error') fgw
					else rep(NA,10)}))
rownames(fracGenomeWgdCompGray) <- names(finalBBGray)



#' Check local CNLOH

v <- finalSnvGray[[34]]
b <- finalBBGray[[34]]
b$major_cn[b$total_cn==2] <- 2
b$minor_cn[b$total_cn==2] <- 0
b$timing_param <- NULL

L <- computeMutCn(v, b, clusters=finalClustersGray[[34]], purity=finalPurityGray[34], xmin=3, gender='male', isWgd=FALSE, n.boot=0)
info(v)[colnames(L$D)] <- L$D
b$timing_param <- L$P
cls <- classifyMutations(v, reclassify='all')

info(v)$CLS <- cls
plotVcf(vcf = v, bb = b, clusters = finalClustersGray[[34]], col = col, ID = names(finalSnvGray)[34],  IS_WGD = FALSE, NO_CLUSTER = FALSE, title=TRUE)



library("Biostrings")
seq = "ATGTCCATGCTCGTGGTCTTTCTCTTGCTGTGGGGTGTCACCTGGGGCCCAGTGACAGAAGCAGCCATATTTTATGAGACGCAGCCCAGCCTGTGGGCAGAGTCC"
seq1 = strsplit(seq,split="")[[1]] # This was my original way of storing the sequences
seq2 = DNAString(seq) # This is the new way of storing the sequences

system.time( for (j in 1:10000) { seq1[10] } ) # 0.006 seconds
system.time( for (j in 1:10000) { seq2[10] } ) # 8.9 seconds
system.time( for (j in 1:10000) { subseq(seq2,10,10) } ) # 4.8 

system.time( for (j in 1:10000) { seq2[10] } ) 


sapply(split(t, donor2type[sample2donor[names(finalSnv)]]), function(x) {if(sum(!is.na(x))>1 & length(unique(x)) > 2) quantile(jitter(x), seq(0,1,0.1), na.rm=TRUE) else 1:10})


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

aggregatePerChromosome2 <- function(bb, isWgd=FALSE){
	s <- split(bb, seqnames(bb))
	r <- t(sapply(s, fractionGenomeWgdCompatible))
	r <- r[c(1:22,"X"),]
	w <- fractionGenomeWgdCompatible(bb)
	r <- rbind(r,WGD=w)
	return(r)
}

allChrAgg <- simplify2array(mclapply(finalBB, aggregatePerChromosome, mc.cores=2))


t <- allChrAgg[1:23,"time",!isWgd]
t[allChrAgg[1:23,"w",!isWgd] < diff(chrOffset)[1:23]*.33] <- NA

s <- split(as.data.frame(t(t)), droplevels(donor2type[sample2donor[names(finalSnv)]])[!isWgd])
n <- 10
#at <- function(x, n){
#	i <- sum(!is.na(x))
#	if(i > 16)
#		m <- n
#	else if(i> 8)
#		m <- n/2
#	else if(i > 4)
#		m <- n/4
#	else
#		m <- n/8
#	t <- table(cut(x, breaks=seq(0,1,1/m), include.lowest=TRUE))
#	t[rep(1:m, each=n/m)]*m/n
#}

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

#col <- colorRampPalette(RColorBrewer::brewer.pal(11, "PRGn")[-c(1,11)])(n)
#col <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired")[c(10,1,4)])(n)
prgn <- RColorBrewer::brewer.pal(11,"PRGn")
set1 <- RColorBrewer::brewer.pal(9,"Set1")

col <- colorRampPalette(set1[c(4,9,3)])(n)

p <- 0
v <- table(droplevels(donor2type[sample2donor[names(finalSnv)]]))
h <- (allChrCancerHist + p)  / rep(v + p, each=prod(dim(allChrCancerHist)[1:2]))
#h[rep(asum(allChrCancerHist,3) < 10, each=n)] <- 0
h <- aperm(h, c(2,3,1))

#ma <- function(x) (c(x[-1], x[length(x)]) + c(x[1], x[-length(x)])+2*x)/4 

a <- colMeans(h[c("All","WGD"),,] * c(23/5,1)) %*% 1:n / asum(h* c(23/5,1), c(1,3))
#a <- h["All",,] %*% 1:n / asum(h["All",,],2) /5 * 23
o <- order(-a)
h <- h[,o,]
w <- v[o]>=20 & apply(h, 2, max) > 0.05*8/n
h <- h[,w,]

#h <- aperm(apply(h, c(1,3), ma), c(2,1,3))

#scl <- 0.75
#l <- cbind(1:nrow(h),rep(1:dim(h)[2], each=nrow(h)))
#mg14:::stars(matrix(sqrt(h)*scl, ncol=n)[,n:1], scale=FALSE, draw.segments=TRUE, col.segments=c(col,NA), init.angle=acos(0), locations=l, plot=TRUE)
#text(x=0, paste0(colnames(h), ", n=", v[o][w]), y=seq_along(colnames(h)), las=2, pos=2)
#text(y=dim(h)[2]+1, 1:nrow(h), rownames(h), pos=3)
#symbols(x=rep(nrow(h)+1, ncol(h)), y=1:ncol(h), circles=sqrt(10/v[o][w])*scl, inches=FALSE, add=TRUE, fg="gray")
##mg14:::stars(matrix(sqrt(h), ncol=n)[,n:1]*1.3, scale=FALSE, draw.segments=TRUE, col.segments=c(col,NA), init.angle=acos(0), locations=l, plot=TRUE, add=TRUE)
#
#dev.copy2pdf(file="clocks.pdf",width=10, height=8, pointsize=8)

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
dev.copy2pdf(file="hist.pdf",width=6, height=6, pointsize=8)


vv <- v[dimnames(h)[[2]]]
vv <- vv/sum(vv)

hh <- matrix(matrix(aperm(h, c(1,3,2)), ncol=length(vv)) %*% vv, nrow=nrow(h))
rownames(hh) <- rownames(h)

par(mfrow=c(1,2), mar=c(3,3,2,2)+.6, mgp=c(2.5,.5,0), xpd=NA)
barplot(hh["WGD",], space=0, col=rev(col), xlab="Time [mutations]", ylab="Relative frequency", width=0.1, ylim=c(0,.065), yaxs='r')
axis(side=1)
barplot(hh["All",], space=0, col=rev(col), xlab="Time [mutations]", ylab="Relative frequency", width=0.1, ylim=c(0,.065), yaxs='r')
axis(side=1)
dev.copy2pdf(file="hist-pancan.pdf",width=8, height=4, pointsize=8)


# Medullo
s <- c(1:22, "X","Y")
l <- as.numeric(width(refLengths[seqnames(refLengths) %in% s]))
names(l) <- s
c <- cumsum(l)
plot(NA,NA, ylab="Sample",xlab="",xlim=c(0,sum(l)), ylim=c(0,146), xaxt="n")
axis(side=1, at=c(0,c), labels=rep('', length(l)+1))
mtext(side=1, line=1, at=chrOffset[1:24] + diff(chrOffset[1:25])/2, text=names(chrOffset[1:24]))
i <- 1
for(bb in finalBB[donor2type[sample2donor[names(finalSnv)]]=="CNS-Medullo"]){
	bb <- bb[which(bb$level %in% c("a","b","c","d") & bb$time.star=="***" & bb$n.snv_mnv > 19)]
	try({
				s <- start(bb) + chrOffset[as.character(seqnames(bb))]
				e <- end(bb) + chrOffset[as.character(seqnames(bb))]
				c <- cut(bb$time, seq(0,1,0.1), include.lowest=TRUE)
				segments(s,i,e,i, col=rev(RColorBrewer::brewer.pal(10, "PRGn"))[c], lwd=4)
			})
	i <- i+1
}

a <- allChrAgg[,"time",]
a["WGD",!isWgd] <- NA
a[1:23,][allChrAgg[1:23,"w",] < diff(chrOffset)[1:23]*.33] <- NA

for(t in c("CNS-GBM", "CNS-Medullo")){
w <- which(donor2type[sample2donor[names(finalBB)]]==t)
a0 <- a[1:23,w]
a0[is.na(a0)] <- -2
o <- rev(c(w[isWgd[w]][hclust(dist(t(a0[,isWgd[w]])))$order], w[!isWgd[w]][hclust(dist(t(a0[,!isWgd[w]])))$order]))
image(z=a[,o], x=1:24, col=rev(col), main=t, useRaster=TRUE)
}


v <- donor2type[sample2donor[names(finalSnv)]] %in% c("CNS-GBM","CNS-Medullo")
allChrAgg2 <- simplify2array(mclapply(finalBB[v], aggregatePerChromosome2, mc.cores=2))

a <- allChrAgg2[,"time.wgd",]
a[allChrAgg2[,"nt.wgd",]==0] <- NA
a[1:23,][allChrAgg2[1:23,"nt.wgd",] < diff(chrOffset)[1:23]*.33] <- NA
a["WGD",!isWgd[v]] <- NA

for(t in c("CNS-GBM", "CNS-Medullo")){
	w <- which(donor2type[sample2donor[names(finalBB)[v]]]==t)
	a0 <- a[1:23,w]
	a0[is.na(a0)] <- -2
	o <- rev(c(w[isWgd[v][w]][hclust(dist(t(a0[,isWgd[v][w]])))$order], w[!isWgd[v][w]][hclust(dist(t(a0[,!isWgd[v][w]])))$order]))
	image(z=a[,o], x=1:24, col=rev(col), main=t, useRaster=TRUE)
}


cancerTiming <- read.table("../final/timed_segments_AD_consensus_clustering_2muts_threshold_chrArms_31-05-2017.txt", header=TRUE, sep="\t")

apc <-
		function(bb, isWgd=FALSE){
	.aggregateSegments <- function(m){
		#m <- mcols(bb)
		t <- weighted.mean(m$pi0, m$N, na.rm=TRUE)
		n <- sum(m$N, na.rm=TRUE)
		sd <- sd(m$pi0, na.rm=TRUE)
		ci <- weighted.mean(m$uCI-m$lCI, m$N, na.rm=TRUE)
		w <- sum(m$end - m$end, na.rm=TRUE)
		c(time=t, n=n, sd=sd, ci=ci,w=w)
	}
#	if(!isWgd){
	s <- split(bb, bb$chromosome)
	r <- t(sapply(s, .aggregateSegments))
	r <- r[c(1:22,"X"),]
#	}else{
	w <- .aggregateSegments(as.data.frame(bb))
	r <- rbind(r,WGD=w)
#	}
	return(r)
}

aca <- simplify2array(mclapply(split(cancerTiming, cancerTiming$Sample), aggregatePerChromosome, mc.cores=2))

par(mfrow=c(2,2), xpd=TRUE)

plot(aca[-24,"time", ], allChrAgg[-24,"time", dimnames(aca)[[3]]], cex=sqrt(aca[-24,"n", ]/100), xlab="CancerTiming", ylab="MutationTime.R")
abline(0,1, col='red')

s <- aca[-24,"n", ] > 10

x <- rep(purityPloidy[dimnames(aca)[[3]], "purity"], each=23)
y <- as.numeric(aca[-24,"time", ] - allChrAgg[-24,"time", dimnames(aca)[[3]]])
plot(x, y, cex=sqrt(aca[-24,"n", ]/100),  ylab="CancerTiming - MutationTime.R", xlab="Purity")
lines(seq(0,1,0.01), predict(loess(y~x, data=data.frame(x=x, y=as.numeric(y))), newdata=data.frame(x=seq(0,1,0.01))), col='red')
abline(h=0, lty=3, col='red')

fClonal <- sapply(finalClusters, function(x) x$n_ssms[1]/sum(x$n_ssms))

z <- rep(fClonal[dimnames(aca)[[3]]], each=23)
plot(z, y, cex=sqrt(aca[-24,"n", ]/100),  ylab="CancerTiming - MutationTime.R", xlab="Clonal fraction")
lines(seq(0,1,0.01), predict(loess(y~z, data=data.frame(z=z, y=y), subset=s), newdata=data.frame(z=seq(0,1,0.01))), col='red')
abline(h=0, lty=3, col='red')

summary(lm(bias ~ purity+f_clonal, data=data.frame(bias=y, purity=x, f_clonal=z)))

m <- as.numeric(aca[-24,"time", ])[s]
n <- as.numeric(allChrAgg[-24,"time", dimnames(aca)[[3]]])[s]
cor(m, n, use='c')

var(m)
var(n, na.rm=TRUE)


i <- paste(cancerTiming$Sample, cancerTiming$chromosome, cancerTiming$start, cancerTiming$end, sep=":")
j <- paste(timing_output$sample, timing_output$seqnames, timing_output$start, timing_output$end, sep=":")

m <- match(i,j)
#plot(timing_output$time[m], cancerTiming$pi0, cex=sqrt(cancerTiming$N/100), col=timing_output$time.star[m], pch=as.numeric(cancerTiming$type)); abline(0,1)
for(w in levels(timing_output$time.star)){
plot(cancerTiming$pi0 ~ timing_output$time[m], subset=timing_output$time.star[m]==w, cex=sqrt(cancerTiming$N/100), col=timing_output$time.star[m],  pch=as.numeric(cancerTiming$type), main=w, xlab="MutationTime.R", ylab="CancerTiming"); abline(0,1)
}

#plot(timing_output$n.snv_mnv[m], cancerTiming$N, log='xy')

d <- timing_output$time[m] - cancerTiming$pi0
boxplot(d ~ cancerTiming$type + timing_output$time.star[m] , las=2, ylab="CancerTiming - MutationTime.R")
abline(h=0, col='red')

summary(lm(d ~ cancerTiming$type*timing_output$time.star[m], weights=cancerTiming$N/100 ))

w <- timing_output$time.star[m] != "*" & cancerTiming$Sample %in% names(which(isWgd))

v <- sapply(split(cancerTiming$pi0[w], cancerTiming$Sample[w]), sd, na.rm=TRUE)
v2 <- sapply(split(timing_output$time[m][w], cancerTiming$Sample[w]), sd, na.rm=TRUE)

qqplot(v,v2, log='xy', xlim=c(0.01,1), ylim=c(0.01,1)); abline(0,1)




#' Cruch all data
allSubsStrand <- simplify2array(mclapply(finalSnv, tableSubsStrandRep, genes=genes, repStrand=repStrand), mc.cores=2)

#' Reorder subs to match rev comp
l <- levels(tncToPyrimidine(vcf))
n <- DNAStringSet(gsub("\\[|\\]|>","",l))
r <- paste0(as.character(complement(DNAStringSet(substr(n,4,4)))), "[", as.character(complement(DNAStringSet(substr(n,2,2)))),">",as.character(complement(DNAStringSet(substr(n,3,3)))),"]",as.character(complement(DNAStringSet(substr(n,1,1)))))
m <- match(c(l,r), dimnames(allSubsStrand)[[1]])

allSubsStrand <- allSubsStrand[m,,,]

#' Condense into matrix
a <- matrix(allSubsStrand, ncol=dim(allSubsStrand)[4])
colnames(a) <- dimnames(allSubsStrand)[[4]]
write.table(data.frame(sub=rep(dimnames(allSubsStrand)[[1]], dim(allSubsStrand)[2]), strand=rep(dimnames(allSubsStrand)[[2]],each=dim(allSubsStrand)[1]),a, check.names=FALSE), file="allSubsStrand.txt", col.names=TRUE,sep="\t", quote=FALSE)
write.table(data.frame(sample=names(finalSnv), type=donor2type[sample2donor[names(finalSnv)]]), col.names=TRUE, sep="\t", quote=FALSE, file="allTumourTypes.txt")

#write.table(as.data.frame.table(allSubsStrand), col.names=TRUE,sep="\t", quote=FALSE,file="allSubsStrand.flat.txt")
d <- as.data.frame.table(allSubsStrand)

d <- data.frame(cbind(sub=dimnames(allSubsStrand)[[1]], ts=rep(dimnames(allSubsStrand)[[2]],each=dim(allSubsStrand)[1]), rs=rep(dimnames(allSubsStrand)[[3]],each=prod(dim(allSubsStrand)[1:2]))),a, check.names=FALSE)
write.table(d, file="allSubsStrand2.txt", col.names=TRUE,sep="\t", quote=FALSE, row.names=FALSE)


#' Per chromosome table
tableSubsChr <- function(vcf){
	s <- getTrinucleotideSubs(vcf = vcf)
	c <- factor(as.character(seqnames(vcf)), levels=c(1:22,"X","Y"))
	table(sub=s, chr=c)
}

allSubsChr <-  simplify2array(mclapply(finalSnv, tableSubsChr, mc.cores=2))
allSubsChr <- allSubsChr[m,,]
a <- matrix(allSubsChr, ncol=dim(allSubsChr)[3])
write.table(data.frame(sub=rep(dimnames(allSubsChr)[[1]], dim(allSubsChr)[2]), chr=rep(dimnames(allSubsChr)[[2]],each=dim(allSubsChr)[1]),a, check.names=FALSE), file="allSubsChr.txt", col.names=TRUE,sep="\t", quote=FALSE)


foo <- allSubsChr[,1:22,which(abs(finalPloidy-2) < 0.025)[1:100]]
foo <- allSubsChr[,1:22,1:100]

s2t <- gsub("\\[|\\]|>.","",rownames(foo))

s <- asum(foo, 2)


r <- BSgenome.Hsapiens.1000genomes.hs37d5#readDNAStringSet(refFile)
tfc <- sapply(seqnames(r)[1:22], function(x) trinucleotideFrequency(r[[x]]))
tf <- asum(tfc,2)

rtfc <- tfc/tf

mu <- sapply(1:dim(s)[2], function(i) s[,i]* rtfc[s2t,], simplify='array')


v <- (foo-mu)^2

w <- which(asum(foo,1:2)>1000)
fit <- glm.nb(as.numeric(foo[,,1:10]) ~ as.numeric(log(mu[,,1:10]+1e-3)))
summary(fit)

cv <- v/mu^2

par(mfrow=c(1,2))

plot(mu[,,w], foo[,,w], log='xy', xlab="Prediction", ylab="Observation")
abline(0,1,col="red")

plot(mu[,,w], cv[,,w], log='xy', xlab="Prediction", ylab="CV")
x <- 10^seq(-2,4,0.01)
lines(x, 1/x, col=2, lty=3)
lines(x, 1/x + 1/fit$theta, col='red')


fractionDiploid <- function(bb){
	w <- which(bb$major_cn==1 & bb$minor_cn==1 & countSubjectHits(findOverlaps(bb,bb)) == 1)
	sum(as.numeric(width(bb)[w]), na.rm=TRUE)/sum(as.numeric(width(bb)))
}

fractionDiploidChr <- function(bb){
	sapply(seqlevels(bb), function(s) fractionDiploid(bb[seqnames(bb)==s]))
}

diploidChrSample <- simplify2array(mclapply(finalBB, fractionDiploidChr, mc.cores=2))
write.table(file="diploidChrSample.txt", t(diploidChrSample), sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)


mu <- sapply((1:dim(allSubsChr)[3]), function(i) {
			w <- diploidChrSample[,i] > 0.99 & rownames(diploidChrSample) %in% 1:22
			si <- rowSums(allSubsChr[,w,i, drop=FALSE])
			tfi <- tfc
			tfi[,!w[1:22]] <- NA
			tf <- mg14:::asum(tfi,2, na.rm=TRUE)
			rtfc <- tfi/tf
			r <-  si* rtfc[s2t,]
			return(r)
		}, simplify='array')


v <- (allSubsChr[,1:22,]-mu)^2

w <- which(asum(allSubsChr[,1:22,],1:2)>50000)[1:100]
#w <- 1:100
chr <- rep(factor(rep(1:22, each=nrow(allSubsChr))), length(w))
sub <- rep(factor(1:nrow(allSubsChr)), 22*length(w))
fit <- glm.nb(as.numeric(allSubsChr[,1:22,w]) ~ as.numeric(log(mu[,,w]+1e-3)) + chr)
summary(fit)
fit$theta

p <- exp(predict(fit))

plot(p, fit$y, log='xy', xlim=c(1e-1, 1e05), xlab="Fit", ylab="Observed")

plot(p, (fit$y - p)^2/p^2, log='xy', xlim=c(1e-1, 1e05), xlab="Fit", ylab="CV")
x <- 10^seq(-2,4,0.01)
lines(x, 1/x, col=2, lty=3)
lines(x, 1/x + 1/fit$theta, col='red')


r <- as.numeric(strsplit("0.3198032
0.3431734
0.3640836
0.3690037
0.3542435
0.3456335
0.3542435
0.3468635
0.3739237
0.3603936
0.3542435
0.3530135
0.3567036
0.3431734
0.3419434
0.3726937
0.3493235
0.3480935
0.3763838
0.3493235
0.3554736
0.3665437
0.3554736
0.3788438
0.3505535
0.3862239
0.3591636
0.3480935
0.3517835
0.3665437
0.3493235
0.3739237
0.3591636
0.3849938
0.3677737
0.3407134
0.3653137
0.3628536
0.3567036
0.3677737
0.3468635
0.3542435
0.3665437
0.3788438
0.3456335
0.3431734
0.3800738
0.3677737
0.3800738
0.3567036", "\n")[[1]])

o <- 0.576
n <- 345+469

hist(r, 10, probability=TRUE, xlim=c(0.25, 0.6))
x <- seq(0.25,0.6,0.001)*n; plot(x,dnorm(x, mean(r)*n, sd(r)*n), type='l')
rug(r*n)
rug(o*n, col='red')
points(o, 0, pch="*")

isDel <- function(bb, isWgd=FALSE) {
	t <- 1+isWgd
	(bb$major_cn < t) + (bb$minor_cn < t)
}

averageDel <- function(bb, isWgd=FALSE){
	sum(width(bb) * isDel(bb, isWgd=isWgd) * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
}

ad <- mapply(averageDel, finalBB, isWgd)

save(finalBB, clinicalData, sample2donor, donor2type, file='PCAWG-lara.RData')

par(mfrow=c(4,4), mar=c(3,3,1,1), mgp=c(2,0.5,0), cex=.5, bty="L")

v <- which(wgdStar=="likely" & !isWgd & fracGenomeWgdComp[,"time.wgd"] <.5)

for(j in 2){
	bb <- finalBB[[v[j]]]
	plotTiming(bb)
	
	w <- which(bb$major_cn > 2 & bb$minor_cn < 2 & bb$clonal_frequency == max(bb$clonal_frequency) & !is.na(bb$time))

	if(length(w) >0 ){
	a <- sapply(bb$timing_param[w], function(x){
				x[x[,"state"]==1, c("T.m.sX","T.m.sX.up","T.m.sX.lo"), drop=FALSE]
			}, simplify=FALSE)
	
	col=rep("black",3)#paste0(RColorBrewer::brewer.pal(5,"Set2")[c(3:5)],"88")
	for(i in seq_along(w)){
		s <- start(bb)[w[i]]
		e <- end(bb)[w[i]]
		x <- chrOffset[as.character(seqnames(bb))[w[i]]]
		time<- a[[i]][-c(1,nrow( a[[i]])),,drop=FALSE]
		t <- rev(cumsum(rev(a[[i]][-1,1]))[-1])
		time <- pmin(time+t,1.1)
		y <- time[,1]
		rect(s+x,time[,3],e+x,time[,2], border=NA, col=col[bb$type[w[i]]], angle = ifelse(bb$time.star[[w[i]]]=="*" | is.na(bb$time.star[[w[i]]]),45,135), density=ifelse(bb$time.star[w[i]] == "***", -1, 72))
		segments(s+x,y,e+x,y)
	}
} else next
}

.ap <- function(bb) {
	u <- !duplicated(bb)
	c <- if(!is.null(bb$copy_number)) bb$copy_number else bb$total_cn
	sum(as.numeric(width(bb)[u] * c[u]) , na.rm=TRUE) / sum(as.numeric(width(bb)[u]), na.rm=TRUE)
}

.ah <- function(bb){
	u <- !duplicated(bb)
	sum(as.numeric(width(bb)[u] * (bb$minor_cn[u] == 0)) , na.rm=TRUE) / sum(as.numeric(width(bb)[u])  , na.rm=TRUE)
}

ap <- sapply(finalBB, .ap)
ah <- sapply(finalBB, .ah)

w <- .classWgd(ap, ah)
plot(ah, finalHom)


write.table(data.frame(sample=names(finalBB), avgHom=finalHom, avgPloidy=finalPloidy, WGD=isWgd),"WGD-final.txt", col.names=TRUE, sep="\t", quote=FALSE, row.names=FALSE)

#' # Gains +2

load("two_gain_times.RData")
doubleGains <- as.data.frame(T.i.F)
t <- aggregate()
m <- paste(doubleGains$sample, doubleGains$cnMaj, doubleGains$cnMin, doubleGains$chromosome, sep="_")
s <- split(doubleGains[,c("sample","tumor_type","T1_raw","T2_raw","n_mutations")], m)
doubleGainsAggregated <- Reduce("rbind",sapply(s, function(x) {
			data.frame(sample=x$sample[1], tumor_type=x$tumor_type[1], T1_raw=weighted.mean(x$T1_raw, x$n_mutations),T2_raw=weighted.mean(x$T2_raw, x$n_mutations), n_mutations=sum(x$n_mutations))
		}, simplify=FALSE))

par(bty="L", mar=c(3,3,0.5,0.5)+.5, mgp=c(2,0.5,0), tcl=-0.1)
x <- doubleGainsAggregated$T1_raw/pmax(1, doubleGainsAggregated$T2_raw)
y <- doubleGainsAggregated$T2_raw/pmax(1, doubleGainsAggregated$T2_raw)
o <- order(doubleGainsAggregated$n_mutations, decreasing=TRUE)
plot(x[o], 
		y[o], 
		col=tissueColors[as.character(donor2type[sample2donor[names(finalSnv)[doubleGainsAggregated$sample[o]]]])], pch=19,
		xlab="Time [mutations], first gain",
		ylab="Time [mutations], second gain",
		cex=sqrt(doubleGainsAggregated$n_mutations[o]/500))
t <- table(doubleGainsAggregated$sample)

par(mfrow=c(5,5))
for(i in as.numeric(names(t)[t>5])[1:25]){
	w <- which(doubleGainsAggregated$sample==i)
	plot(x[w],y[w], col=tissueColors[as.character(donor2type[sample2donor[names(finalSnv)[doubleGainsAggregated$sample[w]]]])], type='p', xlim=c(0,1), ylim=c(0,1), 
			xlab="T1",
			ylab="T2",
			pch=19,
			cex=sqrt(doubleGainsAggregated$n_mutations[w]/500))
	
}


#' Barplot of all drivers
g <- asum(finalGenotypes[,,,,selectedSamples[whiteList]], 2:4)
m <- model.matrix(~ droplevels(donor2type[sample2donor[colnames(g)]])-1)
colnames(m) <- levels(droplevels(donor2type[sample2donor[colnames(g)]]))
g <- g %*% m
g <- g[order(rowSums(g), decreasing=TRUE),]
n <- c("all",sapply(strsplit(rownames(g)[-1],"::"),`[`,3))
rownames(g) <- n
barplot(t(g[-1,])/2583, las=2, border=NA, col=tissueColors[colnames(g)], cex.names=.1, yaxt="n", ylab="Frequency")
mtext(side=2, at=c(0,0.1,0.2,0.3), text=c("0%","10%","20%","30%"), las=2, line=-2)
dev.copy2pdf(file="allDrivers.pdf", width=9, height=4, pointsize=8)



#' Mini TD?
vcf <- finalIndel[[1]]
vcf <- finalIndel[["00aa769d-622c-433e-8a8a-63fb5c41ea42"]]

up=scanFa(file="/lustre/scratch112/sanger/cgppipe/PanCancerReference/genome.fa.gz", shift(resize(granges(vcf), 100,fix="end"),-1))
down <- scanFa(file="/lustre/scratch112/sanger/cgppipe/PanCancerReference/genome.fa.gz", shift(resize(granges(vcf), 100,fix="start"),1))
b <- scanFa(file="/lustre/scratch112/sanger/cgppipe/PanCancerReference/genome.fa.gz", resize(granges(vcf), 1,fix="start"))

a <- unlist(alt(vcf))
r <- substring(a,1,1)
a <- substring(a,2, width(a))
w <- which(width(a)>=1 & width(ref(vcf)) ==1)

ins <- data.frame(
		up2=substring(up, width(up)-2*width(a)+2, width(up)-width(a)+1)[w], 
		up1=substring(up, width(up)-width(a)+1, width(up))[w], 
		insertion=as.character(a[w]), 
		down1=substring(down, 1,width(a))[w],
		down2=substring(down, width(a)+1,2*width(a))[w], 
		stringsAsFactors=FALSE)

table(width=width(ins$insertion),`duplication`=ins$ins == ins$up1 | ins$ins == ins$down1, `repeat`=ins$up1==ins$up2 | ins$down1==ins$down2)


g <- mg14:::asum(finalGenotypes[-555,,,,], 2:4)

u <- apply(g,2,paste, collapse=":")

t <- table(u)
w <- which(t <=27)
sum(t[w])


#' Hotspots

files <- dir("../bam/tumour", pattern=".bam$")
hotspot_counts <- sapply(seq_along(files), function(i){
			e <- new.env()
			load(paste0("../scratch/hotspots/",i,".RData"), envir=e)
			return(e$counts[1,,])
		}, simplify='array')

dimnames(hotspot_counts)[[3]] <- sub(".bam","", files)
dimnames(hotspot_counts)[[2]] <- c("A","T","C","G","-","a","t","c","g","_")

hotspot_counts <- hotspot_counts[,,names(finalSnv)]

save(hotspot_counts, file="hotspot_counts.RData")

c <- mg14:::asum(hotspot_counts,2)

tab <- read.table("Protected_sites_shearwater_normalpanel_after_filteringSNPs.txt", sep="\t", header=FALSE)
hotspots <- VRanges(tab$V1, IRanges(tab$V2, width=1), ref=tab$V3, alt=tab$V4)

load("hotspots_bf.RData")

counts <- aperm(hotspot_counts, c(3,1,2))

counts_h <- sapply(seq_along(tab$V4), function(i) counts[,i,c(as.character(tab$V4[i]), tolower(as.character(tab$V4[i])))], simplify='array')
counts_h <- aperm(counts_h, c(1,3,2))

bfh <- sapply(seq_along(tab$V4), function(i) bf[,i,as.character(tab$V4[i])])
bf_noth <- sapply(seq_along(tab$V4), function(i) bf[,i,!dimnames(bf)[[3]]%in%as.character(tab$V4[i])], simplify='array')


sum(bfh[selectedSamples,] < 0.05)/sum(selectedSamples)

h <- bfh[selectedSamples,] < 0.05

pfh <- 1/(1+1/bfh)

o <- sapply(1:nrow(h), function(i){
			sum(hotspots[which(h[i,])] %over% finalSnv[[which(selectedSamples)[i]]])
		})

w <- which(t(h), arr.ind=TRUE)

op <- sapply(1:nrow(q), function(i){
			sum(hotspots[which(q[i,]<0.1)] %over% finalSnv[[which(selectedSamples)[i]]])
		})


#' All mutations
allSubsTransRepClust <- rhdf5::h5read("allSubsTransRepClust3.h5","allSubsTransRepClust")
X <- t(matrix(allSubsTransRepClust, ncol=dim(allSubsTransRepClust)[[5]]))
a <- allSubsTransRepClust
dim(a) <- c(96,2,dim(a)[-1])
allTable <- as.matrix(read.table("allTable.txt", sep="\t", check.names=FALSE))
Z <- cbind(t(mg14:::asum(a, 2:5)), allTable[,-c(1:96, 331:475)])
w <- which(!is.na(allTable[,330]))
Z <- Z[w,]
any(is.na(Z))

cosineDist <- function(x,y){
	t(x)%*%y/(sqrt(colSums(x^2)  %*% t(colSums(y^2))) )
}

#lrtDist <- function(x,y){
#	z <- x+y 
#	p <- z/sum(z)
#	d <- 2*(dmultinom(x, prob=x/sum(x), log=TRUE) + dmultinom(y, prob=y/sum(y), log=TRUE) - dmultinom(x, prob=p, log=TRUE) - dmultinom(y, prob=p, log=TRUE))
#	#pv <- pchisq(d, df= length(p)-1, lower.tail=FALSE)
#	return(d)
#}
#d <- sapply(1:ncol(z), function(i) sapply(1:ncol(z), function(j) lrtDist(z[,i], z[,j])))



z <- t(cbind(X,allTable[,-c(1:96, 331:475)]))[,w]
r <- colSums(z)
d <- 1-cosineDist(z,z)
d0 <- d
d <- d0 + abs(log10(outer(r,r,"/")))/10
#d <- 1-d
 
set.seed(42)
tsne <- Rtsne::Rtsne(as.dist(d), theta=0.0, pca=FALSE)
t <- tsne$Y


o <- order(colSums(z), decreasing=TRUE)

pdf("tSNE-TRC-MNV-indel-SV.pdf",16,16)
scale <- 5000
b <- ifelse(colMeans(col2rgb(tissueColors)/255) > 0.25, "black","lightgrey")
plot(-t[o,1],t[o,2], bg=tissueColors[as.character(donor2type[sample2donor[colnames(z)[o]]])], pch=21, cex=sqrt(colSums(z)[o]/scale)+0.5, xlab="t-SNE 1", ylab="t-SNE 2", col=b[as.character(donor2type[sample2donor[colnames(z)[o]]])], 
		bty="n", xaxt="n", yaxt="n", xpd=NA)
legend("topleft", pch=21, pt.bg=tissueColors, legend=names(tissueColors), bty="n", ncol=3, pt.cex=2)

Z0 <- Z/rowSums(Z)

col330 <- c(rep(RColorBrewer::brewer.pal(7,"Set1")[-6], each=16), 
		rep("darkgrey",91), 
		rep(RColorBrewer::brewer.pal(8,"Dark2")[7:8], each=20),
		"black", rep("pink",102))

x <- round(range(t[,1])+c(-10,10),-1)
x <- seq(x[1],x[2],5)
y <- round(range(t[,2])+c(-10,10),-1)
y <- seq(y[1],y[2],2.5)
par(mfrow=c(length(y)-1, length(x)-1), mar=c(0.5,0.1,0.1,0.1), xpd=NA)
for(j in rev(seq_along(y)[-1])){
	for(i in rev(seq_along(x))[-1]){
		w <- which(t[,1]>x[i-1] & t[,1]< x[i] & t[,2]> y[j-1] & t[,2]< y[j])
		try({barplot((colMeans(Z0[w,,drop=FALSE])), border=NA, col=col330, ylim=c(0,0.25), xaxt="n", yaxt="n", las=2)->b; mtext(text=colnames(Z), at=b, las=2, line=0,side=1, cex=0.01, xpd=NA )})
	}
}

dev.off()


#' Expression data

fpkm <- read.table("../fpkm/joint_fpkm_uq.tsv.gz", header=TRUE, sep="\t", check.names=FALSE)
g <- fpkm$feature
fpkm <- as.matrix(fpkm[,-1])
rownames(fpkm) <- g

load("../ref/gencode.v19.annotation.gtf.RData")
gencode <- gencode[gencode$type=="gene"]

genes <-  strsplit("MLH1, MLH3, MSH2, MSH3, MSH6, PMS1, PMS2, POLD1, POLE, POLG, POLH, POLK, REV1, APOBEC3A, APOBEC3B, APOBEC3G, BRCA1, BRCA2, RIF1, TP53, MBD4, NEIL1, MGMT", ", ")[[1]]

ensg <- gencode$gene_id[match(genes, gencode$gene_name)]
names(ensg) <- genes

m <- match(ensg, rownames(fpkm))

g2t <- as.character(finalData$tumor_rna_seq_aliquot_id)
names(g2t) <- finalData$sanger_variant_calling_file_name_prefix
g2t[g2t==""] <- NA

scale <- 5000
b <- ifelse(colMeans(col2rgb(tissueColors)/255) > 0.25, "black","lightgrey")
c <- rev(paste(RColorBrewer::brewer.pal(9, "Spectral"),"AA", sep=""))
s <- colnames(z)[o]
f <- fpkm[m,match(g2t[s], colnames(fpkm))]
q <- t(apply(f, 1, function(x) {
#					y <- unlist(sapply(split(x, donor2type[sample2donor[s]]), 
#									function(xx){n <- xx/mean(xx, na.rm=TRUE)
#										#log(n+0.01)
#										rank(n, na.last="keep")/sum(!is.na(n))
#									}
#										))[match(s, unlist(split(s, donor2type[sample2donor[s]])))]
					rank(x,na.last="keep")/sum(!is.na(x))
					#as.numeric(cut(x,quantile(x,seq(0,1,0.01), na.rm=TRUE), include.lowest=TRUE))/100
}))
par(mfrow=c(5,5))
#pdf("tSNE-fpkm.pdf",16,16)
for(i in seq_along(ensg)){
 x <- q[i,]
 plot(-t[o,1],t[o,2], bg=c[cut(x,9)], pch=21, cex=sqrt(colSums(z)[o]/scale)+0.5, xlab="t-SNE 1", ylab="t-SNE 2",bty="n", xaxt="n", yaxt="n", xpd=NA, main=names(ensg)[i], col=NA)
}
#dev.off()
x <- colMeans(q)
#plot(-t[o,1],t[o,2], bg=c[cut(x,9)], pch=21, cex=sqrt(colSums(z)[o]/scale)+0.5, xlab="t-SNE 1", ylab="t-SNE 2",bty="n", xaxt="n", yaxt="n", xpd=NA, main="MMR")

w <- donor2type[sample2donor[s]]=="Skin-Melanoma"
f <- fpkm[,match(g2t[s], colnames(fpkm))]
cr <- cor(t(f[,w]), colSums(z[,s])[w], use="c", method='s')

medFpkm <- t(sapply(split(as.data.frame(t(f)), donor2type[sample2donor[s]]), sapply, median, na.rm=TRUE))
colnames(medFpkm) <- genes

mbd4 <- fpkm[gencode$gene_id[match("MBD4", gencode$gene_name)],]
foo <- rr
names(foo) <- NULL
u <- unlist(foo)

summary(lm(log(u) ~ log(mbd4[g2t[names(u)]]) ))

plot(mbd4[g2t[names(u)]], u, bg=tissueColors[donor2type[sample2donor[names(u)]]], col=tissueBorder[donor2type[sample2donor[names(u)]]], pch=21, log='xy')

plot(fpkm[gencode$gene_id[match("MLH1", gencode$gene_name)],match(g2t[rownames(Z)], colnames(fpkm))], Z[,"delT"], bg=tissueColors[donor2type[sample2donor[rownames(Z)]]], col=tissueBorder[donor2type[sample2donor[rownames(Z)]]], pch=21, log='xy')
plot(fpkm[gencode$gene_id[match("MSH2", gencode$gene_name)],match(g2t[rownames(Z)], colnames(fpkm))], Z[,"delT"], bg=tissueColors[donor2type[sample2donor[rownames(Z)]]], col=tissueBorder[donor2type[sample2donor[rownames(Z)]]], pch=21, log='xy')


KLD <- function(x,y){
	sum(log(x^x) - log(x^y) + log(y^y) - log(y^x))
}

S <- as.matrix(read.table("../ref/PCAWG_signature_patterns.txt", header=TRUE, sep="\t")[,-(1:2)])

f <- function(x, e=0.001) (x+e)/sum(x+e)

k <- sapply(1:ncol(S), function(i) sapply(1:ncol(S), function(j) KLD(f(S[,i]), f(S[,j])))) 
c <- cosineDist(S,S)


t <- do.call("rbind", lapply(timeWgd, `[`, , , "5x"))
library(survival)
s <- Surv(clinicalData$donor_survival_time, grepl("deceased",clinicalData$donor_vital_status))

ss <- s[match(sample2donor[rownames(t)], clinicalData$icgc_donor_id)]
l <- droplevels(donor2type[sample2donor[rownames(t)]])
m <- model.matrix(~l - 1)

coxph(ss~ t[,"hat"] +  age[sample2donor[rownames(t)]] + ridge(m))

for(ll in levels(l)){
	print(ll)
	try(print(summary(coxph(ss~ t[,"hat"] +  age[sample2donor[rownames(t)]], subset=l==ll))))
}


t <- do.call("rbind", lapply(timeSubclones, `[`, , , "5x"))
t <- t[order(t[,"hat"]),]
tt <- tail(na.omit(t), 100)

pdf("early-subclones.pdf", 4,4)
for(i in rownames(tt))
try({
	plotSample(i)
	title(xlab=i)
})
dev.off()

d <- finalDriversAnnotated[grep("DNMT3A",finalDriversAnnotated$ID)][5]
for(s in names(finalSnv)[donor2type[sample2donor[names(finalSnv)]] == "Myeloid-MPN"])
	print(nrow(findOverlaps(finalSnv[[s]], granges(d))))


n <- GRanges(5, IRanges(170814120,170837569 ))
for(s in names(finalSnv)[donor2type[sample2donor[names(finalSnv)]] == "Myeloid-AML"])
	print(nrow(findOverlaps(finalIndel[[s]], n)))

s <- names(finalSnv)[grep("Myeloid",donor2type[sample2donor[names(finalSnv)]] )]


s <- "f92a78d1-90ff-70c8-e040-11ac0d485eca"
bb <- finalBB[["f92a78d1-90ff-70c8-e040-11ac0d485eca"]]
bb[seqnames(bb)=="9"]$total_cn[1:2] <- 2 
bb[seqnames(bb)=="9"]$major_cn[1:2] <- 1 
bb[seqnames(bb)=="9"]$clonal_frequency[1:2] <- 0.940- 0.580
bb <- c(bb, bb[seqnames(bb)=="9"][1:2])
bb$total_cn[length(bb)+c(-1,0)] <- 2
bb$major_cn[length(bb)+c(-1,0)] <- 2
bb$minor_cn[length(bb)+c(-1,0)] <- 0
bb$clonal_frequency[length(bb)+c(-1,0)] <- 0.580
bb <- sort(bb)
bb$timing_param <- NULL


vcf <- finalSnv[["f92a78d1-90ff-70c8-e040-11ac0d485eca"]]
MCN <- computeMutCn(vcf, bb, finalClusters[[s]], finalPurity[[s]], n.boot=0)

bb$timing_param <- MCN$P

# JAK2
j <- GRanges(9, IRanges(5073770, width=1))
subsetByOverlaps(bb,j)

t <- subsetByOverlaps(bb,j)$timing_param[[1]]
t[1:2,"P.m.sX"] <- c(0.5,0.5)
posteriorMutCN(56,56+29, t)


# TET2
t <- subsetByOverlaps(bb,rowData(vcf)[grep("TET2",info(vcf)$DG)])$timing_param[[1]]
posteriorMutCN(info(vcf)$t_alt_count[grep("TET2",info(vcf)$DG)],
		info(vcf)$t_alt_count[grep("TET2",info(vcf)$DG)] + info(vcf)$t_ref_count[grep("TET2",info(vcf)$DG)] , t)

pdf(paste0(s, ".pdf"), 4,4)
plotSample(s, bb=bb)
dev.off()



d <- hotspots[!hotspots %over% finalDrivers]
t <- read.table("../ref/foo.tsv", header=TRUE, sep='\t')


permuteCn <- function(bb, cn.min=8){
	u <- unique(bb)
	w <- which(u$total_cn < cn.min & seqnames(u) %in% 1:22 & !sapply(u$timing_param, is.null))
	s <- sample(w)
	f <- findOverlaps(bb, u[w])
	bbb <- bb
	seqnames(bbb)[queryHits(f)] <- seqnames(u)[s][subjectHits(f)]
	bbb$n.snv_mnv[queryHits(f)] <- round((bbb$n.snv_mnv[queryHits(f)]+0.5) * as.numeric(width(ranges(u)[s][subjectHits(f)]))/as.numeric(width(bbb)[queryHits(f)]))
	ranges(bbb)[queryHits(f)] <- ranges(u)[s][subjectHits(f)]
	bbb <- sort(bbb)
	return(bbb)
}

simulateMutations <- function(bb, purity=max(bb$clonal_frequency, na.rm=TRUE),  n=40, rho=0.01, xmin=3){
	g <- (averagePloidy(bb)*purity + 2*(1-purity))
	V <- list(VRanges())#VRanges()
	for(i in which(!duplicated(bb)))
		if(bb$n.snv_mnv[i]>1 & !is.null( bb$timing_param[[i]]))try({
			cnStates <- bb$timing_param[[i]]
			p <- cnStates[,"pi.s"]* if(!any(is.na(cnStates[,"P.m.sX"]))) cnStates[,"P.m.sX"] else cnStates[,"pi.m.s"]
			pwr <- cnStates[,"power.m.s"]#(cnStates[,"power.s"] * cnStates[,"power.m.s"])
			s <- sample(1:nrow(cnStates), size=pmax(1,ceiling(bb$n.snv_mnv[i] * (p %*% (1/pwr)))), prob=p, replace=TRUE)
			f <- cnStates[s,"f"]
			mu.c <- (bb$total_cn[i]*purity + 2*(1-purity))/g * n
			c <- rnbinom(length(f), size=1/rho, mu=mu.c)
			x <- rbetabinom(n=length(f), size=c, prob=f, rho=rho)
			pos <- round(runif(length(f), min=start(bb)[i], max=end(bb)[i]))
			w <- which(x>=xmin)
			V[[i]] <- VRanges(seqnames=seqnames(bb)[i], IRanges(pos, width=1), ref="N", alt="A", totalDepth=c, altDepth=x)[w]
		})
	V <- do.call("c", V[!sapply(V, is.null)])
	sampleNames(V) <- "SAMPLE"
	v <- as(V, "VCF")
	exptData(v)$header@header$INFO <- rbind(header(v)@header$INFO,info(header(finalSnv[[1]]))[c("t_ref_count","t_alt_count"),])
	info(v)$t_alt_count <- altDepth(V)
	info(v)$t_ref_count <- totalDepth(V) - altDepth(V)
	return(v)
}

w <- which(isWgd)
simSnv <- list()
simBB <- list()
for(ww in w[which(!names(finalBB[isWgd])[1:200] %in% names(simBB))])try({
	print(ww)
	names(ww) <- names(finalSnv)[ww]
	set.seed(42)
	b <- finalBB[[ww]]# permuteCn(finalBB[[ww]])
	v <- simulateMutations(b, n=mean(getTumorDepth(finalSnv[[ww]]), na.rm=TRUE))
	b$timing_param <- NULL
	L <- computeMutCn(v, b, clusters=finalClusters[[ww]], purity=finalPurity[ww], gender=allGender[names(ww), "pred_gender"], isWgd=isWgd[ww], n.boot=0, xmin=0)
	
	i = header(v)@header$INFO
	exptData(v)$header@header$INFO <- rbind(i,mcnHeader())
	info(v) <- cbind(info(v), L$D)
	b$timing_param <- L$P 
	
#' Classify mutations
	cls <- classifyMutations(info(v), reclassify='all')
	info(v)$CLS <- cls
	info(header(v)) <- rbind(info(header(v)), DataFrame(Number="1",Type="String",Description="Mutation classification: {clonal [early/late/NA], subclonal}", row.names="CLS"))
	
#' Timing
	b$n.snv_mnv <- countOverlaps(b, v)
	t <- bbToTime(b)	
	mcols(b)[names(t)] <- DataFrame(t)
	
	simBB[[names(ww)]] <- b
	simSnv[[names(ww)]] <- v
#	plotSample(ww, vcf=v, bb=b)
})

plot(finalBB[[ww]]$time, b$time, cex=sqrt(finalBB[[ww]]$n.snv_mnv/100))

simBB <- simBB[order(names(simBB))]
simSnv <- simSnv[order(names(simSnv))]

save(simBB, simSnv, file="sim200.RData")

for(n in names(simBB)){
	writeVcf(simSnv[[n]], file=paste0("../sim200/vcf/",n,".vcf"))
	write.table(as.data.frame(simBB[[n]])[,c(1:9)], file=paste0("../sim200/cn/",n,".txt"), col.names=TRUE, sep="\t", quote=FALSE)
}

files <- dir("../sim200/annotated/cn", pattern=".RData", full.names=TRUE)
e <- new.env()
for(f in files){
	n <- sub("\\..+","",gsub(".+/","",f))
	load(f, envir=e)
	simBB[[n]] <- e$bb
}

simTime <- do.call("rbind",lapply(simBB, mcols))
realTime <- do.call("rbind", lapply(finalBB[names(simBB)], mcols))

msq <- sqrt((realTime$time - simTime$time)^2)

boxplot(msq ~ cut(simTime$n.snv_mnv, c(0,10^(0:5))), log='y', ylim=c(1e-5, 1))

w <- simTime$time >0.01 & simTime$time <0.99

par(mfrow=c(1,2))
plot(realTime$time, simTime$time, pch=16, cex=sqrt(realTime$n.snv_mnv/500), xlab="Time [simulated]", ylab="Time [estimated]", col="#00000044")
abline(0,1, col='red')

plot(simTime$n.snv_mnv[w],msq[w] + 1e-4, log='xy', pch=16, cex=1, xlab="Number of SNVs", ylab="Error",col="#00000044")
x <- c(10^seq(0,4.5,0.5))
q <- sapply(split( msq[w] +1e-6, cut(simTime$n.snv_mnv[w], x)), quantile, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=TRUE)
lines(x[-length(x)]*2.5, q["50%",], col='red', lwd=2)
lines(x[-length(x)]*2.5, q["25%",], col='red')
lines(x[-length(x)]*2.5, q["75%",], col='red')
lines(x[-length(x)]*2.5, q["2.5%",], col='red', lty=3)
lines(x[-length(x)]*2.5, q["97.5%",], col='red', lty=3)


plot(msq, realTime$time.up - realTime$time.lo)

table(simTime$time.up >= realTime$time & simTime$time.lo <= simTime$time, cut(simTime$n.snv_mnv, 10^(0:5)))

#' A few more checks re power and WCC
load("../final/annotated_010/clusters/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20160830.somatic.snv_mnv.complete_annotation.clusters_purity.RData")
bb <- finalBB[[1]]
bb$timing_param <- NULL
D <- computeMutCn(finalSnv[[1]], bb=bb, purity= purityPloidy[1,"purity"], clusters=clusters, isWgd=TRUE, gender='female', xmin=3, n.boot=0)
v <- finalSnv[[1]]
info(v)[colnames(D$D)] <- D$D

D0 <- computeMutCn(finalSnv[[1]], bb=bb, purity= purityPloidy[1,"purity"], clusters=wccClusters[[1]], isWgd=TRUE, gender='female', xmin=0, n.boot=0)
v0 <- finalSnv[[1]]
info(v0)[colnames(D0$D)] <- D0$D

table(classifyMutations(v), classifyMutations(v0))
table(info(v)$CLS, classifyMutations(v0))


b <- finalBB[[1]]
b$timing_param <- D$P
t <- bbToTime(b)

b0 <- finalBB[[1]]
b0$timing_param <- D0$P
t0 <- bbToTime(b0)

plot(t$time, t0$time)
plot(finalBB[[1]]$time, t$time)

cl <- clusters
cl$proportion <- wccClusters[[1]]$proportion

D3 <- computeMutCn(finalSnv[[1]], bb=bb, purity= purityPloidy[1,"purity"], clusters=cl, isWgd=TRUE, gender='female', xmin=3, n.boot=0)
v3 <- finalSnv[[1]]
info(v3)[colnames(D3$D)] <- D3$D
b0$timing_param <- D3$P
t3 <- bbToTime(b0)

cor(cbind(bb=finalBB[[1]]$time, t=t$time, t0=t0$time, t3=t3$time),use='c')

table(classifyMutations(v3), classifyMutations(v0))
table(classifyMutations(v3), info(v)$CLS)



ySubclone <- do.call("rbind", y)

m <- sapply(rownames(ySubclone),grep, as.character(TiN$tumor_wgs_aliquot_id))
plot(ySubclone[,"hat"], TiN[m, "TiN_donor"], col=tissueColors[donor2type[sample2donor[rownames(ySubclone)]]], pch=19, log='x', xlab="Time [yr]", ylab="TiN")
boxplot(ySubclone[,"hat"]~ TiN[m, "TiN_donor"], notch=TRUE)#, ylim=c(0,10))

t(sapply(split(TiN[m, "TiN_donor"], donor2type[sample2donor[rownames(ySubclone)]]), quantile, na.rm=TRUE))


sapply(split(ySubclone[,"hat"], TiN[m, "TiN_donor"]), quantile,na.rm=TRUE)
#m <- match(rownames(yy),as.character(TiN$tumor_wgs_aliquot_id))
t(sapply(split(TiN[m, "TiN_donor"], donor2type[sample2donor[rownames(ySubclone)]]), quantile, na.rm=TRUE))
sapply(split(TiN[m, "TiN_donor"], donor2type[sample2donor[rownames(ySubclone)]]), median, na.rm=TRUE)

yWgd <- do.call("rbind", y)

plot(yWgd[,"hat"], TiN[match(rownames(yWgd),as.character(TiN$tumor_wgs_aliquot_id)), "TiN_donor"], col=tissueColors[donor2type[sample2donor[rownames(yWgd)]]], pch=19, log='x', xlab="Time [yr]", ylab="TiN")

sapply(split(yWgd[,"hat"], TiN[match(rownames(yWgd),as.character(TiN$tumor_wgs_aliquot_id)), "TiN_donor"]), quantile,na.rm=TRUE)
boxplot(yWgd[,"hat"]~ TiN[match(rownames(yWgd),as.character(TiN$tumor_wgs_aliquot_id)), "TiN_donor"], notch=TRUE)

allClusters <- sapply(names(finalSnv), function(ID) consensusClustersToOld(loadConsensusClusters((ID))), simplify=FALSE)

se <- mg14:::asum(e$finalGenotypesP[,2,,], c(1))

c <- sapply(finalClusters, function(x) 1-x[1,"n_ssms"]/sum(x[,"n_ssms"]))

plot(s[4,]/colSums(s), c[whiteList])

plot(s[4,]/colSums(s), se[4,]/colSums(se))

plot(se[4,]/colSums(se), c[whiteList])


optim(c(1,1), function(x) sum((cbind(1,yy[w]) %*% x - xx[w])^2))

#' ## Signatures
#+ sigTable
sigTable <- simplify2array(mclapply(finalSnv, function(vcf) table(classifyMutations(vcf, reclassify="none"), tncToPyrimidine(vcf)), mc.cores=MC_CORES))
sigTable <- aperm(sigTable, c(2,3,1))

d <- nDeam22 *  t(sapply(finalWgdParam, function(x) if(!is.null(x$P)) x$P[[1]][1:2,"P.m.sX"] else c(NA,NA)))
rownames(d) <- names(finalSnv)[isWgd]
t0 <- 2*finalWgdPi["clonal.2","hat",]/( 2*finalWgdPi["clonal.2","hat",] +  finalWgdPi["clonal.1","hat",])
names(t0) <- dimnames(finalWgdPiAdj)[[5]]
d[rownames(d) %in% remove] <- NA
y <- d[names(t0),]/6
x <- t0
t <- donor2type[sample2donor[names(t0)]]
plot(x,y[,2], bg=tissueColors[t], pch=21,  col=tissueBorder[t], cex=tissueCex[t]*2*sqrt(y[,1]/100), lwd=0.5, xlab="Time", ylab="Early SNVs/Gb", log='')
p <- predict(loess(y~x, span=2), newdata=sort(x, na.last=NA), se=TRUE)
r <- function(x) c(x, rev(x))
polygon(r(sort(x, na.last=NA)), c(p$fit+2*p$se, rev(p$fit-2*p$se)), col="#00000044", border=NA)
lines(sort(x, na.last=NA),p$fit)

plot(x,y[,2], bg=tissueColors[t], pch=21,  col=tissueBorder[t], cex=tissueCex[t]*1, lwd=0.5, xlab="Time", ylab="Early SNVs/Gb", log='')
points(x,y[,1], bg=tissueColors[t], pch=21,  col=tissueBorder[t], cex=tissueCex[t]*1, lwd=0.5)


plot(x,y[,1]+y[,2], bg=tissueColors[t], pch=21,  col=tissueBorder[t], cex=tissueCex[t]*1, lwd=0.5, xlab="Time", ylab="Total SNVs/Gb", log='')
p <- predict(loess(rowSums(y)~x, span=2/3), newdata=sort(x, na.last=NA), se=TRUE)
r <- function(x) c(x, rev(x))
polygon(r(sort(x, na.last=NA)), c(p$fit+2*p$se, rev(p$fit-2*p$se)), col="#00000044", border=NA)
lines(sort(x, na.last=NA),p$fit)

#' Indels, ie MSI in remove?
t <- names(finalSnv)[which(TiN[sample2donor[names(finalSnv)]]> 0.01 | is.na(TiN[sample2donor[names(finalSnv)]]))]
h <- setdiff(remove, t)
i <- mg14:::asum(finalGenotypes[,"indels",,,], 1:3)

boxplot(i + 0.5 ~ I(names(i) %in% h), log='y')

#' Guesstimate acceleration
fmm <- 0.5* fm + 0.25*fmq[2,]+0.25*fmq[1,]
guessAccel <- (100/rbind(hat=fm, fmq) - 0.833)/0.167 #cut(fmm, c(0,0.25,0.5,0.75,1)*100, labels=paste0(rev(accel[-1]),"x"), include.lowest=TRUE)
#names(guessAccel) <- names(fm)
#table(guessAccel)

o <- order(fm)
barplot(pmin(200,guessAccel["hat",o]), col=tissueColors[colnames(guessAccel)][o], border=tissueLines[colnames(guessAccel)][o], las=2,names.arg=rep("",length(fm)) , ylab="Acceleration", log="y", ylim=c(1,200)) -> b
mg14::rotatedLabel(b, labels=names(fm[o]))
segments(b, pmin(200,guessAccel["hat",o]), b, guessAccel["2.5%",o], col=tissueLines[colnames(guessAccel)][o], lwd=2)
segments(b, guessAccel["97.5%",o], b, pmin(200,guessAccel["hat",o]), col=tissueBorder[colnames(guessAccel)][o], lwd=2)

setAccel <- cut(pmin(10,guessAccel["hat",]+1e-6), accel, include.lowest=TRUE, labels=paste0(accel[-1],"x"))
names(setAccel) <- colnames(guessAccel)

u <- setdiff(names(finalSnv)[uniqueSamples], remove)
qGuessWgd <- sapply(names(timeWgd), function(n){
			x <- timeWgd[[n]][,"hat",as.character(setAccel[n])]
			quantile(x[names(x) %in% u], na.rm=TRUE)
		})


y <- Reduce("c",deamRate)
y <- y[!names(y) %in% remove]
x <- age[sample2donor[names(y)]]
y <- y*x
t <- donor2type[sample2donor[names(y)]]
d <- data.frame(x,y,t)

f <- lme4::lmer(y ~ (x-1|t) + t-1, data=d)
r <- lme4::ranef(f)$t

fmr <- pmax(r[,2],0)*ma[rownames(r)]/(pmax(r[,2],0)*ma[rownames(r)] + pmax(0,r[,1]+lme4::fixef(f)))*100
