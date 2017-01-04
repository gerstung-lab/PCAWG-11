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
PCAWG_signatures <- read.table("../ref/PCAWG_signatures.tsv", header=TRUE, sep="\t")
PCAWG_S <- as.matrix(PCAWG_signatures[,-2:-1])
rownames(PCAWG_S) <- paste0(substr(PCAWG_signatures$Trinucelotide,1,1), "[", PCAWG_signatures$Mutation.Type, "]",substr(PCAWG_signatures$Trinucelotide,3,3))
PCAWG_S <- PCAWG_S[match(rownames(sigTable),rownames(PCAWG_S)),-35]
PCAWG_activity <- read.table("../ref/PCAWG_signatures_in_samples.tsv", header=TRUE, sep="\t", skip=1)
PCAWG_A <- as.matrix(PCAWG_activity[,-2:-1])
rownames(PCAWG_A) <- PCAWG_activity[,2]


PCAWG_sigDecomp <- array(0, dim=c(dim(sigTable)[-1], dim(PCAWG_S)[2]), dimnames=c(dimnames(sigTable)[-1], list(colnames(PCAWG_S))))

for(sample in rownames(PCAWG_sigDecomp)){
	if(class(try({
				w <- PCAWG_A[sample,]!=0
				E <- nmSolve(sigTable[,sample,], PCAWG_S[,w], maxIter=10000, tol=1e-4)
				PCAWG_sigDecomp[sample,,w] <- t1(E)
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
	
	correct.t <- function(pi, ta, a){
		pi <- pi[1:2]/sum(pi[1:2])
		t <- pi[1]/(pi[1]+2*pi[2])
		if(t < ta)
			crct <- c(1+(a-1)*(ta-t)/(1-t),1)
		else
			crct <- c(a,1+(a-1)*(t-ta)/t)
		return(crct)
	}
	
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

driVers <- VRanges(seqnames = drivers$chr, ranges=IRanges(drivers$pos, width =  width(r)), ref=r, alt=a, sampleNames  = drivers$sample_id)
mcols(driVers) <- drivers[,-c(1,3,4,5,6)]
driVers$MCN <- driVers$MNCN <- driVers$MJCN <- as.numeric(rep(NA,nrow(drivers)))
driVers$CLS <- factor(NA, levels = levels(info(allVcfNew[[1]])$CLS))

for(s in levels(sampleNames(driVers))){
	if(! s %in%  sampleIds) next
	for(x in c(allVcfNew[[s]], indelAnnotated[[s]])){
		w <- sampleNames(driVers) == s
		f <- findOverlaps(driVers[w] , x, select='first')
		mcols(driVers)[w[!is.na(f)],c("MCN","MNCN","MJCN","CLS")] <- info(x)[f[!is.na(f)], c("MCN","MNCN","MJCN","CLS")]
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