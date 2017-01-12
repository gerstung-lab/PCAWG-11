# TODO: Add comment
# 
# Author: mg14
###############################################################################

vcfPath <- '/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/final_consensus_12oct_passonly/snv_mnv'
basePath <-  '/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/dp/20161213_vanloo_wedge_consSNV_prelimConsCNAallStar'
dpPath <- paste0(basePath,'/2_subclones/')
cancerGenes <- read.table('/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/ref/cancer_genes.txt')$V1
purityPloidy <- read.table(paste0(basePath,'/1_purity_ploidy/purity_ploidy.txt'), header=TRUE, row.names=1)
colnames(purityPloidy) <- c("purity","ploidy")
cnPath <- paste0(basePath,'/4_copynumber/')
bbPath <- paste0(basePath,'/4_copynumber/')

gender <- read.table('/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/gender/2016_12_09_inferred_sex_all_samples.txt', header=TRUE, sep='\t')
gender <- gender[!duplicated(gender$tumourid) & gender$tumourid != 'tumourid',]
#write.table(gender,'/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/gender/2016_12_09_inferred_sex_all_samples_CORRECTED_MG.txt', row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
rownames(gender) <- gender$tumourid

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

# Density fucntions - truncated beta-binomial
#L <- matrix(sapply(pmin(cnStates[1:k,"f"],1), function(pp) dbetabinom(altCount[hh],tumDepth[hh],pp, 0.01) + .Machine$double.eps), ncol=k)
dtrbinom <- function(x, size, prob, xmin=0) dbinom(x,size,prob) / pbinom(xmin-1, size, prob, lower.tail=FALSE)
pbetabinom <- function(x, size, prob, rho){
	if(rho!=0)
		VGAM::pbetabinom(x, size, prob, rho)
	else
		pbinom(x, size, prob)
}
dbetabinom <- function(x, size, prob, rho){
	if(rho!=0)
		VGAM::dbetabinom(x, size, prob, rho)
	else
		dbinom(x, size, prob)
}
dtrbetabinom <- function(x, size, prob, rho, xmin=0) {y <- dbetabinom(x,size,prob,rho) / (1-pbetabinom(xmin-1, size, prob, rho))
	y[x<xmin] <- 0
	return(y)}
ptrbetabinom <- function(x, size, prob, rho, xmin=0) {
	pmin <- pbetabinom(xmin-1, size, prob, rho)
	pmax(0,(pbetabinom(x,size,prob,rho) - pmin) / (1-pmin))}


computeMutCn <- function(vcf, bb, clusters=allClusters[[meta(header(vcf))["ID",]]], purity =purityPloidy[meta(header(vcf))["ID",], 'purity'], gender='female', xmin=0, rho=0.001){
	n <- nrow(vcf)
	D <- DataFrame(MutCN=rep(NA,n), MutDeltaCN=rep(NA,n), MajCN=rep(NA,n), MinCN=rep(NA,n), MajDerCN=rep(NA,n), MinDerCN=rep(NA,n), CNF=rep(NA,n), CNID =as(vector("list", n),"List"), pMutCN=rep(NA,n), pGain=rep(NA,n),pSingle=rep(NA,n),pSub=rep(NA,n), pMutCNTail=rep(NA,n))	
	P <- vector(mode='list', length(bb))
	if(n==0)
		return(list(D=D, P=P))
	
	altCount <- getAltCount(vcf)
	tumDepth <- getTumorDepth(vcf)
	names(altCount) <- names(tumDepth) <- NULL
	ID <- meta(header(vcf))["ID",]
	
	# Match VCF and CN
	overlaps <- findOverlaps(vcf, bb)
	D[["CNID"]] <- as(overlaps, "List")
	majorCN <- split(bb$major_cn[subjectHits(overlaps)], queryHits(overlaps))
	m <- bb$minor_cn #hack: minor_cn > 0 in male samples - Battenberg bug?
	if(gender=='male')
		m[as.character(seqnames(bb)) %in% c('X','Y')] <- 0
	minorCN <- split(m[subjectHits(overlaps)], queryHits(overlaps))	
	h <- selectHits(overlaps, "first")
	H <- selectHits(overlaps, "last")
	
	cnNormal <- 2 - (gender=='male' & seqnames(vcf)=="X" | seqnames(vcf)=="Y")
	
	# Fix cluster and purity discrepancies
	clusters$proportion[which.max(clusters$proportion)] <- purity
	
	cloneFreq <- split(bb$clonal_frequency[subjectHits(overlaps)], queryHits(overlaps))
	cnStates <- matrix(0, nrow=10000, ncol=6)
	colnames(cnStates) <- c("state","m","f","n.m.s","pi.m.s","s")
	
	power.c <- rep(0, nrow(clusters))
	
	for(ii in 1:2)
	for( i in which(diff(c(-1, h)) !=0 | is.na(diff(c(-1, h)) !=0) )){
		if(!i %in% names(majorCN)) next
		if(is.na(h[i])) next
		if(if(i>1) h[i] != h[i-1] | is.na(h[i-1]) else TRUE){ #ie. new segment
			majcni <- majorCN[[as.character(i)]]
			mincni <- minorCN[[as.character(i)]]
			cfi <- cloneFreq[[as.character(i)]]
			effCnTumor <- sum((majcni + mincni)*cfi)
			effCnNormal <- as.numeric(cnNormal[i]) * (1-purity)
			
			if(any(is.na(majcni))) next
			
			mixFlag <- FALSE
			subclonalGainFlag <- FALSE
			clonalFlag <- TRUE
			majdelta <- 0
			mindelta <- 0
			
			majanc <- majder <- majcni
			minanc <- minder <- mincni
			
			if(length(cfi)>1){ # multiple (subclonal) CN states, if so add clonal option (ie. mixture of both states), subclonal states only change by 1..delta(CN)
				d <- colSums(abs(rbind(majcni, mincni) - c(1,1) * (1+ (purityPloidy[ID,2] > 2.7))))
				derived <- d == max(d) ## derived state further from diploid/tetraploid
				if(all(derived)) next 
				majanc <- majcni[!derived]
				minanc <- mincni[!derived]
				majder <- majcni[derived]
				minder <- mincni[derived]
				majdelta <- majcni[derived] - majcni[!derived]
				mindelta <- mincni[derived] - mincni[!derived]
				majcni <- c(min(majcni), # clonal, sub on allele that doesn't change
						(majcni[!derived] + cfi[derived]/purity*majdelta), # clonal, sub on allele that does change
						(majcni[derived] >0) + (majdelta > 0)) # subclonal, subs after or before CNA, m=1,delta
				mincni <- c(min(mincni), (mincni[!derived] + cfi[derived]/purity*mindelta), (mincni[derived] >0) + (mindelta > 0))
				majdelta <- c(0, cfi[derived]/purity*majdelta, majdelta)
				mindelta <- c(0, cfi[derived]/purity*mindelta, mindelta)
				cfi <- c(purity, purity,  cfi[derived])
				mixFlag <- c(FALSE, TRUE, FALSE)
				clonalFlag <- c(TRUE,TRUE,FALSE)
				subclonalGainFlag <- c(FALSE, FALSE, TRUE)
			}
			
			a <- sapply(clusters$proportion, function(p) all(abs(p-cfi) > 0.05)) # subclone(s) not coinciding with CN change
			if(any(a)){ # assume subclones have derived from most abundant CN state
				majcni <- c(majcni, rep(majcni[which.max(cfi)]>0, sum(a))+0)
				mincni <- c(mincni, rep(mincni[which.max(cfi)]>0, sum(a))+0)
				cfi <- c(cfi, clusters$proportion[a])
				mixFlag <- c(mixFlag, rep(FALSE, sum(a)))
				clonalFlag <- c(clonalFlag, rep(FALSE, sum(a)))
				subclonalGainFlag <- c(subclonalGainFlag, rep(FALSE, sum(a)))
				majdelta <- c(majdelta, rep(0,sum(a)))
				mindelta <- c(mindelta, rep(0,sum(a)))
			}
			totcni <- majcni+mincni
			if(all(totcni==0)) next
			
			k <- 0
			for( j in seq_along(majcni)){
				if(majcni[j]==0 & mincni[j]==0) {
					f <- m <- 0 # allele frequency
					pi.m.s <- n.m.s <- 1
					l <- 1
				}else{
					l <- 1:max(majcni[j], mincni[j]) # mincni>majcni can occur if minor allele changes subclonally
					m <- l
					n.m.s <- rep(1, length(l)) #multiplicity of cn states
					if(!mixFlag[j]){ # single subclone, or no subclonal cn
						f <- m * cfi[j] / (effCnTumor + effCnNormal)
						if(mincni[j] > 0)
							n.m.s[1:min(majcni[j], mincni[j])] <- 2
						pi.m.s <- n.m.s/sum(n.m.s)
					}else{ # coexisting cn subclones, use mixture
						d <- if(majdelta[j] != 0) majdelta[j] else mindelta[j] 
						M <- max(mincni[j]*(mindelta[j]!=0), majcni[j]*(majdelta[j]!=0))
						m <- (d > 0):M + (M-floor(M))
						l <- seq_along(m)
						f <- m *cfi[j] / (effCnTumor + effCnNormal) # Variant on major allele
						n.m.s <- rep(1, length(l))
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
			hh <- which(h==h[i] & !is.na(altCount) &! is.na(tumDepth))
			if(length(hh)==0) next
			whichStates <- (1:k)[cnStates[1:k,"f"]>0]
			
			# State probabilities - based on cancer cell fractions
			pi.s <- sapply(unique(cfi), function(p) ifelse(min(abs(clusters$proportion - p)) < 0.05, clusters$n_ssms[which.min(abs(clusters$proportion - p))], 1))
			pi.s <- pi.s/sum(pi.s)
			
			which.c <- sapply(unique(cfi), function(p) ifelse(min(abs(clusters$proportion - p)) < 0.05, which.min(abs(clusters$proportion - p)), NA)) # map to cluster
			
			cnStates[1:k,"s"] = as.numeric(factor(cfi, levels=unique(cfi)))[cnStates[1:k,"state"]]
						
			# Likelihood
			L <- matrix(sapply(pmin(cnStates[whichStates,"f"],1), function(pp) dtrbetabinom(altCount[hh],tumDepth[hh],pp, rho=rho, xmin=pmin(altCount[hh],0)) + .Machine$double.eps), ncol=length(whichStates))

			# Power
			power.sm <- colMeans(matrix(sapply(pmin(cnStates[whichStates,"f"],1), function(pp) 1-ptrbetabinom(pmin(altCount[hh],xmin),tumDepth[hh],pp, rho=rho, xmin=0)), ncol=length(whichStates)), na.rm=TRUE)
			if(ii==2){
				P.m.sX <- P[[h[i]]][,"P.m.sX"]
				power.s <- sapply(split(power.sm * P.m.sX, cnStates[whichStates,"s"]), sum) # Power for state
				power.s[!is.na(which.c)] <- power.c[which.c]
				power.m.s <- power.sm # Relative power of m states (local) conditioned on s (global).
		        for(s in unique(cnStates[whichStates,"s"])) power.m.s[cnStates[whichStates,"s"]==s] <- power.m.s[cnStates[whichStates,"s"]==s] / max(power.m.s[cnStates[whichStates,"s"]==s]) #power.s[s]
			}else{ # First iteration
				power.m.s <- power.sm
				power.s <- rep(1,length(whichStates))
			}
	
			# EM algorithm (mixture fitting) for pi
			P.m.sX <- cnStates[whichStates,"pi.m.s"]
			for(em.it in 1:100){
				P.xsm <- L * rep(pi.s[cnStates[whichStates,"s"]] * P.m.sX / power.m.s / power.s[cnStates[whichStates,"s"]], each=nrow(L)) # P(X,s,m)
				P.sm.x <- P.xsm/rowSums(P.xsm) # P(s,m|Xi)
				P.sm.X <- colMeans(P.sm.x) # P(s,m|X) / piState[cnStates[1:k,"state"]] / cnStates[1:k,"pi.m.s"]
				P.s.X <- sapply(split(P.sm.X, cnStates[whichStates,"s"]), sum)
				P.m.sX <- P.sm.X / P.s.X[cnStates[whichStates,"s"]]
			}
			
			if(ii==1){
				p <- (sapply(split(power.sm * P.m.sX, cnStates[whichStates,"s"]), sum) * nrow(L)/sum(!is.na(h) & !is.na(altCount) &! is.na(tumDepth)))[match(round(clusters$proportion,3), round(unique(cfi),3))]
				if(!any(is.na(p) | is.nan(p)))
					power.c <- power.c + p 
			}
			# Tail probs
			pMutCNTail <- matrix(sapply(pmin(cnStates[whichStates,"f"],1), function(pp) ptrbetabinom(altCount[hh],tumDepth[hh],pp, rho=0.01, xmin=pmin(altCount[hh],xmin))), ncol=length(whichStates)) #%*% c(pi.s[cnStates[whichStates,"state"]] * P.m.sX)
			
			
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
			P[[h[i]]] <- cbind(cnStates[whichStates,,drop=FALSE], cfi=cfi[cnStates[whichStates,"state"]], pi.s=pi.s[cnStates[whichStates,"s"]], P.m.sX=P.m.sX, power.s=power.s[cnStates[whichStates,"s"]], power.m.s = power.m.s,
					majCN=majcni[cnStates[whichStates,"state"]], minCN=mincni[cnStates[whichStates,"state"]], 
					majDelta = majdelta[cnStates[whichStates,"state"]], minDelta = mindelta[cnStates[whichStates,"state"]], 
					clonalFlag=clonalFlag[cnStates[whichStates,"state"]], subclonalGainFlag=subclonalGainFlag[cnStates[whichStates,"state"]], mixFlag=mixFlag[cnStates[whichStates,"state"]])
			if(H[i] != h[i]) P[[H[[i]]]] <- P[[h[i]]]

			w <- apply(P.sm.x, 1, function(x) if(any(is.na(x))) NA else which.max(x) )
			if(all(is.na(w))) next
			
			D[hh, "pSub"] <- rowSums(P.sm.x[, !cnStates[whichStates,"state"] %in% which(clonalFlag), drop=FALSE])
			D[hh, "pGain"] <- rowSums(P.sm.x[, cnStates[whichStates,"state"] %in% which(clonalFlag) & cnStates[whichStates,"m"] > 1.00001 + majdelta[cnStates[whichStates,"state"]] + mindelta[cnStates[whichStates,"state"]], drop=FALSE])
			#D[hh, "pSingle"] <- rowSums(P.sm.x[, cnStates[1:k,"state"] %in% which(clonalFlag) & cnStates[1:k,"m"]<=1, drop=FALSE])
			D[hh, "pSingle"] <-  1 - D[hh, "pSub"] - D[hh, "pGain"]			
			
			D[hh,"MutCN"]  <- cnStates[whichStates[w],"m"]
			D[hh,"MutDeltaCN"]  <- majdelta[cnStates[whichStates[w],"state"]] + mindelta[cnStates[whichStates[w],"state"]]
			D[hh,"MinCN"] <- minanc
			D[hh,"MajCN"] <- majanc
			D[hh,"MinDerCN"] <- minder
			D[hh,"MajDerCN"] <- majder
			
			D[hh,"CNF"]  <- cfi[cnStates[whichStates[w],"state"]] 
			D[hh,"pMutCN"] <- sapply(seq_along(w), function(i) P.sm.x[i,w[i]])
			D[hh,"pMutCNTail"] <- sapply(seq_along(w), function(i) pMutCNTail[i,w[i]])
		}
		
		
		#L <- dbinom(altCount[i],tumDepth[i],pmin(cnStates[1:k,"f"],1)) + .Machine$double.eps
		#post <- L * piState[cnStates[1:k,"state"]] * cnStates[1:k,"pi.m.s"]
		#post <- post/sum(post)
#		post <- P.sm.x[hh==i,]
#		if(any(is.nan(post) | is.na(post))) next
#		D[i,"pSub"] <- sum(post[cnStates[1:k,"state"]!=which.max(cfi)])
#		D[i,"pGain"] <- sum(post[cnStates[1:k,"state"]==which.max(cfi) & cnStates[1:k,"m"]>1])
#		D[i,"pSingle"] <- sum(post[cnStates[1:k,"state"]==which.max(cfi) & cnStates[1:k,"m"]<=1])
#		
#		w <- which.max(post)
#		#idx <- as.numeric(strsplit(names(prob[w]), ":")[[1]])
#		#names(prob) <- NULL
#		D[i,"MutCN"] <- cnStates[w,"m"]
#		D[i,"MinCN"] <- mincni[cnStates[w,"state"]]
#		D[i,"MajCN"] <- majcni[cnStates[w,"state"]]
#		D[i,"CNF"] <- cfi[cnStates[w,"state"]] 
#		D[i,"pMutCN"] <- post[w]
		
	}
	return(list(D=D,P=P, power.c=power.c))
}

mcnHeader <- function() {
	DataFrame(Number=c(1,1,1,1,1,1,1,1,".",1,1,1,1),Type=c("Float","Float","Integer","Integer","Integer","Integer","Float","Float","Integer","Float","Float","Float","Float"), 
			Description=c("Mutation copy number","Change in MutCN between ancestral and derived state","Major copy number (ancestral)","Minor copy number (ancestral)","Major copy number (derived)","Minor copy number (derived)","Copy number frequency (relative to all cancer cells)", "MutCN probability","BB segment ids","Posterior prob: Early clonal","Posterior prob: Late clonal","Posterior prob: Subclonal", "Tail prob of mixture model"),
			row.names=c("MutCN","MutDeltaCN","MajCN","MinCN","MajDerCN","MinDerCN","CNF","pMutCN","CNID","pGain","pSingle","pSub", "pMutCNTail"))
}

addMutCn <- function(vcf, bb=allBB[[meta(header(vcf))["ID",]]], clusters=allClusters[[meta(header(vcf))["ID",]]]){
	i = header(vcf)@header$INFO
	exptData(vcf)$header@header$INFO <- rbind(i, mcnHeader())
	info(vcf) <- cbind(info(vcf), computeMutCn(vcf, bb, clusters)$D)
	return(vcf)
}

classifyMutations <- function(vcf, reclassify=c("missing","all","none")) {
	reclassify <- match.arg(reclassify)
	if(nrow(vcf) ==0 )
		return(factor(NULL, levels=c("clonal [early]", "clonal [late]", "clonal [NA]", "subclonal")))
	i <- info(vcf)
	.clsfy <- function(i) {
		cls <- i$CLS
		if(reclassify %in% c("missing", "none")){
			if(all(unique(cls) %in% c("early", "late", "clonal", "subclonal")))
				cls <- factor(cls, levels=c("early", "late", "clonal", "subclonal"), labels=c("clonal [early]", "clonal [late]", "clonal [NA]", "subclonal"))
			cls <- as.character(cls)
			cls[cls=="NA"] <- NA
			if(reclassify=="missing" & any(is.na(cls)))
				cls[is.na(cls)] <- paste(factor(apply(as.matrix(i[is.na(cls), c("pGain","pSingle","pSub")]), 1, function(x) if(all(is.na(x))) NA else which.max(x)), levels=1:3, labels=c("clonal [early]", "clonal [late]","subclonal"))) ## reclassify missing
		}else{
			cls <- paste(factor(apply(as.matrix(i[, c("pGain","pSingle","pSub")]), 1, function(x) if(all(is.na(x))) NA else which.max(x)), levels=1:3, labels=c("clonal [early]", "clonal [late]","subclonal"))) ## reclassify missing
			
		}
		cls[i$pGain==0 & cls!="subclonal"] <- "clonal [NA]"
		if(!is.null(i$MajCN))
			cls[cls!="subclonal" & (i$MajCN == 1 | i$MinCN == 1) & abs(i$MutCN - i$MutDeltaCN -1) <= 0.0001] <- "clonal [NA]"
		cls <- factor(cls, levels=c("clonal [early]", "clonal [late]", "clonal [NA]", "subclonal"))
	}
	cls <- .clsfy(i = i)
	return(cls)
}

posteriorMutCN <- function(x,n, cnStates, xmin=3, rho=0.01){
	whichStates <- 1:nrow(cnStates)
	L <- matrix(sapply(pmin(cnStates[whichStates,"f"],1), function(pp) dtrbetabinom(x,n,pp, rho=rho, xmin=pmin(x,xmin)) + .Machine$double.eps), ncol=length(whichStates))
	P.xsm <- L * rep(cnStates[whichStates,"pi.s"] * cnStates[whichStates,"P.m.sX"], each=nrow(L)) # P(X,s,m)
	P.sm.x <- P.xsm/rowSums(P.xsm) # P(s,m|Xi)
	return(P.sm.x)
}

getGenotype <- function(vcf, reclassify='missing', ...){
	cls <- classifyMutations(vcf = vcf, reclassify=reclassify)
	t <- info(vcf)$TCN
	if(is.null(t))
		t <- info(vcf)$MinCN + info(vcf)$MajCN
	hom <- factor(info(vcf)$MutCN==t, levels=c(TRUE,FALSE))
	dg <- factor(unlist(info(vcf)$DG), levels=as.character(cancerGenes))
	table(gene=dg, class=cls, homozygous=hom, ...)
}

tryExceptNull <- function(x) if(class(x)=="try-error") GRanges() else x

loadVcf <- function(ID){
	file <- dir(vcfPath, pattern=paste0(ID, ".+somatic.snv_mnv.TNC.vcf.bgz$"), full.names=TRUE)
	pos <- loadPositions(ID)
	vcf <- readVcf(file, genome="GRCh37") #, param=ScanVcfParam(which=pos))
	f <- findOverlaps(pos, vcf, select="first")
	vcf <- vcf[na.omit(f)]
	vcf <- addDriver(vcf, cancerGenes)
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


bbplot <- function(bb){
	s <- c(1:22, "X","Y")
	l <- as.numeric(width(refLengths[seqnames(refLengths) %in% s]))
	names(l) <- s
	plot(NA,NA, xlab="",ylab="",xlim=c(0,sum(l)), ylim=c(0,max(max(bb$copy_number, na.rm=TRUE))), xaxt="n")
	c <- cumsum(l)[-length(l)]
	axis(side=1, at=c, labels=rep('', length(l)-1))
	mtext(side=1, at= cumsum(l) - l/2, text=names(l), line=2)
	abline(v=c, lty=3)
	x0 <- start(bb) + cumsum(l)[as.character(seqnames(bb))] - l[as.character(seqnames(bb))]
	x1 <- end(bb) + cumsum(l)[as.character(seqnames(bb))] - l[as.character(seqnames(bb))]
	segments(x0=x0, bb$major_cn ,x1, bb$major_cn, col=2, lwd=5* bb$clonal_frequency)
	segments(x0=x0, bb$minor_cn,x1, bb$minor_cn, col=3, lwd=5* bb$clonal_frequency)
	segments(x0=x0, bb$copy_number,x1, bb$copy_number, col=1, lwd=5* bb$clonal_frequency)
	cv <- coverage(bb)
	cv <- cv[s[s%in%names(cv)]]
	par(xpd=NA)
	for(n in names(cv)){
		cc <- cv[[n]]
		segments(start(cc) + cumsum(l)[n] - l[n] ,-runValue(cc)/2,end(cc)+ cumsum(l)[n] - l[n], -runValue(cc)/2, col=4)
	}
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
library(VariantAnnotation)
drivers <- read.table("/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/final/driver/pcawg_whitelist_coding_drivers_v1_sep302016.txt", header=TRUE, sep="\t")
r <- DNAStringSet(drivers$ref)
a <- DNAStringSet(drivers$alt)
driVers <- VRanges(seqnames = drivers$chr, ranges=IRanges(drivers$pos, width =  width(r)), ref=r, alt=a, sampleNames  = drivers$sample_id)
mcols(driVers) <- drivers[,-c(1,3,4,5,6)]


addFinalDriver <- function(vcf, driVers){
	i = header(vcf)@header$INFO
	exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="String",Description="Driver gene", row.names="DG"))
	if(nrow(vcf)==0){
		info(vcf)$DG <- factor(rep(NA,nrow(vcf)), levels = levels(d$gene))
		return(vcf)
	}
	ID <- meta(header(vcf))["ID",1]
	d <- driVers[sampleNames(driVers)==ID]
	overlaps <- findOverlaps(vcf, d)
	g <- factor(rep(NA,nrow(vcf)), levels = levels(d$gene))
	g[queryHits(overlaps)] <- d$gene[subjectHits(overlaps)]
	info(vcf)$DG <- g
	return(vcf)
}