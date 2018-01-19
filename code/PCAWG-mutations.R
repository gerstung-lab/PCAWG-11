#' # More signatures
#' Get gene models
library(rtracklayer)
g <- "../ref/gencode.v19.annotation.gtf.gz"
gencode <- import(g, format="gtf")
save(gencode, file=sub(".gz",".RData",g))

genes <- gencode[gencode$type=="gene"]
seqlevels(genes) <- sub("chr","", seqlevels(genes))

transStrand <- reduce(genes)
transStrand <- setdiff(transStrand, intersect(transStrand[strand(transStrand)=="+"],transStrand[strand(transStrand)=="-"], ignore.strand=TRUE))
gr <- scanFaIndex(refFile)[1:24]
transStrand <- sort(c(transStrand, setdiff(gr, transStrand, ignore.strand=TRUE)))
transStrand <- transStrand[seqnames(transStrand) %in% c(1:22, "X","Y")]
seqlevels(transStrand) <- c(1:22, "X","Y")

getTranscriptStrandOld <- function(vcf, genes){
	f <- findOverlaps(vcf, genes)
	s <- as.factor(strand(genes))[selectHits(f, select="first")]
	s[is.na(s) | selectHits(f,"count") > 1] <- "*"
	return(s)
}

getTranscriptStrand <- function(vcf, transStrand){
	f <- findOverlaps(vcf, transStrand)
	s <- as.factor(strand(transStrand))[selectHits(f, select="first")]
	s[is.na(s) | selectHits(f,"count") > 1] <- "*"
	return(s)
}

getTrinucleotideSubs <- function(vcf) {
	tnc <- DNAStringSet(info(vcf)$TNC)
	s <- paste0(substr(tnc,1,1),"[",ref(vcf), ">",unlist(alt(vcf)),"]", substr(tnc,3,3))
	n <- c("A","C","G","T")
	f <- paste0(rep(n, each=4), "[", rep(n, each=96/2), ">", c(rep(c("C","G","T"), each=48/3),rep(c("A","G","T"), each=48/3),rep(c("A","C","T"), each=48/3), rep(c("A","C","G"), each=48/3)), "]", n)
	s <- factor(s, levels=f)
}

isClustered0 <- function(vcf, max.dist=1000){
	d <-  diff(start(vcf))
	w <- (d > 1 & d < max.dist)
	c(FALSE, w) | c(w, FALSE)
}


dbetageom <- function(x, alpha, beta, mu=NULL, log=FALSE) {
	if(!is.null(mu))
		beta <- alpha * (1-mu)/mu
	lp <- lbeta(alpha+1,beta + x -1) - lbeta(alpha,beta)
	if(log)
		return(lp)
	else
		exp(lp)
}
	
isClustered <- function(vcf, p=1e-3, q=0.1, r=100, alpha=1){
	d <-  diff(start(vcf))
	w <- d > 1 & diff(as.numeric(seqnames(vcf))) == 0
#	p <- 1e-3 # P N>Kat
#	q <- 0.05 # P Kat>N
	P <- matrix(c(1-p,p,q, 1-q), ncol=2, byrow=TRUE) # Transition matrix
	p0 <- c(1,0)
	s <- c(mean(d[w]), r)
	dw <- d[w]
	l <- length(dw)
	T1 <- T2 <- matrix(0,ncol=l, nrow=2)
	T1[,1] <- log(c(q/(q+p), p/(q+p)))
	lP <- log(P)
	dg <- rbind(dbetageom(dw, alpha=alpha, mu = 1/s[1], log=TRUE),#dgeom(dw, prob=1/s[1], log=TRUE), 
			dgeom(dw, prob=1/s[2], log=TRUE))
	for(i in 2:l){
		x <- ((T1[,i-1] + lP) + dg[,i])
		T2[1,i] <- (x[1,1] < x[2,1])+1
		T2[2,i] <- (x[1,2] < x[2,2])+1
		T1[1,i] <- x[T2[1,i],1]
		T1[2,i] <- x[T2[2,i],2]
	}
	z <- numeric(l)
	z[l] <- 1
	for(i in l:2){
		z[i-1] <- T2[z[i],i]
	}
	k <- numeric(nrow(vcf))
	k[-1][w][-1] <- z[-l]-1
	k[-nrow(vcf)][w][-1] <- (z[-l]-1) | k[-nrow(vcf)][w][-1]
	
	# Other clustered
	pbg <- cumsum(dbetageom(seq(1,max(dw)), alpha=alpha, mu=1/s[1]))
	pc <- pbg[dw]
	#pc <- pgeom(dw, prob=1/s[1], lower.tail=TRUE)
	qc <- p.adjust(pc, "BH") < 0.1
	
	cl <- numeric(nrow(vcf))
	cl[-1][w] <- qc
	cl[-nrow(vcf)][w] <- qc | cl[-nrow(vcf)][w]
	
	clk <- factor(pmin(cl + 2*k,2), levels=0:2, labels=c("None","Clustered","Kataegis"))
	return(clk)
}

#' Replication times
repTimeGm12878 <- import("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq//wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig", format="BigWig")
repTimeK569 <- import("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig", format="BigWig")
repTimeHela <- import("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHelas3WaveSignalRep1.bigWig", format="BigWig")
repTimeHuvec <- import("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig", format="BigWig")
repTimeHepg2 <- import("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig", format="BigWig")

repTime <- granges(repTimeGm12878)
mcols(repTime) <- DataFrame(mg12878=repTimeGm12878$score, hela=repTimeHela$score, k562=repTimeK569$score, huvec=repTimeHuvec$score, hepg2=repTimeHepg2$score)
ndiff <- function(x) (c(NA,x[-length(x)]) - c(x[-1],NA))
t <- sapply(mcols(repTime), ndiff)
m <- rowMeans(t)
s <- rowMeans(sqrt((t-m)^2))
table(abs(m)-2*s > 0)

repStrand <- granges(repTime)
strand(repStrand)[which(m<0)] <- "+" 
strand(repStrand)[which(m>0)] <- "-" 
strand(repStrand)[which(abs(m)-2*s < 0)] <- "*" 
seqlevels(repStrand) <- sub("chr","", seqlevels(repStrand))

save(repStrand, file="../ref/repStrand.RData")

#' Some functions
getRepStrand <- function(vcf, repStrand){
	f <- findOverlaps(vcf, repStrand)
	s <- as.factor(strand(repStrand))[selectHits(f, select="first")]
	return(s)
}

tableSubsStrandRep <- function(vcf, genes, repStrand){
	s <- getTrinucleotideSubs(vcf = vcf)
	t <- getTranscriptStrand(vcf, genes)
	r <- getRepStrand(vcf, repStrand)
	table(sub=s, ts=t, rs=r)
}

isMNV <- function(vcf) {
	d <- diff(start(vcf)) == 1 & abs(diff(getAltCount(vcf) )) <= 1
	w <- c(FALSE, d) | c(d, FALSE)
	return(w)
}

tableSubsTransRepClust <- function(vcf, ts, repStrand){	
	w <- !isMNV(vcf)
	s <- getTrinucleotideSubs(vcf = vcf)
	t <- getTranscriptStrand(vcf, ts)
	c <- isClustered(vcf)
	r <- getRepStrand(vcf, repStrand)
	t <- table(sub=s[w], ts=t[w], rs=r[w], cl=c[w])
	l <- levels(tncToPyrimidine(vcf))
	n <- DNAStringSet(gsub("\\[|\\]|>","",l))
	r <- paste0(as.character(complement(DNAStringSet(substr(n,4,4)))), "[", as.character(complement(DNAStringSet(substr(n,2,2)))),">",as.character(complement(DNAStringSet(substr(n,3,3)))),"]",as.character(complement(DNAStringSet(substr(n,1,1)))))
	t[c(l,r),,,]
}

{
allSubsTransRepClust <- sapply(finalSnv, function(vcf) tableSubsTransRepClust(vcf, transStrand, repStrand), simplify='array')
save(allSubsTransRepClust, file="allSubsTransRepClust3.RData")

library(rhdf5)
h5createFile("allSubsTransRepClust3.h5")
h5write(allSubsTransRepClust, "allSubsTransRepClust3.h5","allSubsTransRepClust")
}

#' A few plots
pdf("Kataegis.pdf", 12,4)
for(i in 1:100)
{
	k <- isClustered(finalSnv[[i]])
	d <- diff(start(finalSnv[[i]]))
	o <- d < 1000 
	plot(d, log='y', col=k[-1], main=names(finalSnv)[i])
	legend("bottomleft", legend=table(k), col=1:3, pch=1)
	print(i)
}
dev.off()

#' Some checks
finalDrivers$samples[which(finalDrivers$gene=="POLE")]

vcf <- finalSnv[["2df02f2b-9f1c-4249-b3b4-b03079cd97d9"]]
s <- getTrinucleotideSubs(vcf = vcf)
r <- getRepStrand(vcf, repStrand)
t <- table(s,r)
m <- match(c(l,r), rownames(t))
barplot(t[m,1:2], beside=TRUE)

#' # TNC context
rs <- repStrand[strand(repStrand)!="*"]
rs <- sort(c(rs, setdiff(gr, rs, ignore.strand=TRUE)))

reference <- scanFa(refFile)
names(reference)[1:24] <- c(1:22,"X","Y")

.tnc <- function(ts, reference){
	rowSums(sapply(seqlevels(ts), function(i)
						colSums(trinucleotideFrequency(Views(reference[[i]], ranges(ts[seqnames(ts)==i]))))))
}


tncTransRep <- sapply(c("+","-","*"), function(x) sapply(c("+","-","*"), function(y) {
						z <- subsetByOverlaps(transStrand[strand(transStrand)==y], rs[strand(rs)==x], ignore.strand=TRUE)
						.tnc(z, reference)
					}), simplify='array')
m <- match(gsub("\\[|\\]|>.","",a$sub[1:192]), rownames(tncTransRep))
write.table(matrix(tncTransRep[m,,]), file="tncTransRep.txt",  col.names=TRUE,sep="\t", quote=FALSE, row.names=FALSE)

#' Sanity check that dims are right.
x <- matrix(mg14:::asum(tncTransRep,1), nrow=9)
y <- matrix(mg14:::asum(allSubsStrand, 1), nrow=9)
quantile(cor(x,y))
x <- matrix(t(mg14:::asum(tncTransRep,1)), nrow=9)
quantile(cor(x,y))


#' # Other mutation types
NUCLEOTIDES <- c("A","C","G","T")
d <- as.character(t(outer(c("C","T","A","G"),NUCLEOTIDES,paste, sep="")))
e <- as.character(reverseComplement(DNAStringSet(d)))

DINUCLEOTIDES <- character()
while(length(d)>0){
	dd <- d[1]
	ee <- e[1]
	DINUCLEOTIDES <- c(DINUCLEOTIDES,dd)
	w <- c(1, which(d==ee))
	d <- d[-w]
	e <- e[-w]
}

MNV_LEVELS <- as.character(sapply(strsplit(DINUCLEOTIDES,""), function(d) paste(paste(d, collapse=""), outer(setdiff(NUCLEOTIDES,d[1]), setdiff(NUCLEOTIDES,d[2]), paste, sep=""), sep=">")))


getMNV <- function(vcf){
	u <- reduce(granges(vcf[which(isMNV(vcf))]))
	u <- u[width(u)>1]
	if(length(u)==0) return(table(factor(NULL, levels=MNV_LEVELS), useNA='a'))
	r <- as.character(ref(vcf)[vcf %over% u])
	h <- subjectHits(findOverlaps(granges(vcf),u))
	a <- as.character(unlist(alt(vcf))[vcf %over% u])
	rr <- sapply(split(r, h), paste, collapse="")
	aa <- sapply(split(a, h), paste, collapse="")
	
	w <- which(!rr %in% DINUCLEOTIDES)
	rr[w] <- as.character(reverseComplement(DNAStringSet(rr[w])))
	aa[w] <- as.character(reverseComplement(DNAStringSet(aa[w])))
	t <- table(factor(paste(rr,aa,sep=">"), levels=MNV_LEVELS), useNA='a')
	names(t)[length(t)] <- "MNV(other)"
	return(t)
}

mnvTable <- simplify2array(mclapply(finalSnv, getMNV, mc.cores=1))

write.table(t(mnvTable), file="mnvTable.txt", sep="\t")

getIndels <- function(vcf){
	r <- ref(vcf)
	r <- substr(r,2,width(r))
	a <- unlist(alt(vcf))
	a <- substr(a,2,width(a))
	c <- r != "" & a != ""
	w <- which(!r %in% c("C","T", DINUCLEOTIDES))
	r[w] <- as.character(reverseComplement(DNAStringSet(r[w])))
	w <- which(width(r)>2)
	l <- cut(width(r)[w],breaks=c(2:9, seq(10,100,10), Inf))
	levels(l) <- c(3:10, levels(l)[-(1:8)])
	r[w] <- as.character(l)
	w <- which(!a %in% c("C","T", DINUCLEOTIDES))
	a[w] <- as.character(reverseComplement(DNAStringSet(a[w])))
	w <- which(width(a)>2)
	l <- cut(width(a)[w],breaks=c(2:9, seq(10,100,10), Inf))
	levels(l) <- c(3:10, levels(l)[-(1:8)])
	a[w] <- as.character(l)
	indel <- as.character(factor((r=="") + 2*(a==""), levels=0:3, labels=c("indel","ins","del","")))
	indel[indel=="ins"] <- paste("ins", a[indel=="ins"], sep="")
	indel[indel=="del"] <- paste("del", r[indel=="del"], sep="")
	i <- c(t(outer(c("del","ins"), c("C","T", DINUCLEOTIDES, levels(l)), paste, sep="")), "indel")
	t <- table(factor(indel, levels=i), useNA='a')
	names(t)[lenght(t)] <- "indel(other)"
	t
}

indelTable <- sapply(finalIndel, getIndels)

write.table(t(indelTable), file="indelTable.txt", sep="\t")

getSNV <- function(vcf){
	w <- which(!isMNV(vcf))
	a <- unlist(alt(vcf))[w]
	r <- ref(vcf)[w]
	tnc <- DNAStringSet(info(vcf)$TNC)[w]
	rc <- grepl("A|G", r)
	tnc[rc] <- reverseComplement(tnc[rc])
	a[rc] <- reverseComplement(a[rc])
	t <- paste0(substr(tnc,1,1), "[",substr(tnc,2,2), ">",a, "]", substr(tnc,3,3))
	n <- c("A","C","G","T")
	f <- paste0(rep(n, each=4), "[", rep(c("C","T"), each=96/2), ">", c(rep(c("A","G","T"), each=48/3), rep(c("A","C","G"), each=48/3)), "]", n)
	table(factor(t, levels=f))
}

snvTable <- sapply(finalSnv, getSNV)

write.table(t(snvTable), file="snvTable.txt", sep="\t")

svTable <- t(read.table("../data/SV_nmf_input_and_results.20170827.txt", header=TRUE, row.names=1, sep="\t"))
id <- gsub("\\.","-",sub("^[A-Z]+\\.[A-Z]+\\.","", rownames(svTable)))
rownames(svTable) <- id

nrow(svTable)
w <- names(finalSnv)[whiteList] %in% rownames(svTable) | (! sample2donor[names(finalSnv)[whiteList]] %in% sample2donor[rownames(svTable)] & ! duplicated(sample2donor[names(finalSnv)[whiteList]]))

s <- matrix(NA, ncol =ncol(svTable), nrow=ncol(snvTable), dimnames=list(colnames(snvTable), colnames(svTable)))
s[which(w),] <- 0
s[rownames(svTable),] <- svTable

getCNA <- function(bb){
	u <- !duplicated(bb)
	M <- cut(bb$major_cn[u], c(-1:4, 10, Inf))
	m <- cut(bb$minor_cn[u], c(-1:4, 10, Inf))
	levels(M) <- levels(m) <- c(0:4, "5-10",">10")
	c <- cut(width(bb)[u],10^seq(1,6))
	d <- as.data.frame(table(m,M,c))
	f <- d$Freq
	names(f) <- paste(d$m, d$M, d$c, sep=":")
	f
}

cnaTable <- sapply(finalBB, getCNA)
r <- rownames(cnaTable)
c <- sub("\\:\\(.+", "",r)
t <- table(c[rowSums(cnaTable)==0])
w <- !c %in% names(t)[t==5] &! c=="1:1"
sum(cnaTable[w,])

allTable <- cbind(t(snvTable), t(mnvTable), t(indelTable), s, t(cnaTable[w,]))

write.table(allTable, file="allTable.txt", sep="\t")

