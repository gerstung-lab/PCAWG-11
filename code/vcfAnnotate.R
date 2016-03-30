# TODO: Add comment
# 
# Author: mg14
###############################################################################


args <- commandArgs(trailingOnly = TRUE)

source("functions.R")

vcfFileName <- args[1]

library(VariantAnnotation)
library(Matrix)
library(CoxHD)
library(igraph)

r = "/lustre/scratch112/sanger/cgppipe/PanCancerReference/genome.fa.gz" #meta(header(v))["reference",]

dpFiles <- dir(dpPath, pattern="_subclonal_structure.txt", recursive=TRUE)

sampleIds <- sub("_mutation_assignments.txt.gz","",dir(dpPath, pattern="_mutation_assignments.txt.gz"))
sampleIds <- intersect(sampleIds, sub("\\..+","",dir(vcfPath, pattern=".bgz$")))

s <- strsplit(vcfFileName,"/")[[1]]
ID <- sub("\\..+", "", s[length(s)])


print(ID)

clusters <- loadClusters(ID)

if(all(is.na(purityPloidy[ID,]))) # Missing purity
	purityPloidy[ID,] <- c(max(clusters$proportion),NA)


# Load BB
bb <- loadBB(ID)
if(length(bb)==0){ # Missing BB, use consensus CN
	cnPath <- "/nfs/users/nfs_c/cgppipe/pancancer/workspace/mg14/cn/consensus_001/all"
	file <- paste0(cnPath, "/",ID,"_segments.txt")
	tab <- read.table(file, header=TRUE, sep='\t')
	bb <- GRanges(tab$chromosome, IRanges(tab$start, tab$end), strand="*", tab[-3:-1])
	#m <- read.table(gzfile(paste0(basePath, "/0_multiplicity/",ID,"_multiplicity.txt.gz")), header=TRUE)
	#bb <- GRanges(m$chr, IRanges(m$pos, width=1), copy_number=m$tumour_copynumber, major_cn=m$nMaj1, minor_cn=m$nMin1, clonal_frequency=purityPloidy[ID,'purity'])
	#meta(header(vcf)) <- rbind(meta(header(vcf)), DataFrame(Value="False", row.names="Battenberg"))
}
	

# Load vcf
vcf <- readVcf(vcfFileName, genome="GRCh37") #, param=ScanVcfParam(which=pos))

# Load assignments
pos <- loadPositions(ID)
f <- findOverlaps(pos, vcf, select="first")
vcf <- vcf[na.omit(f)]
i = header(vcf)@header$INFO
exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="Integer",Description="DP cluster", row.names="DPC"))
i = header(vcf)@header$INFO
exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="Float",Description="DP cluster probability", row.names="DPP"))
info(vcf)$DPC <- pos$cluster[!is.na(f)]
info(vcf)$DPP <- pos$likelihood[!is.na(f)]	

# Add driver genes
vcf <- addDriver(vcf, cancerGenes)

# Add ID
meta(header(vcf)) <- rbind(meta(header(vcf)), DataFrame(Value=ID, row.names="ID"))

# Add TNC
if(!"TNC" %in% rownames(header(vcf)@header$INFO)){
    tnc=scanFa(file=r, resize(granges(vcf), 3,fix="center"))
    i = header(vcf)@header$INFO
    exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="String",Description="Trinucleotide context", row.names="TNC"))
    info(vcf)$TNC <- as.character(tnc)
}

#' Add mutation copy numbers
vcf <-  addMutCn(vcf, bb, clusters)


#' Remove spurious clusters
info(vcf)$DPC[!info(vcf)$DPC %in% clusters$cluster[clusters$proportion < 1] ] <- NA

#' Classify mutations
class <- rep(2,nrow(vcf))
class[info(vcf)$DPC < max(clusters$cluster[clusters$proportion < 1])] <- 3
class[info(vcf)$MCN > 1 & class==2] <- 1
class <- factor(class, levels=1:3, labels=c("early","late","subclonal"))
class[is.na(info(vcf)$DPC) | is.na(info(vcf)$MCN)] <- NA
info(vcf)$CLS <- class
info(header(vcf)) <- rbind(info(header(vcf)), DataFrame(Number="1",Type="String",Description="Mutation classification {early, late, subclonal}", row.names="CLS"))

#' Save output
fnew <- sub(".vcf",".complete_annotation.vcf",vcfFileName)
writeVcf(vcf, file=fnew)
bgzip(fnew, overwrite=TRUE)
