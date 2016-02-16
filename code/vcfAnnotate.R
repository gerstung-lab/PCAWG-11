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


dpFiles <- dir(dpPath, pattern="_subclonal_structure.txt", recursive=TRUE)

sampleIds <- sub("_mutation_assignments.txt.gz","",dir(dpPath, pattern="_mutation_assignments.txt.gz"))
sampleIds <- intersect(sampleIds, sub("\\..+","",dir(vcfPath, pattern=".bgz$")))

s <- strsplit(vcfFileName,"/")[[1]]
ID <- sub("\\..+", "", s[length(s)])


print(ID)

clusters <- loadClusters(ID)

bb <- loadBB(ID) 

# Load
vcf <- loadVcf(ID)


# Add ID
meta(header(vcf)) <- rbind(meta(header(vcf)), DataFrame(Value=ID, row.names="ID"))

#' Add mutation copy numbers
vcf <-  addMutCn(vcf, bb)


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
bgzip(fnew)