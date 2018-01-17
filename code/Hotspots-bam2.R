# TODO: Add comment
# 
# Author: mg14
###############################################################################

library(deepSNV)

args <- commandArgs(trailingOnly = TRUE)

setwd("/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/mg14/bam/tumour")


i <- as.numeric(args[1])

tab <- read.table("Protected_sites_shearwater_normalpanel_after_filteringSNPs.txt", sep="\t", header=FALSE)
sites <- GRanges(tab$V1, IRanges(tab$V2, tab$V2))

files = dir(".", pattern=".bam$")
f = files[i]

counts = loadAllData(files=f, regions=sites)

save(counts, file=paste0("../../scratch/hotspots/",i,".RData"))

