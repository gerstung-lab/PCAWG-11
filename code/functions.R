# TODO: Add comment
# 
# Author: mg14
###############################################################################

library(VariantAnnotation)

vcfPath <- '/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/variant_calling_pilot_64/OICR_Sanger_Core'
dpPath <- '/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/pilot_64/sanger_analysis/deliverables/dirichlet_clustering_002_calibration_final/sanger/data/'

sampleIds <- sub("_DPou.+","",dir(dpPath))

loadTree <- function(ID){
 	file <- paste0(dpPath,"/",ID,"_DPoutput_1250iters_250burnin/",ID,"_optimaInfo.txt")
	read.table(file, header=TRUE, sep="\t")
}

loadVcf <- function(ID){
	file <- dir(vcfPath, pattern=paste0(ID, ".+somatic.snv_mnv.vcf.gz$"), full.names=TRUE)
	readVcf(file, genome="GRCh37")
}

load
