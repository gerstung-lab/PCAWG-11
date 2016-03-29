library(VariantAnnotation)
args <- commandArgs(trailingOnly = TRUE)
vcfFiles <- args[1]

r = "/lustre/scratch112/sanger/cgppipe/PanCancerReference/genome.fa.gz" #meta(header(v))["reference",]

for(f in vcfFiles){
	vcf=readVcf(f, genome="GRCh37")
	if(!"TNC" %in% rownames(header(vcf)@header$INFO)){ 
		tnc=scanFa(file=r, resize(granges(vcf), 3,fix="center"))
		i = header(vcf)@header$INFO
		exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=1,Type="String",Description="Trinucleotide context", row.names="TNC"))
		info(vcf)$TNC <- as.character(tnc)
		fnew <- sub(".vcf",".TNC.vcf",f)
		writeVcf(vcf, file=fnew)
		bgzip(fnew)
		print(paste("Added TNC to", f))
	}
}
