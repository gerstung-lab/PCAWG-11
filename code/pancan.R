vcfFiles = dir("/nfs/users/nfs_c/cgppipe/pancancer/workspace/sd11/variant_calling_pilot_50/deliverables/2014_10_27_train1_results_Moritz_bugfix/vcf_annotated/", pattern=".vcf.gz$", full.names=TRUE)

r = "/lustre/scratch112/sanger/cgppipe/PanCancerReference/genome.fa.gz" #meta(header(v))["reference",]

for(f in vcfFiles){
 v=readVcf(f, genome="GRCh37")
 s=scanFa(file=r, resize(granges(v), 3,fix="center"))
 t=paste(as.character(unlist(alt(v))),as.character( s),sep="@")
 ct=grepl("(A.CG)|(T..CG)", t)
 print(table(info(v)$DPCL, ct))
}

vcfFiles = dir("vagrent", pattern="TNC.vcf.gz$", full.names=TRUE)
vcf <- readVcf(vcfFiles[3], genome="GRCh37")

vcfList <- sapply(vcfFiles, readVcf, genome="GRCh37")

genes <- lapply(vcfList, function(vcf){
			na.omit(sapply(strsplit(grep("non_synonymous_codon",info(vcf)$VW, value=TRUE), "\\|"), `[`,1))
		})

for(vcf in vcfList){
			print(table(vcf@info$DPCL, grepl("(A.CG)|(T..CG)", paste(as.character(unlist(alt(vcf))),vcf@info$TNC))))
		}
		
t <- sapply(vcfList, function(vcf){
			w <- grepl("TP53", info(vcf)$VW) & sapply(info(vcf)$VC, function(x) if(length(x) ==0) FALSE  else x %in% c('nonsense','missense'))
			if(sum(w)>0)
			info(vcf)$DPCL[w]
		else 0
		})

p <- sapply(vcfList, function(vcf){
			w <- grepl("PTEN", info(vcf)$VW) & sapply(info(vcf)$VC, function(x) if(length(x) ==0) FALSE  else x %in% c('nonsense','missense'))
			if(sum(w)>0)
				info(vcf)$DPCL[w]
			else 0
		})

trunk <- sapply(vcfList, function(vcf){
			which.max(table(info(vcf)$DPCL))
		})

data.frame(trunk=trunk, TP53=t, PTEN=p, row.names=NULL)