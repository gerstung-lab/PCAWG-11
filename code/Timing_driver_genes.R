library(VariantAnnotation)
library(ggplot2)
library(cowplot)
library(abind)
library(boot)
library(reshape2)
library(parallel)

classifyMutations <- function(x, reclassify=c("missing","all","none")) {
  reclassify <- match.arg(reclassify)
  if(nrow(x) ==0 )
    return(factor(NULL, levels=c("clonal [early]", "clonal [late]", "clonal [NA]", "subclonal")))
  if(class(x)=="CollapsedVCF")
    x <- info(x)
  .clsfy <- function(x) {
    cls <- x$CLS
    if(reclassify %in% c("missing", "none") &! is.null(cls)){
      if(all(unique(cls) %in% c("early", "late", "clonal", "subclonal")))
        cls <- factor(cls, levels=c("early", "late", "clonal", "subclonal"), labels=c("clonal [early]", "clonal [late]", "clonal [NA]", "subclonal"))
      cls <- as.character(cls)
      cls[cls=="NA"] <- NA
      if(reclassify=="missing" & any(is.na(cls)))
        cls[is.na(cls)] <- paste(factor(apply(as.matrix(x[is.na(cls), c("pGain","pSingle","pSub")]), 1, function(x) if(all(is.na(x))) NA else which.max(x)), levels=1:3, labels=c("clonal [early]", "clonal [late]","subclonal"))) ## reclassify missing
    }else{
      cls <- paste(factor(apply(as.matrix(x[, c("pGain","pSingle","pSub")]), 1, function(x) if(all(is.na(x))) NA else which.max(x)), levels=1:3, labels=c("clonal [early]", "clonal [late]","subclonal"))) ## reclassify missing
      
    }
    cls[x$pGain==0 & cls!="subclonal"] <- "clonal [NA]"
    if(!is.null(x$MajCN))
      cls[cls!="subclonal" & (x$MajCN == 1 | x$MinCN == 1) & abs(x$MutCN - x$MutDeltaCN -1) <= 0.0001] <- "clonal [NA]"
    cls <- factor(cls, levels=c("clonal [early]", "clonal [late]", "clonal [NA]", "subclonal"))
  }
  cls <- .clsfy(x)
  return(cls)
}

getGenotype <- function(vcf, reclassify='missing', ...){
  w <- c(TRUE,diff(start(vcf)) != 1)
  cls <- classifyMutations(vcf, reclassify=reclassify)
  hom <- factor(paste(info(vcf)$MajCN, info(vcf)$MinCN, sep = ":"), levels=c("1:0","1:1", "2:0", "2:1", "2:2", "3:1", "3:2", "3:3"))
  dg <- factor(unlist(info(vcf)$DG), levels=as.character(CANCERGENES))
  table(gene=dg[w], class=cls[w], CN=hom[w], ...)
}

Genotypes <- function(snv_path, indel_path){
  finalSnv <- list()
  for(f in dir(snv_path, pattern="*.vcf.bgz", full.names=TRUE)){
    vcf <- readVcf(f)
    finalSnv[[f]] <- vcf
  }
  names(finalSnv) <- sub(".conse.+","",dir(snv_path, pattern="*.vcf.bgz", full.names=FALSE))
  finalIndel <- list()
  for( f in dir(indel_path, pattern="*.vcf.bgz", full.names=TRUE)){
    vcfIndel <- tryCatch({readVcf(f)}, error = function(cond) { return(finalIndel[[1]][0,])})
    finalIndel[[f]] <- vcfIndel
  }
  names(finalIndel) <- sub(".conse.+","",dir(indel_path, pattern="*.vcf.RData", full.names=FALSE))
  finalDriversAnnotated <- finalDrivers
  d <- info(finalSnv[[3]])[seq_along(finalDriversAnnotated),19:32]
  
  mcols(finalDriversAnnotated)[colnames(d)] <- d
  finalDriversAnnotated <- finalDriversAnnotated[finalDriversAnnotated$sample %in% names(finalSnv),]
  for(i in seq_along(finalDriversAnnotated)){
    if(finalDriversAnnotated$mut_type[i] %in% c("snv","mnv")){
      v <- finalSnv[[as.character(finalDriversAnnotated$sample[i])]]
    }else{
      v <- finalIndel[[as.character(finalDriversAnnotated$sample[i])]]
    }
    j <- findOverlaps(finalDriversAnnotated[i], v, select='first')
    if(!is.na(j))
      mcols(finalDriversAnnotated)[i,colnames(d)] <- info(v)[j, colnames(d)]
    else
      mcols(finalDriversAnnotated)[i,colnames(d)] <- NA
  }
  
  MC_CORES=1
  
  finalGenotypesSnv <- simplify2array(mclapply(finalSnv, getGenotype, mc.cores=MC_CORES, useNA="always"))
  finalGenotypesIndel <- simplify2array(mclapply(finalIndel, getGenotype, mc.cores=MC_CORES, useNA="always"))
  finalGenotypes <- aperm(abind::abind(subs=finalGenotypesSnv,indels=finalGenotypesIndel, along=5), c(1,5,2,3,4))
  rm(finalGenotypesSnv,finalGenotypesIndel)
  
  #save(finalGenotypes, file ="finalGenotypes.RData")
  return(finalGenotypes)
}

asum <- function(x, dim) apply(x, setdiff(seq_along(dim(x)), dim), sum)

timing_drivers <- function(finalGenotypes, purity_ploidy, tissueColors, sample2donor, specimenData, n_top){
  filter_patients <- purity_ploidy$samplename[purity_ploidy$purity > 0.3]
  filter_patients <- intersect(filter_patients, dimnames(finalGenotypes)[[5]])
  
  topgenes <- asum(finalGenotypes[,,,,filter_patients],c(2,4,5))
  topgenes <- cbind(topgenes, total = apply(topgenes, 1, function(x) sum(x[1],x[2],x[3],x[4])))
  
  gene <- strsplit(dimnames(finalGenotypes)[[1]], "::")
  
  gene <- sapply(gene, function(x) paste(x[2],x[3], sep = "::"))
  gene <- unique(gene[1:(length(gene)-1)])
  topgenes <- t(sapply(as.character(gene), function(x) data.frame(sum(topgenes[grep(paste("::",x,"::", sep = ""), rownames(topgenes)), 1]),
                                                                  sum(topgenes[grep(paste("::",x,"::", sep = ""), rownames(topgenes)), 2]),
                                                                  sum(topgenes[grep(paste("::",x,"::", sep = ""), rownames(topgenes)), 3]),
                                                                  sum(topgenes[grep(paste("::",x,"::", sep = ""), rownames(topgenes)), 4]),
                                                                  sum(topgenes[grep(paste("::",x,"::", sep = ""), rownames(topgenes)), "total"]))))
  
  colnames(topgenes) <- c("clonal [early]", "clonal [late]", "clonal [NA]", "subclonal", "total")
  topgenes <- topgenes[-grep("gencode::RMRP", rownames(topgenes)),]
  topgenes <- data.frame(topgenes[order(as.numeric(topgenes[,"total"]), decreasing = TRUE ),][1:n_top,])
  ttype <- data.frame(patient = names(sample2donor), stringsAsFactors = FALSE)
  ttype$tumour <- sapply(ttype$patient, function(x) unique(as.character(specimenData$histology_abbreviation[specimenData$icgc_donor_id == sample2donor[[x[1]]]]))[1])
  ttype$color <- sapply(ttype$tumour, function(x) tissueColors[x[1]])
  
  count_genes <- data.frame()
  
  r_early <- list()
  r_clonal <- list()
  for (gene in rownames(topgenes)){
    
    lgene <- grep(paste("::",gene,"::", sep = ""), dimnames(finalGenotypes)[[1]], value = TRUE)
    Emut <- asum(finalGenotypes[lgene,,"clonal [early]",,filter_patients], 1:((length(lgene) > 1 )+2))
    Emut <- names(Emut[Emut>0])
    
    Lmut <- asum(finalGenotypes[lgene,, "clonal [late]",,filter_patients], 1:((length(lgene) > 1 )+2))
    Lmut <- names(Lmut[Lmut>0])
    
    if( length(union(Emut,Lmut)) > 5){
      odds <- boot(data = union(Emut,Lmut), statistic = function (x1, x2){
                     Xemut <- sum(finalGenotypes[lgene,,"clonal [early]",,x1[x2]])
                     if(Xemut == 0) if(rbinom(1,1,0.5) == 1) Xemut <- 1
                     
                     Xlmut <- sum(finalGenotypes[lgene,,"clonal [late]",,x1[x2]])
                     if(Xlmut == 0) if(rbinom(1,1,0.5) == 1) Xlmut <- 1
                     
                     Xback <- matrix(numeric(18), nrow = 2)
                     colnames(Xback) <- c("1:0", "1:1", "2:0", "2:1", "2:2", "3:1", "3:2", "3:3", "NA")
                     rownames(Xback) <- c("muts", "patients")
                     
                     for(c_CN in colnames(Xback)){
                       if(c_CN == "NA") c_CN <- 9
                       c_patients <- asum(finalGenotypes[lgene,,c("clonal [early]","clonal [late]"),c_CN,x1[x2]], 1:((length(lgene) > 1 )+2))
                       c_patients <- names(c_patients[c_patients > 0])
                       if(length(c_patients) == 0 ) next
                       if(length(c_patients) == 1 ){
                         Xback["muts",c_CN] <- mean(sum(finalGenotypes[,,"clonal [early]",c_CN,c_patients]) / sum(finalGenotypes[,,c("clonal [early]","clonal [late]"),c_CN,c_patients]))
                       }else{
                           Xback["muts",c_CN] <- mean(asum(finalGenotypes[,,"clonal [early]",c_CN,c_patients],c(1,2)) / asum(finalGenotypes[,,c("clonal [early]","clonal [late]"),c_CN,c_patients],c(1,2,3)))
                       }
                       Xback["patients",c_CN] <- length(c_patients)
                     }
                     ((Xemut / (Xemut + Xlmut)) / (1 - (Xemut / (Xemut + Xlmut)))) / (weighted.mean(Xback["muts",], Xback["patients",]) / (1- weighted.mean(Xback["muts",], Xback["patients",])))
                   },R = 1000)
      
      r_early[[gene]] <- list(length(Emut), length(Lmut), odds$t0, as.numeric(quantile(odds$t, 0.025, na.rm =TRUE)), as.numeric(quantile(odds$t, 0.975, na.rm =TRUE)), mean(odds$t, na.rm = TRUE), median(odds$t, na.rm = TRUE), sum(odds$t <= 1, na.rm = TRUE))
    }
    
    Emut <- asum(finalGenotypes[lgene,, c("clonal [early]","clonal [late]", "clonal [NA]"),,filter_patients], 1:((length(lgene) > 1 )+3))
    Emut <- names(Emut[Emut>0])
    
    Lmut <- asum(finalGenotypes[lgene,, "subclonal",,filter_patients], 1:((length(lgene) > 1 )+2))
    Lmut <- names(Lmut[Lmut>0])
    
    if( length(union(Emut,Lmut)) > 5){
      odds <- boot(data = union(Emut,Lmut), statistic = function (x1, x2){
                     Xemut <- sum(finalGenotypes[lgene,, c("clonal [early]","clonal [late]", "clonal [NA]"),,x1[x2]])
                     if(Xemut == 0) if(rbinom(1,1,0.5) == 1) Xemut <- 1
                     
                     Xlmut <- sum(finalGenotypes[lgene,, "subclonal",,x1[x2]])
                     if(Xlmut == 0) if(rbinom(1,1,0.5) == 1) Xlmut <- 1
                     
                     Xback <- matrix(numeric(18), nrow = 2)
                     colnames(Xback) <- c("1:0", "1:1", "2:0", "2:1", "2:2", "3:1", "3:2", "3:3", "NA")
                     rownames(Xback) <- c("muts", "patients")
                     
                     for(c_CN in colnames(Xback)){
                       if(c_CN == "NA") c_CN <- 9
                       c_patients <- asum(finalGenotypes[lgene,,c("clonal [early]","clonal [late]", "clonal [NA]","subclonal"),c_CN,x1[x2]], 1:((length(lgene) > 1 )+2))
                       c_patients <- names(c_patients[c_patients > 0])
                       if(length(c_patients) == 0) next
                       if(length(c_patients) == 1){
                         Xback["muts",c_CN] <- mean(sum(finalGenotypes[,,c("clonal [early]","clonal [late]", "clonal [NA]"),c_CN,c_patients]) / sum(finalGenotypes[,,c("clonal [early]","clonal [late]","clonal [NA]","subclonal"),c_CN,c_patients]))
                       }else{
                         Xback["muts",c_CN] <- mean(asum(finalGenotypes[,,c("clonal [early]","clonal [late]", "clonal [NA]"),c_CN,c_patients],c(1,2,3)) / asum(finalGenotypes[,,c("clonal [early]","clonal [late]","clonal [NA]","subclonal"),c_CN,c_patients],c(1,2,3)))
                       }
                       Xback["patients",c_CN] <- length(c_patients)
                     }
                     ((Xemut / (Xemut + Xlmut)) / (1 - (Xemut / (Xemut + Xlmut)))) / (weighted.mean(Xback["muts",], Xback["patients",]) / (1- weighted.mean(Xback["muts",], Xback["patients",])))
                   },R = 1000)
      
      r_clonal[[gene]] <- list(length(Emut), length(Lmut), odds$t0, as.numeric(quantile(odds$t, 0.025, na.rm =TRUE)), as.numeric(quantile(odds$t, 0.975, na.rm =TRUE)), mean(odds$t, na.rm = TRUE), median(odds$t, na.rm = TRUE), sum(odds$t <= 1, na.rm = TRUE))
    }
    
    count_genes <- rbind(count_genes, data.frame(Gene = gene, data.frame(table(ttype$color[ttype$patient %in% names(asum(finalGenotypes[lgene,,,,filter_patients], 1:((length(lgene) > 1 )+3))[asum(finalGenotypes[lgene,,,,filter_patients], 1:((length(lgene) > 1 )+3)) > 0])]))))
  }
  
  r_all <- data.frame()
  r_all <- rbind(r_all, data.frame(t(sapply(r_early, function(x) c(x[[3]],x[[4]],x[[5]],x[[8]]))), state = "early", Gene = names(r_early), stringsAsFactors = FALSE))
  r_all <- rbind(r_all, data.frame(t(sapply(r_clonal, function(x) c(x[[3]],x[[4]],x[[5]],x[[8]]))), state = "clonal", Gene = names(r_clonal), stringsAsFactors = FALSE))
  
  r_all$Gene <- gsub("gencode::","",r_all$Gene)
  r_all$Gene <- gsub("::NA","",r_all$Gene)
  rownames(topgenes) <- gsub("gencode::","",rownames(topgenes))
  rownames(topgenes) <- gsub("::NA","",rownames(topgenes))
  count_genes$Gene <- gsub("gencode::","",count_genes$Gene)
  count_genes$Gene <- gsub("::NA","",count_genes$Gene)
  
  
  colnames(r_all)[1] <- "Odds"
  colnames(r_all)[2] <- "Odds_err1"
  colnames(r_all)[3] <- "Odds_err2"
  colnames(r_all)[4] <- "Smaller_1"
  
  r_all$Pvalue <- vapply(r_all$Smaller_1, function(x) as.numeric(x)/10000, as.numeric(1))
  r_all$Pvalue[r_all$Pvalue > 0.5] <- 1 - r_all$Pvalue[r_all$Pvalue > 0.5]
  r_all$Pvalue <- 2*r_all$Pvalue
  r_all$BH <- p.adjust(r_all$Pvalue, method = "fdr")
  
  #### Plotting
 
  r_all$Gene <- factor(r_all$Gene, levels = rownames(topgenes)[order(as.numeric(topgenes[,"total"]), decreasing = TRUE)])
  r_all[is.infinite(r_all$Odds),"Odds"] <- 50
  
  r_genes <- data.frame(Gene = rownames(topgenes))
  for (i in colnames(topgenes)) r_genes[,i] <- unlist(topgenes[,i])
  r_genes <- melt(r_genes, id.vars = c("Gene", "total"))
  r_genes$Gene <- factor(r_genes$Gene, levels = rownames(topgenes)[order(as.numeric(topgenes[,"total"]), decreasing = TRUE)])
  r_genes$variable <- factor(r_genes$variable, levels = c("subclonal", "clonal..NA.", "clonal..late.", "clonal..early."))
  
  count_genes$Gene <-  factor(count_genes$Gene, levels = rownames(topgenes)[order(as.numeric(topgenes[,"total"]), decreasing = TRUE)])
  count_genes$Tissue <- sapply(as.character(count_genes$Var1), function(x) names(tissueColors[tissueColors == x[1]]))
  count_genes$Var1 <- factor(as.character(count_genes$Var1), levels = unique(count_genes$Var1))
  count_genes$Percent <- apply(count_genes, 1, function(x) as.numeric(x["Freq"])/sum(as.numeric(count_genes[count_genes$Gene == x["Gene"],"Freq"])) *100 )
  
  cbPalette <- c("#E41A1C", "#377EB8", "#984EA3", "#4DAF4A")
  
  q1 <- ggplot(r_genes) + theme(panel.background = element_blank(), text = element_text(family = "mono"),
                                plot.title = element_text(hjust = 0.5),
                                axis.title.x=element_blank(),
                                axis.title.y=element_blank(),
                                axis.text.x = element_blank(), #element_text(angle = 45, vjust = 1, hjust = 1),
                                axis.ticks.x = element_blank(),
                                legend.position="none")
  q1 <- q1 + aes(x = Gene, y = value, fill = variable) + geom_bar(stat = "identity") + scale_fill_manual(values=cbPalette)
  q1 <- q1 + scale_y_continuous(expand = c(0.005,0))
  q1
  
  r_all$Odds_err1[r_all$Odds_err1 < 0.1] <- 0.1
  r_all$Odds_err2[r_all$Odds_err2 > 50] <- 50
  
  r_all$Color <- unlist(apply(r_all,1, function(x){
    if(x["state"] == "early"){
      if(as.numeric(x["Odds_err1"]) > 1){
        "#4DAF4A"
      }else if(as.numeric(x["Odds_err2"]) < 1){
        "#984EA3"
      }else{
        "#A5A4A4"
      }
    } else{
      if(as.numeric(x["Odds_err1"]) > 1){
        "#377EB8"
      }else if(as.numeric(x["Odds_err2"]) < 1){
        "#E41A1C"
      }else{
        "#A5A4A4"
      }
    }
  }))
  r_all$Color <- factor(r_all$Color)
  
  q2 <- ggplot(r_all) + theme(panel.background = element_blank(), text = element_text(family = "mono"), 
                              plot.title = element_text(hjust = 0.5),
                              axis.title.x=element_blank(),
                              axis.title.y=element_blank(),
                              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                              legend.position="none")
  q2 <- q2 + aes(x = Gene, y = Odds, ymax = Odds_err2, ymin = Odds_err1, group = factor(state, levels = c("early","clonal")), color = Color)
  q2 <- q2 + geom_linerange(position = position_dodge(width = 0.5), colour = "#A5A4A4") + geom_hline( yintercept = 1) + geom_point(position = position_dodge(width = 0.5), size = 3) 
  q2 <- q2 + scale_y_log10(breaks = c(seq(0.1, 0.9, 0.1), seq(1, 9, 1), seq(10, 100, 10)), limits = c(0.09,60), labels = c(0.1,rep("",8),1,rep("",8),10,rep("",3),50,rep("",4),100), expand = c(0,0))
  q2 <- q2 + scale_color_manual(breaks = levels(r_all$Color), values = levels(r_all$Color))
  q2
  
  tissuePalette <- unique(as.character(count_genes$Var1))
  
  q3 <- ggplot(count_genes) + theme(panel.background = element_blank(), text = element_text(family = "mono"),
                                    plot.title = element_text(hjust = 0.5),
                                    axis.title.x=element_blank(),
                                    axis.title.y=element_blank(),
                                    axis.text.x = element_blank(),
                                    axis.ticks.x = element_blank(),
                                    axis.line.x = element_blank(),
                                    legend.position="none")
  q3 <- q3 + aes(x = Gene, y = Percent, fill = Var1) +  geom_bar(stat = "identity") 
  q3 <- q3 + scale_fill_manual(values= tissuePalette, labels = unique(count_genes$Tissue))
  q3 <- q3 + scale_y_continuous(breaks = c(0,100), labels = c(0,1), expand = c(0,0))
  q3
  
  qall <- plot_grid(q1, q3, q2, nrow = 3, ncol = 1, align = "v", rel_heights = c(5,1,5))
  
  #ggsave("2018_08_odds_driver.pdf", plot = qall, width = 12, height = 7)
  return(list(data = r_all, plot = qall, count_genes = count_genes))
}


timing_tp53 <- function(finalGenotypes, purity_ploidy, tissueColors, sample2donor, specimenData, count_genes){
  gene <- "TP53"
  r_tp53 <- data.frame()
  lgene <- grep(paste("::",gene,"::", sep = ""), dimnames(finalGenotypes)[[1]], value = TRUE)
 
  ttype <- data.frame(patient = names(sample2donor), stringsAsFactors = FALSE)
  ttype$tumour <- sapply(ttype$patient, function(x) unique(as.character(specimenData$histology_abbreviation[specimenData$icgc_donor_id == sample2donor[[x[1]]]]))[1])
  ttype$color <- sapply(ttype$tumour, function(x) tissueColors[x[1]])
 
  r_early <- list()
  r_clonal <- list()
  for (tissue in names(tissueColors)){
    patients <- ttype[ttype$tumour == tissue,"patient"]
    patients <- patients[patients %in% dimnames(finalGenotypes)[[5]]]
    if (length(patients) < 5) next
    Emut <- asum(finalGenotypes[lgene,,"clonal [early]",,patients], 1:((length(lgene) > 1 )+2))
    Emut <- names(Emut[Emut>0])
    Lmut <- asum(finalGenotypes[lgene,,"clonal [late]",,patients], 1:((length(lgene) > 1 )+2))
    Lmut <- names(Lmut[Lmut>0])
    
    if( length(union(Emut,Lmut)) > 5){
      odds <- boot(data = union(Emut,Lmut), statistic = function (x1, x2){ 
                     Xemut <- sum(finalGenotypes[lgene,,"clonal [early]",,x1[x2] ])# + 0.5
                     if(Xemut == 0) if(rbinom(1,1,0.5) == 1) Xemut <- 1
                     Xlmut <- sum(finalGenotypes[lgene,,"clonal [late]",,x1[x2] ])# + 0.5
                     if(Xlmut == 0) if(rbinom(1,1,0.5) == 1) Xlmut <- 1
                     
                     Xback <- matrix(numeric(18), nrow = 2)
                     colnames(Xback) <- c("1:0", "1:1", "2:0", "2:1", "2:2", "3:1", "3:2", "3:3", "NA")
                     rownames(Xback) <- c("muts", "patients")
                     
                     for(c_CN in colnames(Xback)){
                       if(c_CN == "NA") c_CN <- 9
                       c_patients <- asum(finalGenotypes[lgene,,c("clonal [early]","clonal [late]"),c_CN,x1[x2]], 1:((length(lgene) > 1 )+2))
                       c_patients <- names(c_patients[c_patients > 0])
                       if(length(c_patients) == 0 ) next
                       if(length(c_patients) == 1 ){
                           Xback["muts",c_CN] <- mean(sum(finalGenotypes[,,"clonal [early]",c_CN,c_patients]) / sum(finalGenotypes[,,c("clonal [early]","clonal [late]"),c_CN,c_patients]))
                       }else{
                           Xback["muts",c_CN] <- mean(asum(finalGenotypes[,,"clonal [early]",c_CN,c_patients],c(1,2)) / asum(finalGenotypes[,,c("clonal [early]","clonal [late]"),c_CN,c_patients],c(1,2,3)))
                       }
                       Xback["patients",c_CN] <- length(c_patients)
                     }
                     ((Xemut / (Xemut + Xlmut)) / (1 - (Xemut / (Xemut + Xlmut)))) / (weighted.mean(Xback["muts",], Xback["patients",]) / (1- weighted.mean(Xback["muts",], Xback["patients",])))
                   } ,R = 1000)
      r_early[[tissue]] <- list(length(Emut), length(Lmut), odds$t0, as.numeric(quantile(odds$t, 0.025, na.rm = TRUE)), as.numeric(quantile(odds$t, 0.975, na.rm =TRUE)), mean(odds$t, na.rm =TRUE), median(odds$t, na.rm =TRUE))
    }
    
    Emut <- asum(finalGenotypes[lgene,, c("clonal [early]","clonal [late]", "clonal [NA]"),,patients],  1:((length(lgene) > 1 )+3))
    Emut <- names(Emut[Emut>0])
    Lmut <- asum(finalGenotypes[lgene,, "subclonal",,patients],  1:((length(lgene) > 1 )+2))
    Lmut <- names(Lmut[Lmut>0])
    
    if( length(union(Emut,Lmut)) > 5){
      odds <- boot(data = union(Emut,Lmut),
                   statistic = function (x1, x2){ 
                     Xemut <- sum(finalGenotypes[lgene,, c("clonal [early]","clonal [late]", "clonal [NA]"),,x1[x2]]) #+ 0.5
                     if(Xemut == 0) if(rbinom(1,1,0.5) == 1) Xemut <- 1
                     Xlmut <- sum(finalGenotypes[lgene,, "subclonal",,x1[x2]]) #+ 0.5
                     if(Xlmut == 0) if(rbinom(1,1,0.5) == 1) Xlmut <- 1
                     
                     Xback <- matrix(numeric(18), nrow = 2)
                     colnames(Xback) <- c("1:0", "1:1", "2:0", "2:1", "2:2", "3:1", "3:2", "3:3", "NA")
                     rownames(Xback) <- c("muts", "patients")
                     
                     for(c_CN in colnames(Xback)){
                       #for(c_CN in "all"){
                       if(c_CN == "NA") c_CN <- 9
                       c_patients <- asum(finalGenotypes[lgene,,c("clonal [early]","clonal [late]", "clonal [NA]","subclonal"),c_CN,x1[x2]], 1:((length(lgene) > 1 )+2))
                       c_patients <- names(c_patients[c_patients > 0])
                       if(length(c_patients) == 0) next
                       if(length(c_patients) == 1){
                           Xback["muts",c_CN] <- mean(sum(finalGenotypes[,,c("clonal [early]","clonal [late]", "clonal [NA]"),c_CN,c_patients]) / sum(finalGenotypes[,,c("clonal [early]","clonal [late]","clonal [NA]","subclonal"),c_CN,c_patients]))
                       }else{
                           Xback["muts",c_CN] <- mean(asum(finalGenotypes[,,c("clonal [early]","clonal [late]", "clonal [NA]"),c_CN,c_patients],c(1,2,3)) / asum(finalGenotypes[,,c("clonal [early]","clonal [late]","clonal [NA]","subclonal"),c_CN,c_patients],c(1,2,3)))
                       }
                       Xback["patients",c_CN] <- length(c_patients)
                     }
                     ((Xemut / (Xemut + Xlmut)) / (1 - (Xemut / (Xemut + Xlmut)))) / (weighted.mean(Xback["muts",], Xback["patients",]) / (1- weighted.mean(Xback["muts",], Xback["patients",])))
                   } ,R = 1000)
      r_clonal[[tissue]] <- list(length(Emut), length(Lmut), odds$t0, as.numeric(quantile(odds$t, 0.025, na.rm =TRUE)), as.numeric(quantile(odds$t, 0.975, na.rm =TRUE)), mean(odds$t), median(odds$t, na.rm =TRUE))
    }
  }
  
  r_tp53 <- rbind(r_tp53, data.frame(t(sapply(r_early, function(x) c(x[[1]],x[[2]],x[[3]],x[[4]],x[[5]]))), state = "early", Gene = names(r_early), stringsAsFactors = FALSE))
  r_tp53 <- rbind(r_tp53, data.frame(t(sapply(r_clonal, function(x) c(x[[1]],x[[2]],x[[3]],x[[4]],x[[5]]))), state = "clonal", Gene = names(r_clonal), stringsAsFactors = FALSE))
  
  colnames(r_tp53)[1] <- "Pre"
  colnames(r_tp53)[2] <- "Pos"
  colnames(r_tp53)[3] <- "Odds"
  colnames(r_tp53)[4] <- "Odds_err1"
  colnames(r_tp53)[5] <- "Odds_err2"
  
  #### Plotting
  r_tp53$Gene <- paste(r_tp53$Gene, " (", sapply(r_tp53$Gene, function(x) count_genes$Freq[count_genes$Gene == "TP53" & count_genes$Tissue == x]), ")", sep = "")
  r_tp53$Gene <- factor(r_tp53$Gene, levels = unique(r_tp53$Gene)[order(sapply(unique(r_tp53$Gene), function(x) sum(r_tp53[r_tp53$Gene == x, c("Pre", "Pos")])), decreasing = TRUE)])
  r_tp53[is.infinite(r_tp53$Odds),"Odds"] <- 50
  
  r_tp53$Color <- unlist(apply(r_tp53,1, function(x){
    if(x["state"] == "early"){
      if(as.numeric(x["Odds_err1"]) > 1){
        "#4DAF4A"
      }else if(as.numeric(x["Odds_err2"]) < 1){
        "#984EA3"
      }else{
        "#A5A4A4"
      }
    } else{
      if(as.numeric(x["Odds_err1"]) > 1){
        "#377EB8"
      }else if(as.numeric(x["Odds_err2"]) < 1){
        "#E41A1C"
      }else{
        "#A5A4A4"
      }
    }
  }))
  r_tp53$Color <- factor(r_tp53$Color)
  
  r_tp53$Odds_err1[r_tp53$Odds_err1 < 0.1] <- 0.1
  r_tp53$Odds_err2[r_tp53$Odds_err2 > 50] <- 50
  
  cbPalette <- c("#E41A1C", "#377EB8", "#984EA3", "#4DAF4A")
  
  q5 <- ggplot(r_tp53) + theme(panel.background = element_blank(), text = element_text(family = "mono"), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),legend.position="none")
  q5 <- q5 + aes(x = Gene, y = Odds, ymax = Odds_err2, ymin = Odds_err1, group = factor(state), color = Color) + ggtitle(gene)
  q5 <- q5 + geom_linerange(position = position_dodge(width = 0.5), color = "#A5A4A4") + geom_point(position = position_dodge(width = 0.5), size = 3) + geom_hline(yintercept = 1)
  q5 <- q5 + scale_y_log10(breaks = c(seq(0.1, 0.9, 0.1), seq(1, 9, 1), seq(10, 100, 10)), limits = c(0.09,60), labels = c(0.1,rep("",8),1,rep("",8),10,rep("",3),50,rep("",4),100), expand = c(0,0))
  q5 <- q5 + scale_color_manual(breaks = levels(r_tp53$Color), values = levels(r_tp53$Color))
  q5
  
  #ggsave("2018_08_tp53_odds.pdf", plot = q5, width = 7.5, height = 4.5)
  return(list(data = r_tp53, plot = q5))
}
