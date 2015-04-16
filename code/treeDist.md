---
title: "Tufte Handout"
author: "Moritz Gerstung"
date: "August 13th, 2014"
output: rmarkdown::tufte_handout
---
Tree distance
=============


```r
library(igraph)
library(CoxHD)
```

```r
toTree <- function(edgelist, branch.length, labels){
	g <- graph.edgelist(edgelist)
	o <- order(get.edgelist(g)[,2])
	E(g)$weight[o] <- branch.length
	E(g)$name[o] <- labels
	return(g)	
}

toPhylo <- function(edge, edge.length, labels){
	edge <- as.matrix(edge)
	dimnames(edge) <- NULL
	all <- unique(c(edge[,1], edge[,2]))
	tips <- setdiff(edge[,2],edge[,1])
	nodes <- setdiff(all, tips)
	order <- c(tips, nodes)
	o <- order(order)
	l <- list(edge=cbind(o[edge[,1]],o[edge[,2]]), edge.length=edge.length, Nnode=length(nodes), tip.label=labels[tips], node.labels=labels[nodes])
	class(l) <- "phylo"
	attr(l, "order") <- "cladewise"
	return(l)
}

#tree <- list(`1`=2,`2`=c(3,4),`3`=5:7)
```

# Read data


```r
data <- read.table("BRCA.txt", header=TRUE, sep="\t")
kable(data)
```



|Sample | Parent| Child|Annotation                                       | Length| Total|
|:------|------:|-----:|:------------------------------------------------|------:|-----:|
|PD9775 |      0|     1|TP53                                             |   0.95|  3104|
|PD9775 |      1|     2|                                                 |   0.05|    NA|
|PD9773 |      0|     1|TP53;ARID1B                                      |   0.95|  2036|
|PD9773 |      1|     2|                                                 |   0.05|    NA|
|PD9776 |      0|     1|PTEN;ARID1B;MYC;FGFR2;17p-;BRCA1                 |   0.95| 13219|
|PD9776 |      1|     2|MYC                                              |   0.05|    NA|
|PD9776 |      1|     3|                                                 |   0.05|    NA|
|PD9777 |      0|     1|TP53;CBL;JAK2;17p-                               |   0.80|  4851|
|PD9777 |      1|     2|                                                 |   0.05|    NA|
|PD9777 |      2|     3|FGFR2;MYC                                        |   0.05|    NA|
|PD9777 |      2|     4|FGFR2                                            |   0.05|    NA|
|PD9777 |      4|     5|                                                 |   0.05|    NA|
|PD9769 |      0|     1|TP53;CTNNB1;PTEN;MYC;CCNE1;17p-                  |   0.70|  6582|
|PD9769 |      1|     2|RUNX1                                            |   0.20|    NA|
|PD9769 |      2|     3|                                                 |   0.10|    NA|
|PD9769 |      1|     4|IGF2R;RUNX1                                      |   0.20|    NA|
|PD9771 |      0|     1|TP53;FGFR1;WHSC1L1;ERLIN2;IGF1R;CDK4;AR;MET;17p- |   0.70| 11479|
|PD9771 |      1|     2|NRAS;CDKN2B                                      |   0.05|    NA|
|PD9771 |      1|     3|CIC                                              |   0.10|    NA|
|PD9771 |      3|     4|RAD21;NCOR                                       |   0.20|    NA|
|PD9771 |      3|     5|                                                 |   0.05|    NA|
|PD9771 |      5|     6|                                                 |   0.05|    NA|
|PD9770 |      0|     1|TP53;PIK3CA;EGFR                                 |   0.40|  7624|
|PD9770 |      1|     2|EGFR;NRG1;CDC7                                   |   0.45|    NA|
|PD9770 |      1|     3|CDK6;CUX1;3p-                                    |   0.40|    NA|
|PD9770 |      3|     4|CUX1                                             |   0.20|    NA|
|PD9770 |      3|     5|ERBB2                                            |   0.20|    NA|
|PD9849 |      0|     1|TP53;PTEN;FGFR2;ERLIN2;WHSC1L1;FGFR1;17p-        |   0.90|  4212|
|PD9849 |      1|     2|CBFB                                             |   0.05|    NA|
|PD9849 |      1|     3|                                                 |   0.10|    NA|
|PD9849 |      1|     4|                                                 |   0.05|    NA|
|PD9852 |      0|     1|TP53;PIK3CA;PBRM1;NF1;EGFR;ERBB2;17p-            |   0.85|  6231|
|PD9852 |      1|     2|                                                 |   0.10|    NA|
|PD9852 |      1|     3|PBRM1                                            |   0.15|    NA|
|PD9694 |      0|     1|SF3B1;CREBBP;PTEN                                |   0.20|  3371|
|PD9694 |      1|     2|BMX;FGFR1;AKT3                                   |   0.40|    NA|
|PD9694 |      2|     3|PTEN                                             |   0.40|    NA|
|PD9694 |      2|     4|                                                 |   0.10|    NA|
|PD9694 |      1|     5|TERT                                             |   0.20|    NA|
|PD9694 |      5|     6|                                                 |   0.60|    NA|

Split


```r
dataList <- split(data, data$Sample)
```

# Create tree structurs


```r
allTrees <- lapply(dataList, function(x){	
			toTree(as.matrix(x[,c("Parent","Child")])+1, branch.length = x$Length, labels=as.character(x$Annotation))
})
names(allTrees) <- names(dataList)
i <- 1; for(tree in allTrees){
	plot(tree, layout=layout.reingold.tilford(tree), main=names(allTrees)[i], edge.label=E(tree)$name); i<-i+1}
```

![plot of chunk trees](figure/trees1.png) ![plot of chunk trees](figure/trees2.png) ![plot of chunk trees](figure/trees3.png) ![plot of chunk trees](figure/trees4.png) ![plot of chunk trees](figure/trees5.png) ![plot of chunk trees](figure/trees6.png) ![plot of chunk trees](figure/trees7.png) ![plot of chunk trees](figure/trees8.png) ![plot of chunk trees](figure/trees9.png) ![plot of chunk trees](figure/trees10.png) 

# Compute distance of nodes


```r
shortest.paths(allTrees$PD9694, mode="out") # example
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
## [1,]    0  0.2  0.6  1.0  0.7  0.4  1.0
## [2,]  Inf  0.0  0.4  0.8  0.5  0.2  0.8
## [3,]  Inf  Inf  0.0  0.4  0.1  Inf  Inf
## [4,]  Inf  Inf  Inf  0.0  Inf  Inf  Inf
## [5,]  Inf  Inf  Inf  Inf  0.0  Inf  Inf
## [6,]  Inf  Inf  Inf  Inf  Inf  0.0  0.6
## [7,]  Inf  Inf  Inf  Inf  Inf  Inf  0.0
```

```r
allDist <- sapply(allTrees, function(tree){ #all
			e <- get.edgelist(tree)
			w <-  c(0,E(tree)$weight/2)
			d <- shortest.paths(tree, mode='out') - rep(w, each=length(w))
			diag(d) <- NA
			d
		})
```

# Get all mutations


```r
u <- unlist(strsplit(as.character(data$Annotation), ";"))
table(u)
```

```
## u
##    17p-     3p-    AKT3      AR  ARID1B     BMX   BRCA1    CBFB     CBL 
##       6       1       1       1       2       1       1       1       1 
##   CCNE1    CDC7    CDK4    CDK6  CDKN2B     CIC  CREBBP  CTNNB1    CUX1 
##       1       1       1       1       1       1       1       1       2 
##    EGFR   ERBB2  ERLIN2   FGFR1   FGFR2   IGF1R   IGF2R    JAK2     MET 
##       3       2       2       3       4       1       1       1       1 
##     MYC    NCOR     NF1    NRAS    NRG1   PBRM1  PIK3CA    PTEN   RAD21 
##       4       1       1       1       1       2       2       5       1 
##   RUNX1   SF3B1    TERT    TP53 WHSC1L1 
##       2       1       1       8       2
```

```r
allMutations <- sort(unique(u))
allMutationsGermline <- c("germline",allMutations)
```

# Put into large array


```r
allDistMutations <- array(0, c(rep(length(allMutationsGermline),2), length(dataList)), dimnames=list(allMutationsGermline, allMutationsGermline, names(dataList)))
for(n in names(dataList)){
	m <- c(germline=1,sapply(allMutations, function(x) {g <- grep(x, dataList[[n]]$Annotation); if(length(g)==0) NA else g[1]+1}))
	d <- allDist[[n]][m,m]
	#d[is.infinite(d)] <- NA
	allDistMutations[,,n] <- d
}

weights <- sapply(dataList, function(x){
				m <- sapply(allMutations, function(y) {g <- grep(y, x$Annotation); if(length(g)==0) NA else g[1]})
				1/x$Length[m]
		})
rownames(weights) <- allMutations
```

# Fit model


```r
observed <- !is.na(allDistMutations) & !is.infinite(allDistMutations) & c(TRUE,rep(FALSE, length(allMutations))) # only distance to root for the moment..
w <- which(observed, arr.ind = TRUE)
y <- allDistMutations[observed]
x <- MakeInteger(factor(w[,2], levels=seq_along(allMutationsGermline), labels=allMutationsGermline)) - MakeInteger(factor(w[,1], levels=seq_along(allMutationsGermline), labels=allMutationsGermline))

fit <- lm(y ~ . -1 -germline, data=x , weights=1/mapply(function(i,j) weights[i,j], w[,2]-1,w[,3]))
s <- summary(fit)
kable(s$coefficients)
```



|        | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|
|:-------|--------:|----------:|-------:|------------------:|
|`17p-`  |   0.4138|     0.0355| 11.6585|             0.0000|
|`3p-`   |   0.6000|     0.1242|  4.8301|             0.0000|
|AKT3    |   0.4000|     0.1242|  3.2201|             0.0032|
|AR      |   0.4413|     0.0487|  9.0582|             0.0000|
|ARID1B  |   0.4750|     0.0570|  8.3339|             0.0000|
|BMX     |   0.4000|     0.1242|  3.2201|             0.0032|
|BRCA1   |   0.4750|     0.0806|  5.8930|             0.0000|
|CBFB    |   0.9250|     0.3513|  2.6327|             0.0134|
|CBL     |   0.4000|     0.0878|  4.5539|             0.0001|
|CCNE1   |   0.3500|     0.0939|  3.7273|             0.0008|
|CDC7    |   0.6250|     0.1171|  5.3366|             0.0000|
|CDK4    |   0.3500|     0.0939|  3.7273|             0.0008|
|CDK6    |   0.6000|     0.1242|  4.8301|             0.0000|
|CDKN2B  |   0.7250|     0.3513|  2.0635|             0.0481|
|CIC     |   0.7500|     0.2484|  3.0188|             0.0052|
|CREBBP  |   0.1000|     0.1757|  0.5692|             0.5736|
|CTNNB1  |   0.3500|     0.0939|  3.7273|             0.0008|
|CUX1    |   0.6000|     0.1242|  4.8301|             0.0000|
|EGFR    |   0.3530|     0.0703|  5.0235|             0.0000|
|ERBB2   |   0.5155|     0.0767|  6.7233|             0.0000|
|ERLIN2  |   0.4062|     0.0621|  6.5408|             0.0000|
|FGFR1   |   0.4050|     0.0556|  7.2903|             0.0000|
|FGFR2   |   0.4737|     0.0570|  8.3108|             0.0000|
|IGF1R   |   0.3500|     0.0939|  3.7273|             0.0008|
|IGF2R   |   0.8000|     0.1757|  4.5539|             0.0001|
|JAK2    |   0.4000|     0.0878|  4.5539|             0.0001|
|MET     |   0.3500|     0.0939|  3.7273|             0.0008|
|MYC     |   0.4353|     0.0603|  7.2241|             0.0000|
|NCOR    |   0.9000|     0.1757|  5.1231|             0.0000|
|NF1     |   0.4250|     0.0852|  4.9874|             0.0000|
|NRAS    |   0.7250|     0.3513|  2.0635|             0.0481|
|NRG1    |   0.6250|     0.1171|  5.3366|             0.0000|
|PBRM1   |   0.4250|     0.0852|  4.9874|             0.0000|
|PIK3CA  |   0.3530|     0.0703|  5.0235|             0.0000|
|PTEN    |   0.4077|     0.0474|  8.6063|             0.0000|
|RAD21   |   0.9000|     0.1757|  5.1231|             0.0000|
|RUNX1   |   0.8000|     0.1757|  4.5539|             0.0001|
|SF3B1   |   0.1000|     0.1757|  0.5692|             0.5736|
|TERT    |   0.3000|     0.1757|  1.7077|             0.0984|
|TP53    |   0.4094|     0.0314| 13.0276|             0.0000|
|WHSC1L1 |   0.4062|     0.0621|  6.5408|             0.0000|

# Plot


```r
r <- rank(coef(fit), ties.method = 'random')
c <- coef(fit)
par(bty="n", xpd=NA, mar=c(5,3,1,5))
plot(c, r, xlim=c(0,1), xlab="relative time", yaxt='n', pch=19, cex=sqrt(table(u)[gsub("`","",names(r))]/2))
v1 <- s$coefficients[,2]^2
d1 <- colSums(abs(x[,-1]))
d0 <- 1 # prior df
v0 <- 1 # prior variance, i.e., any time
v <- (d1*v1 + d0 * v0)/(d0+d1)  
sd <- sqrt(v)
segments(pmax(0,c-sd),r, pmin(1,c+sd),r)
text(pmin(1,c+sd),r, paste0(names(c), ", n=",table(u)), pos=4)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 

