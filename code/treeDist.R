#' ---
#' title: "Tufte Handout"
#' author: "Moritz Gerstung"
#' date: "August 13th, 2014"
#' output: rmarkdown::tufte_handout
#' ---
#' Tree distance
#' =============
#+ libs
library(igraph)
library(CoxHD)

#+ function
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

#' # Read data
#+ data, results='asis'
data <- read.table("BRCA.txt", header=TRUE, sep="\t")
kable(data)

#' Split
dataList <- split(data, data$Sample)

#' # Create tree structurs
#+ trees, fig.width=4, fig.height=4
allTrees <- lapply(dataList, function(x){	
			toTree(as.matrix(x[,c("Parent","Child")])+1, branch.length = x$Length, labels=as.character(x$Annotation))
})
names(allTrees) <- names(dataList)
i <- 1; for(tree in allTrees){
	plot(tree, layout=layout.reingold.tilford(tree), main=names(allTrees)[i], edge.label=E(tree)$name); i<-i+1}

#' # Compute distance of nodes
shortest.paths(allTrees$PD9694, mode="out") # example
allDist <- sapply(allTrees, function(tree){ #all
			e <- get.edgelist(tree)
			w <-  c(0,E(tree)$weight/2)
			d <- shortest.paths(tree, mode='out') - rep(w, each=length(w))
			diag(d) <- NA
			d
		})

#' # Get all mutations
u <- unlist(strsplit(as.character(data$Annotation), ";"))
table(u)
allMutations <- sort(unique(u))
allMutationsGermline <- c("germline",allMutations)

#' # Put into large array
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

#' # Fit model
#+ fit, results='asis'
observed <- !is.na(allDistMutations) & !is.infinite(allDistMutations) & c(TRUE,rep(FALSE, length(allMutations))) # only distance to root for the moment..
w <- which(observed, arr.ind = TRUE)
y <- allDistMutations[observed]
x <- MakeInteger(factor(w[,2], levels=seq_along(allMutationsGermline), labels=allMutationsGermline)) - MakeInteger(factor(w[,1], levels=seq_along(allMutationsGermline), labels=allMutationsGermline))

fit <- lm(y ~ . -1 -germline, data=x , weights=1/mapply(function(i,j) weights[i,j], w[,2]-1,w[,3]))
s <- summary(fit)
kable(s$coefficients)

#' # Plot
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
