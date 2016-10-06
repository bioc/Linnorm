#' Linnorm-gene correlation network analysis.
#'
#' This function first performs Linnorm transformation on the dataset. Then, it will perform correlation network analysis on the dataset.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets. If a Linnorm transfored dataset is being used, please set the "input" argument into "Linnorm".
#' @param input	Character. "Raw" or "Linnorm". In case you have already transformed your dataset with Linnorm, set input into "Linnorm" so that you can input the Linnorm transformed dataset into the "datamatrix" argument. Defaults to "Raw".
#' @param method	Character. "pearson", "kendall" or "spearman". Method for the calculation of correlation coefficients. Defaults to "pearson"
#' @param showinfo	Logical. Show lambda value calculated. Defaults to FALSE.
#' @param minZeroPortion	Double >=0, <= 1. For example, setting minZeroPortion as 0.5 will remove genes with more than half data values being zero in the calculation of normalizing parameter. Since this test is based on correlation coefficient, which requires more non-zero values, it is suggested to set it to a larger value. Defaults to 2/3.
#' @param perturbation	Integer >=2. To search for an optimal minimal deviation parameter (please see the article), Linnorm uses the iterated local search algorithm which perturbs away from the initial local minimum. The range of the area searched in each perturbation is exponentially increased as the area get further away from the initial local minimum, which is determined by their index. This range is calculated by 10 * (perturbation ^ index).
#' @param sig.q	Double >=0, <= 1. Only gene pairs with q values less than this threshold will be included in the "Results" data frame. Defaults to 0.05.
#' @param plotNetwork	Logical. Should the program output the network plot to a file? An "igraph" object will be included in the output regardless. Defaults to TRUE. 
#' @param plotNumPairs	Integer >= 50. Number of gene pairs to be used in the network plot. Defaults to 5000.
#' @param plotdegree	Integer >= 0. In the network plot, genes (vertices) without at least this number of degree will be removed. Defaults to 0.
#' @param plotname	Character. Name of the network plot. File extension will be appended to it. Defaults to "networkplot".
#' @param plotformat	Character. "pdf" or "png". Network plot output format. Defaults to "png".
#' @param plotVertexSize	Double >0. Controls vertex Size in the network plot. Defaults to 1.
#' @param plotFontSize	Double >0. Controls font Size in the network plot. Defaults to 1.
#' @param plot.Pos.cor.col	Character. Color of the edges of positively correlated gene pairs. Defaults to "red".
#' @param plot.Neg.cor.col	Character. Color of the edges of negatively correlated gene pairs. Defaults to "green".
#' @param vertex.col	Character. "cluster" or a color. This controls the color of the vertices. Defaults to "cluster".
#' @param plotlayout	Character. "kk" or "fr". "kk" uses Kamada-Kawai algorithm in igraph to assign vertex and edges. It scales edge length with correlation strength. However, it can cause overlaps between vertices. "fr" uses Fruchterman-Reingold algorithm in igraph to assign vertex and edges. It prevents overlatps between vertices better than "kk", but edge lengths are not scaled to correlation strength. Defaults to "kk".
#' @param clusterMethod	Character. "cluster_edge_betweenness", "cluster_fast_greedy", "cluster_infomap", "cluster_label_prop", "cluster_leading_eigen", "cluster_louvain", "cluster_optimal", "cluster_spinglass" or "cluster_walktrap". These are clustering functions from the igraph package. Defaults to "cluster_edge_betweenness".
#' @details  This function performed gene correlated study in the dataset by using Linnorm transformation.
#' @return This function will output a list with the following objects:
##' \itemize{
##'  \item{Results:}{ A data frame containing the results of the analysis, showing only the significant results determined by "sig.q" (see below).}
##'  \item{Cor.Matrix:}{ The resulting correlation matrix between each gene. }
##'  \item{q.Matrix:}{ A matrix of q values of each of the correlation coefficient from Cor.Matrix. }
##'  \item{Cluster:}{ A data frame that shows which gene belongs to which cluster.}
##'  \item{igraph:}{ The igraph object for users who want to draw the network plot manually. }
##'  \item{Linnorm:}{ Linnorm transformed and filtered data matrix.}
##' }
#' @return The "Results" data frame has the following columns:
##' \itemize{
##'  \item{Gene1:}{ Name of gene 1.}
##'  \item{Gene2:}{ Name of gene 2.}
##'  \item{XPM1:}{ Gene 1 average expression level in XPM. If input is raw counts or CPM, this column is in CPM unit. If input is RPKM, FPKM or TPM, this column is in the TPM unit.}
##'  \item{XPM2:}{ Gene 2 average expression level in XPM. If input is raw counts or CPM, this column is in CPM unit. If input is RPKM, FPKM or TPM, this column is in the TPM unit.}
##'  \item{Cor:}{ Correlation coefficient between the two genes.}
##'  \item{p.value:}{ p value of the correlation coefficient.}
##'  \item{q.value:}{ q value of the correlation coefficient.}
##' }
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric correlation coefficient kendall pearson spearman
#' @export
#' @examples
#' data(Islam2011)
#' #Analysis on Islam2011 embryonic stem cells
#' results <- Linnorm.Cor(Islam2011[,1:48])

Linnorm.Cor <- function(datamatrix, input="Raw", method = "pearson", showinfo = FALSE, perturbation=10, minZeroPortion=2/3, sig.q=0.05, plotNetwork=TRUE, plotNumPairs=5000, plotdegree=0, plotname="networkplot", plotformat = "png", plotVertexSize=1, plotFontSize=1, plot.Pos.cor.col="red", plot.Neg.cor.col="green", vertex.col="cluster", plotlayout="kk", clusterMethod = "cluster_edge_betweenness") {
	keepAll <- FALSE
	if (input != "Raw" && input != "Linnorm") {
		stop("input argument is not recognized.")
	}
	if (method != "pearson" && method != "spearman"&& method != "kendall") {
		stop("method is not recognized.")
	}
	if (sig.q < 0 || sig.q > 1) {
		stop("Invalid sig.q value.")
	}
	if (plotNumPairs <50) {
		stop("plotNumPairs is too small.")
	}
	if (plotdegree < 0) {
		stop("Invalid plotdegree value.")
	}
	if (plotformat != "pdf" && plotformat != "png") {
		stop("plotformat is not recognized.")
	}
	if (plotVertexSize <= 0) {
		stop("Invalid plotVertexSize.")
	}
	if (plotFontSize <= 0) {
		stop("Invalid plotFontSize.")
	}
	if (plotlayout != "kk" && plotlayout != "fr") {
		stop("Invalid plotlayout.")
	}
	if (!areColors(vertex.col) && vertex.col != "cluster") {
		stop("Invalid vertex.col.")
	}
	igraphFun <- c("cluster_edge_betweenness", "cluster_fast_greedy", "cluster_infomap", "cluster_label_prop", "cluster_leading_eigen", "cluster_louvain", "cluster_optimal", "cluster_spinglass", "cluster_walktrap")
	if (!(clusterMethod %in% igraphFun)) {
		stop("Invalid clusterMethod.")
	}
	
	filtered <- 0
	expdata <- 0
	if (input == "Raw") {
		#Linnorm transformation
		expdata <- Linnorm(datamatrix, showinfo = showinfo, method="internal",perturbation=perturbation, minZeroPortion = minZeroPortion, keepAll = keepAll)
		datamatrix <- expdata[[1]]
		expdata <- log1p(datamatrix * expdata[[2]])
	} 
	if (input == "Linnorm"){
		filtered <- expdata[rowSums(expdata != 0) < ncol(expdata) * minZeroPortion,]
		datamatrix <- datamatrix[rowSums(datamatrix != 0) >= ncol(datamatrix) * minZeroPortion,]
		expdata <- datamatrix
		datamatrix <- exp(datamatrix) - 1
		for (i in seq_along(datamatrix[1,])) {
			datamatrix[,i] <- (datamatrix[,i])/sum(datamatrix[,i])
		}
	}
	correlation <- cor(t(expdata), method=method)
	datamatrix <- datamatrix[rownames(correlation),]
	XPM <- rowMeans(datamatrix * 1000000)
	correlations <- correlation[upper.tri(correlation,diag=FALSE)]
	
	#Index for locating genes with index in "correlations"
	#Note that if you reverse column 1 and column 2 in index, it becomes lower index.
	index <- createUpperIndex(ncol(correlation), length(correlations))
	
	pvalues <- r.sig(correlations, ncol(expdata))
	qvalues <- p.adjust(pvalues,"BH")
	
	qvaluematrix <- UpperToMatrix(qvalues,index)
	
	#Result matrix
	wanted <- which(qvalues <= sig.q)
	resultmatrix <- data.frame(Gene1=rownames(correlation)[index[wanted,1]],Gene2=rownames(correlation)[index[wanted,2]],XPM1=XPM[index[wanted,1]],XPM2=XPM[index[wanted,2]],Cor=correlations[wanted],p.value=pvalues[wanted],q.value=qvalues[wanted])
	
	
	#Network for the top "plotNumPairs" significant positively correlated gene pairs.
	AllGene1 <- as.character(resultmatrix[,1])
	AllGene2 <- as.character(resultmatrix[,2])
	
	if (plotNumPairs > length(resultmatrix[,6])) {
		plotNumPairs <- length(resultmatrix[,6])
	}
	
	orderbyCor <- order(resultmatrix[,6],decreasing=FALSE)
	nodes <- vector(mode="character",length(plotNumPairs) * 2)
	index <- 1
	for (i in 1:plotNumPairs) {
		nodes[index] <- AllGene1[orderbyCor[i]]
		index <- index + 1
		nodes[index] <- AllGene2[orderbyCor[i]]
		index <- index + 1
	}
	
	g1 <- graph(nodes,directed=FALSE)
	Thislayout <- 0
	
	positive <- which(resultmatrix[orderbyCor[1:plotNumPairs],5] > 0)
	negative <- which(resultmatrix[orderbyCor[1:plotNumPairs],5] < 0)
	E(g1)[positive]$color <- "red"
	E(g1)[negative]$color <- "green"
	if (plotlayout == "kk") {
		E(g1)$weight <- (resultmatrix[orderbyCor[1:plotNumPairs],5] - 3)^2
		Thislayout <- layout_with_kk(g1)
	}
	if (plotlayout == "fr") {
		E(g1)$weight <- resultmatrix[orderbyCor[1:plotNumPairs],5] + 1
		Thislayout <- layout_with_fr(g1)
	}
	
	#Clustering
	clustering <- 0
	if (clusterMethod == "cluster_edge_betweenness") { clustering <- cluster_edge_betweenness(g1)}
	if (clusterMethod == "cluster_fast_greedy") { clustering <- cluster_fast_greedy(g1)}
	if (clusterMethod == "cluster_infomap") { clustering <- cluster_infomap(g1)}
	if (clusterMethod == "cluster_label_prop") { clustering <- cluster_label_prop(g1)}
	if (clusterMethod == "cluster_leading_eigen") { clustering <- cluster_leading_eigen(g1)}
	if (clusterMethod == "cluster_louvain") { clustering <- cluster_louvain(g1)}
	if (clusterMethod == "cluster_optimal") { clustering <- cluster_optimal(g1)}
	if (clusterMethod == "cluster_spinglass") { clustering <- cluster_spinglass(g1)}
	if (clusterMethod == "cluster_walktrap") { clustering <- cluster_walktrap(g1)}	
	
	if (vertex.col == "cluster") {
		vertex.col <- clustering$membership
	}
	
	#Cluster results
	Clust.res <- data.frame(Gene=clustering$names,Cluster=clustering$membership)
	
	
	if (plotNetwork) {
		if(plotformat == "pdf") {
			pdf(paste(plotname,".pdf",sep=""),width = 10, height = 10)
		}
		if (plotformat == "png"){
			png(paste(plotname,".png",sep=""),res=2000, width = 10, height = 10, units = 'in')
		}
		plot1 <- plot(g1, vertex.color=vertex.col, vertex.size=0.8 * plotVertexSize,vertex.frame.color="transparent", vertex.label.color="black",vertex.label.cex=0.03 * plotFontSize, edge.curved=0.1,edge.width=0.05 * plotVertexSize, layout=Thislayout, margin=0)
		print(plot1)
		dev.off()
	}
	if (keepAll) {
		expdata <- rbind(expdata, filtered)
	}
	listing <- list(resultmatrix,correlation,qvaluematrix,Clust.res,g1,expdata)
	result <- setNames(listing, c("Results", "Cor.Matrix", "q.Matrix", "Cluster", "igraph", "Linnorm"))
	return (result)
}
