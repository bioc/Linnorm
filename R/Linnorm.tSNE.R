#' Linnorm t-SNE Clustering pipeline for subpopulation Analysis
#'
#' This function first performs Linnorm transformation on the dataset. Then, it will perform t-distributed stochastic neighbor embedding (t-SNE) dimensionality reduction on the dataset and use k-means clustering to identify subpopulations of cells.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets.
#' @param input	Character. "Raw" or "Linnorm". In case you have already transformed your dataset with Linnorm, set input into "Linnorm" so that you can put the Linnorm transformed dataset into the "datamatrix" argument. Defaults to "Raw".
#' @param num_PC	Integer >= 2. Number of principal componenets to be used in K-means clustering. Defaults to 3.
#' @param  num_center	Numeric vector. Number of clusters to be tested for k-means clustering. fpc, vegan, mclust and apcluster packages are used to determine the number of clusters needed. If only one number is supplied, it will be used and this test will be skipped. Defaults to c(1:20).
#' @param  Group	Character vector with length equals to sample size. Each character in this vector corresponds to each of the columns (samples) in the datamatrix. In the plot, the shape of the points that represent each sample will be indicated by their group assignment. Defaults to NULL.
#' @param Coloring	Character. "kmeans" or "Group". If Group is not NA, coloring in the plot will reflect each sample's group. Otherwise, coloring will reflect k means clustering results. Defaults to "Group".
#' @param  kmeans.iter	Numeric. Number of iterations in k-means clustering. Defaults to 2000.
#' @param plot.title	Character. Set the title of the plot. Defaults to "t-SNE K-means clustering".
#' @param ... arguments that will be passed into Linnorm's transformation function.
#' @details  This function performs t-SNE K-means clustering using Linnorm transformation.
#' @return It returns a list with the following objects:
##' \itemize{
##'  \item{k_means:}{ Output of kmeans(for K-means clustering) from the stat package. Note: It contains a "cluster" object that indicates each sample's cluster assignment.}
##'  \item{tSNE:}{ Output from Rtsne.}
##'  \item{plot:}{ Plot of t-SNE K-means clustering.}
##'  \item{Linnorm:}{ Linnorm transformed and filtered data matrix.}
##' }
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric PCA Principal Component Analysis k-means K-means kmeans Clustering t-distributed stochastic neighbor embedding t-SNE
#' @export
#' @examples
#' #Obtain example matrix:
#' data(Islam2011)
#' #Example:
#' tSNE.results <- Linnorm.tSNE(Islam2011)
Linnorm.tSNE <- function(datamatrix, input = "Raw", num_PC=2, num_center=c(1:20), Group=NULL, Coloring="kmeans", kmeans.iter=2000, plot.title="t-SNE K-means clustering",...) {
	if (input != "Raw" && input != "Linnorm") {
		stop("input argument is not recognized.")
	}
	if (Coloring != "Group" && Coloring != "kmeans") {
		stop("Coloring argument is not recognized.")
	}
	if (num_PC < 2) {
		stop("num_PC is too small.")
	}
	if (!is.numeric(num_center)) {
		stop("num_center must be a vector of integers.")
	}
	if (length(Group) > 0) {
		if (length(Group) != length(datamatrix[1,])) {
			stop("Group must be a vector with the same length as sample size.")
		}
	}
	if (kmeans.iter < 10) {
		stop("kmeans.iter is too small.")
	}
	expdata <- 0
	if (input == "Raw") {
		#Linnorm transformation
		expdata <- Linnorm(datamatrix, ...)
	}
	x <- list(...)
	if (input == "Linnorm"){
		x <- list(...)
		if (sum(x$Filter == TRUE) == 1  && is.numeric(x$minZeroPortion)) {
			if (x$minZeroPortion > 1 || x$minZeroPortion < 0) {
				stop("Invalid minZeroPortion.")
			}
			datamatrix <- datamatrix[rowSums(datamatrix != 0) >= ncol(datamatrix) * x$minZeroPortion,]
		}
		expdata <- datamatrix
	}
	
	
	#Principal Component Analysis
	res.pca <- Rtsne(t(expdata),dims=num_PC)
	
	#Extract Principal Components for k means clustering.
	data <- res.pca$Y
	
	num_clust <- c()
	if (length(num_center) == 1) {
		num_clust <- num_center
	} else {
		#Automatically determine number of centers by fpc, vegan, mclust and apcluster packages.
		#fpc
		pamk.best <- pamk(data,num_center)

		num_clust <- c(pamk.best$nc)

		#vegan
		fit <- cascadeKM(scale(data, center = TRUE,  scale = TRUE), min(num_center), max(num_center), iter = 500)
		calinski.best <- as.numeric(which.max(fit$results[2,]))

		num_clust <- c(num_clust, calinski.best)
		
		#mclust
		d_clust <- Mclust(as.matrix(data), G=num_center)
		m.best <- dim(d_clust$z)[2]
		
		num_clust <- c(num_clust, m.best)
		
		#apcluster
		d.apclus <- apcluster(negDistMat(r=2), data)
		num_clust <- c(num_clust, length(d.apclus@clusters))
		
		#Mode of num_clust
		ux <- unique(num_clust)
		if (length(ux) == 4) {
			num_clust <- sort(ux)[2]
		} else {
			num_clust <- ux[which.max(tabulate(match(num_clust, ux)))]
		}
		
		
		if (sum(x$showinfo == TRUE) == 1) {
			cat("Number of clusters from the fpc package:", pamk.best$nc, "\n")
			cat("Number of clusters from the vegan package:", calinski.best, "\n")
			cat("Number of clusters from the mclust package:", m.best, "\n")
			cat("Number of clusters from the apcluster package:", length(d.apclus@clusters), "\n")
			cat("Final number of clusters:", num_clust, "\n")
		}
	}
	
	#K-means clustering	
	results <- kmeans(data, num_clust, iter.max = kmeans.iter)

	#Plotting
	PC1 <- c()
	PC2 <- c()
	Cluster <- c()
	x <- c()
	y <- c()
	
	for (i in 1:length(data[,1])) {
		PC1 <- c(PC1, data[i,1])
		PC2 <- c(PC2, data[i,2])
		Cluster <- c(Cluster, paste("Cluster",results[[1]][i] ))
	}
	
	
	
	if (length(Group) == 0) {
		plotdata <- data.frame(PC1=PC1,PC2=PC2,Cluster=Cluster)
		#Find cluster elipse
		df_ell <- data.frame()
		for(g in plotdata$Cluster){
			df_ell <- rbind(df_ell, cbind(as.data.frame(with(plotdata[plotdata$Cluster==g,], ellipse(cor(PC1, PC2),scale=c(sd(PC1),sd(PC2)),centre=c(mean(PC1),mean(PC2))))),Cluster=g))
		}

		render_plot <- ggplot_build(ggplot(plotdata, aes(x=PC1, y=PC2, color=Cluster)) + geom_point(aes(shape=Cluster), size = 2) + geom_path(data=df_ell, aes(x=x, y=y,colour=Cluster), size=0.5, linetype=2) + scale_x_continuous("PC1") + scale_y_continuous("PC2") + ggtitle(plot.title) + theme(aspect.ratio=1))
	} else {
		shaping <- c()
		numShape <- length(unique(Group))
		shapeindex <- 1
		for (i in 1:numShape) {
			if (shapeindex <= 25) {
				shaping <- c(shaping, shapeindex)
				shapeindex <- shapeindex + 1
			} else {
				shapeindex <- 1
				shaping <- c(shaping, shapeindex)
			}
		}
		if (Coloring == "kmeans") {
			plotdata <- data.frame(PC1=PC1,PC2=PC2,Cluster=Cluster,Group=Group)
			#Find cluster elipse
			df_ell <- data.frame()
			for(g in plotdata$Cluster){
				df_ell <- rbind(df_ell, cbind(as.data.frame(with(plotdata[plotdata$Cluster==g,], ellipse(cor(PC1, PC2),scale=c(sd(PC1),sd(PC2)),centre=c(mean(PC1),mean(PC2))))),Cluster=g))
			}

			render_plot <- ggplot_build(ggplot(plotdata, aes(x=PC1, y=PC2, color=Cluster)) + geom_point(aes(shape=Group), size = 2) + scale_shape_manual(values=shaping) + geom_path(data=df_ell, aes(x=x, y=y,colour=Cluster), size=0.5, linetype=2) + scale_x_continuous("PC1") + scale_y_continuous("PC2") + ggtitle(plot.title) + theme(aspect.ratio=1))
		} else if (Coloring == "Group") {
			plotdata <- data.frame(PC1=PC1,PC2=PC2,Cluster=Group,Group=Group)
			#Find cluster elipse
			df_ell <- data.frame()
			for(g in plotdata$Group){
				df_ell <- rbind(df_ell, cbind(as.data.frame(with(plotdata[plotdata$Group==g,], ellipse(cor(PC1, PC2),scale=c(sd(PC1),sd(PC2)),centre=c(mean(PC1),mean(PC2))))),Group=g))
			}

			render_plot <- ggplot_build(ggplot(plotdata, aes(x=PC1, y=PC2, color=Group)) + geom_point(aes(shape=Group), size = 2) + scale_shape_manual(values=shaping) + geom_path(data=df_ell, aes(x=x, y=y,colour=Group), size=0.5, linetype=2) + scale_x_continuous("PC1") + scale_y_continuous("PC2") + ggtitle(plot.title) + theme(aspect.ratio=1))
		}
	}
	listing <- list(results, res.pca, render_plot, expdata)
	results <- setNames(listing, c("k_means", "tSNE", "plot", "Linnorm"))
	return (results)
}
