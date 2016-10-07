#' Linnorm-hierarchical clustering analysis.
#'
#' This function first performs Linnorm transformation on the dataset. Then, it will perform hierarchical clustering analysis.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets. If a Linnorm transfored dataset is being used, please set the "input" argument into "Linnorm".
#' @param input	Character. "Raw" or "Linnorm". In case you have already transformed your dataset with Linnorm, set input into "Linnorm" so that you can input the Linnorm transformed dataset into the "datamatrix" argument. Defaults to "Raw".
#' @param showinfo	Logical. Show information about the computing process. Defaults to FALSE.
#' @param minZeroPortion	Double >=0, <= 1. For example, setting minZeroPortion as 0.5 will remove genes with more than half data values being zero in the calculation of normalizing parameter. Defaults to 0.
#' @param keepAll	Logical. After applying minZeroPortion filtering, should Linnorm keep all genes in the results? Defualts to TRUE.
#' @param perturbation	Integer >=2. To search for an optimal minimal deviation parameter (please see the article), Linnorm uses the iterated local search algorithm which perturbs away from the initial local minimum. The range of the area searched in each perturbation is exponentially increased as the area get further away from the initial local minimum, which is determined by their index. This range is calculated by 10 * (perturbation ^ index).
#' @param method_hclust	Charcter. Method to be used in hierarchical clustering. (From hclust {fastcluster}: the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".) Defaults to "ward.D2".
#' @param method_dist	Charcter. Method to be used in hierarchical clustering. (From Dist {amap}: the distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "correlation", "spearman" or "kendall". Any unambiguous substring can be given.) Defaults to "pearson".
#' @param  Group	Character vector with length equals to sample size. Each character in this vector corresponds to each of the columns (samples) in the datamatrix. If this is provided, sample names will be colored according to their group. Defaults to NULL.
#' @param num_Clust	Integer >= 0. Number of clusters in hierarchical clustering. No cluster will be highlighted if this is set to 0. Defaults to 4.
#' @param ClustRect	Logical. If num_Clust > 0, should a rectangle be used to highlight the clusters? Defaults to TRUE.
#' @param RectColor	Character. If ClustRect is TRUE, this controls the color of the rectangle. Defaults to "red".
#' @param fontsize	Numeric. Font size of the texts in the figure. Defualts to 0.5.
#' @param linethickness	Numeric. Controls the thickness of the lines in the figure. Defaults to 0.5.
#' @details  This function performs PCA clustering using Linnorm transformation.
#' @return It returns a list with the following objects:
##' \itemize{
##'  \item{Results:}{ If num_Clust > 0, this outputs a named vector that contains the cluster assignment information of each sample. Else, this outputs a number 0.}
##'  \item{plot:}{ Plot of hierarchical clustering.}
##'  \item{Linnorm:}{ Linnorm transformed and filtered data matrix.}
##' }
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric hierarchical Clustering
#' @export
#' @examples
#' #Obtain example matrix:
#' data(Islam2011)
#' #Example:
#' HClust.results <- Linnorm.HClust(Islam2011, Group=c(rep("ESC",48), rep("EF",44), rep("NegCtrl",4)), num_Clust=3, fontsize=2)

Linnorm.HClust <- function(datamatrix, showinfo = FALSE, input="Raw", perturbation=10, minZeroPortion=0, keepAll=TRUE, method_hclust="ward.D2", method_dist="pearson", Group=NULL, num_Clust=4, ClustRect=TRUE, RectColor="red", fontsize=0.5, linethickness=0.5) {
	if (input != "Raw" && input != "Linnorm") {
		stop("input argument is not recognized.")
	}
	if (length(Group) != 0) {
		if (length(Group) != length(datamatrix[1,])) {
			stop("Group must be a vector with the same length as sample size.")
		}
	}
	if (num_Clust < 0) {
		stop("Invalid number of clusters.")
	}
	
	expdata <- 0
	if (input == "Raw") {
		#Linnorm transformation
		expdata <- Linnorm(datamatrix, showinfo = showinfo, method="default",perturbation=perturbation, minZeroPortion = minZeroPortion, keepAll = keepAll)
	} 
	if (input == "Linnorm"){
		if (!keepAll) {
			datamatrix <- datamatrix[rowSums(datamatrix != 0) >= ncol(datamatrix) * minZeroPortion,]
		}
		expdata <- datamatrix
	}
	expdata <- t(expdata)
	
	#Clustering
	hc <- hclust(Dist(expdata, method = method_dist), ,method = method_hclust)
	dendr <- dendro_data(hc, type = "rectangle")
	
	#plot object
	render_plot <- 0
	
	#Render Color of plot
	if (length(unique(Group)) > num_Clust) {
		colorCode <- c("grey60", rainbow(length(unique(Group))))
	} else {
		colorCode <- c("grey60", rainbow(num_Clust))
	}
	
	#Cluster object
	clust <- 0
	
	if (num_Clust > 0) {
		clust <- cutree(hc, k = num_Clust)
		# Split dendrogram into upper grey section and lower coloured section
		height <- unique(dendr$segments$y)[order(unique(dendr$segments$y), decreasing = TRUE)]
		cut.height <- mean(c(height[num_Clust], height[num_Clust-1]))
		dendr$segments$line <- ifelse(dendr$segments$y == dendr$segments$yend &
		   dendr$segments$y > cut.height, 1, 2)
		dendr$segments$line <- ifelse(dendr$segments$yend  > cut.height, 1, dendr$segments$line)

		# Number the clusters
		dendr$segments$cluster <- c(-1, diff(dendr$segments$line))
		change <- which(dendr$segments$cluster == 1)
		for (i in 1:num_Clust) dendr$segments$cluster[change[i]] = i + 1
		dendr$segments$cluster <-  ifelse(dendr$segments$line == 1, 1, 
					 ifelse(dendr$segments$cluster == 0, NA, dendr$segments$cluster))
		dendr$segments$cluster <- na.locf(dendr$segments$cluster)
		
		#plotting
		render_plot <- ggplot() + 
		geom_segment(data = segment(dendr),
			aes(x=x, y=y, xend=xend, yend=yend, size=factor(line), colour=factor(cluster)), size=linethickness, lineend = "square", show.legend = FALSE
		) +
		scale_colour_manual(values = colorCode) +
		scale_y_reverse(expand = c(0.2, 0)) + 
		labs(x = NULL, y = NULL) +
		coord_flip() +
		theme(axis.line.y = element_blank(),
			axis.ticks.y = element_blank(),
			axis.text.y = element_blank(),
			axis.title.y = element_blank(),
			panel.background = element_rect(fill = "white"),
			panel.grid = element_blank()
		)
	} else {
		#Plotting without clusters
		render_plot <- ggplot() + 
		geom_segment(data = segment(dendr),
			aes(x=x, y=y, xend=xend, yend=yend, size=factor(line)), size=linethickness, lineend = "square", show.legend = FALSE
		) +
		scale_y_reverse(expand = c(0.2, 0)) + 
		labs(x = NULL, y = NULL) +
		coord_flip() +
		theme(axis.line.y = element_blank(),
			axis.ticks.y = element_blank(),
			axis.text.y = element_blank(),
			axis.title.y = element_blank(),
			panel.background = element_rect(fill = "white"),
			panel.grid = element_blank()
		)
	}
	
	if (length(Group) != 0) {
		#Set label color
		labcol <- as.numeric(label(dendr)[,3])
		unilab <- as.character(unique(as.character(Group)))
		for (i in 1:length(unilab)) {
			labcol[which(as.character(dendr$label$label) %in% rownames(expdata)[which(Group == unilab[i])])] <- i+1
		}
		render_plot <- render_plot + geom_text(data = label(dendr), 
			aes(x, y, label = label, colour=factor(labcol)),
			hjust = -0.2, size = fontsize, show.legend = FALSE
		)
	} else {
		render_plot <- render_plot + geom_text(data = label(dendr), 
			aes(x, y, label = label),
			hjust = -0.2, size = fontsize, show.legend = FALSE
		)
	}
	
	if (ClustRect) {
		if (length(rownames(expdata)) > length(unique(rownames(expdata)))) {
			warning("Duplicate sample names found. Rectangle not drawn.")
		} else {
			#rectangle
			clust.df <- data.frame(label=rownames(expdata), cluster=factor(clust))
			dendr2 <- merge(dendr[["labels"]],clust.df, by="label")
			rect <- aggregate(x~cluster,dendr2,range)
			rect <- data.frame(rect$cluster,rect$x)
			ymax <- mean(hc$height[length(hc$height)-((num_Clust-2):(num_Clust-1))])
			render_plot <- render_plot + geom_rect(data=rect, 
			aes(xmin=X1-.5, xmax=X2+.5, ymin=0, ymax=ymax), 
			color=RectColor, fill=NA, size=linethickness)
		}
	}
	render_plot <- ggplot_build(render_plot)

	listing <- list(clust, render_plot, t(expdata))
	results <- setNames(listing, c("Results", "plot", "Linnorm"))
	return (results)
}
