#' Linnorm-gene correlation pipeline for gene correlation study.
#'
#' This function first performs Linnorm transformation on the dataset. Then, it will calculate correlation coefficient for all pairs of genes.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets. If a Linnorm transfored dataset is being used, please set the "input" argument into "Linnorm".
#' @param input	Character. "Raw" or "Linnorm". In case you have already transformed your dataset with Linnorm, set input into "Linnorm" so that you can input the Linnorm transformed dataset into the "datamatrix" argument. Defaults to "Raw".
#' @param method	Character. "pearson", "kendall" or "spearman". Method for the calculation of correlation coefficients. Defaults to "pearson"
#' @param showinfo	Logical. Show lambda value calculated. Defaults to FALSE.
#' @param minZeroPortion	Double >=0, <= 1. For example, setting minZeroPortion as 0.5 will remove genes with more than half data values being zero in the calculation of normalizing parameter. Since this test is based on variance, which requires more non-zero values, it is suggested to set it to a larger value. Defaults to 0.5.
#' @param keepAll	Boolean. After applying minZeroPortion filtering, should Linnorm keep all genes in the results? Defualts to FALSE.
#' @param perturbation	Integer >=2. To search for an optimal minimal deviation parameter (please see the article), Linnorm uses the iterated local search algorithm which perturbs away from the initial local minimum. The range of the area searched in each perturbation is exponentially increased as the area get further away from the initial local minimum, which is determined by their index. This range is calculated by 10 * (perturbation ^ index).
#' @param sig.q	Double >=0, <= 1. Only gene pairs with q values less than this threshold will be included in the results. Defaults to 0.2.
#' @details  This function performed gene correlated study in the dataset by using Linnorm transformation.
#' @return This function will output a data frame with the following columns:
##' \itemize{
##'  \item{Gene1:}{ Name of gene 1.}
##'  \item{Gene2:}{ Name of gene 2.}
##'  \item{XPM1:}{ Gene 1 average expression level in XPM. If input is raw coutns or CPM, this column is in CPM unit. If input is RPKM, FPKM or TPM, this column is in the TPM unit.}
##'  \item{XPM2:}{ Gene 2 average expression level in XPM. If input is raw coutns or CPM, this column is in CPM unit. If input is RPKM, FPKM or TPM, this column is in the TPM unit.}
##'  \item{Cor:}{ Correlation coefficient between the two genes.}
##'  \item{p.value:}{ p value of the correlation coefficient.}
##'  \item{q.value:}{ q value of the correlation coefficient.}
##' }
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric correlation coefficient kendall pearson spearman
#' @export
#' @examples
#' data(Islam2011)
#' #Correlation of randomly chosen 100 genes
#' results <- Linnorm.Cor(Islam2011[sample(1:nrow(Islam2011),100),])

Linnorm.Cor <- function(datamatrix, input="Raw", method = "pearson", showinfo = FALSE, perturbation=10, minZeroPortion=0.5, keepAll=FALSE, sig.q=0.2) {
	if (input != "Raw" && input != "Linnorm") {
		stop("input argument is not recognized.")
	}
	if (method != "pearson" && method != "spearman"&& method != "kendall") {
		stop("method is not recognized.")
	}
	if (sig.q < 0 || sig.q > 1) {
		stop("Invalid sig.q value.")
	}
	
	expdata <- 0
	if (input == "Raw") {
		#Linnorm transformation
		expdata <- Linnorm(datamatrix, showinfo = showinfo, method="internal",perturbation=perturbation, minZeroPortion = minZeroPortion, keepAll = keepAll)
		datamatrix <- expdata[[1]]
		expdata <- log1p(datamatrix * expdata[[2]])
	} 
	if (input == "Linnorm"){
		expdata <- datamatrix
		datamatrix <- exp(datamatrix)
		for (i in seq_along(datamatrix[1,])) {
			datamatrix[,i] <- (datamatrix[,i])/sum(datamatrix[,i])
		}
	}

	correlation <- cor(t(expdata), method=method)
	datamatrix <- datamatrix[rownames(correlation),]
	XPM <- rowMeans(datamatrix * 1000000)
	correlations <- correlation[upper.tri(correlation,diag=FALSE)]
	#Index for locating genes with index in "correlations"
	index <- createUpperIndex(ncol(correlation), length(correlations))
	
	pvalues <- r.sig(correlations, ncol(expdata))
	qvalues <- p.adjust(pvalues,"BH")
	
	wanted <- which(qvalues <= sig.q)
	
	return (data.frame(Gene1=rownames(correlation)[index[wanted,1]],Gene2=rownames(correlation)[index[wanted,2]],XPM1=XPM[index[wanted,1]],XPM2=XPM[index[wanted,2]],Cor=correlations[wanted],p.value=pvalues[wanted],q.value=qvalues[wanted]))
}


