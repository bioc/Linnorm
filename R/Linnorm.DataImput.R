#' Linnorm Data Imputation Function.
#'
#' This function performs data imputation for (sc)RNA-seq expression data or large scale count data. It will treat every zero count in the dataset as missing data and replace them with predicted values.
#' @param datamatrix	The matrix or data frame that contains your dataset. It is only compatible with log transformed datasets.
#' @param DIMZP	Double >=0, <= 1. Genes not satisfying this threshold will be removed. For exmaple, if set to 0.3, genes without at least 30 percent of the samples being non-zero will be removed. Defaults to 0.5.
#' @param method	Character. Method for calculating the distance matrix. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "correlation", "spearman" or "kendall". Any unambiguous substring can be given. Defaults to "euclidean".
#' @details  This function performs data imputation on the dataset. It first separates the genes into numPortion portions. Then, for each portion, it performs linear regression between the non zero values of each sample and its corresponding weighted mean experssion to predict missing values. The linear regerssion lines are continuous throughout all genes. If the number of genes of each portion is less than 100, the missing values will simply be replaced by weighted means.
#' @return This function returns a data matrix.
#' @keywords Linnorm single cell RNA-seq data imputation missing value
#' @export
#' @examples
#' #Obtain example matrix:
#' data(Islam2011)
#' #Transformation:
#' transformedExp <- Linnorm.DataImput(Islam2011)
#' @import
#' Rcpp
#' RcppArmadillo

Linnorm.DataImput <- function(datamatrix, DIMZP=0.5, method="euclidean") {
	if (DIMZP > 1 || DIMZP < 0) {
		stop("Invalid DIMZP.")
	}
	datamatrix <- as.matrix(datamatrix)
	DM <- as.matrix(Dist(t(datamatrix),method = method))
	DM <- 1/(DM/rowSums(DM))
	DM[is.infinite(DM)] <- 0
	DM <- DM/rowSums(DM, na.rm=TRUE)
	datamatrix <- datamatrix[rowSums(datamatrix != 0) >= ncol(datamatrix) * DIMZP,]
	Answer <- datamatrix
	for (i in 1:ncol(datamatrix)) {
		zeroes <- which(datamatrix[,i] == 0)
		yvec <- WNZrowMeans(datamatrix[zeroes,],DM[i,])
		Answer[zeroes,i] <- yvec
	}
	return (Answer)
}
