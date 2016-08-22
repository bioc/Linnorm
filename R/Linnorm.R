#' Linnorm Transformation Function
#'
#' This function performs the Linear model and normality based transformation method (Linnorm) for RNA-seq expression data or large scale count data.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets.
#' @param showinfo Logical. Show lambda value calculated. Defaults to FALSE.
#' @param method	"default" or "lambda" The program will output the transformed matrix if the method is "default". If the method is "lambda", the program will output a lambda value.
#' @param minZeroPortion Double >=0, <= 1. For example, setting minZeroPortion as 0.5 will remove genes with more than half data values being zero in the calculation of normalizing parameter. It is strongly suggested to change this to 0 for single cell RNA-seq data. Defaults to 2/3.
#' @param keepAll	Boolean. After applying minZeroPortion filtering, should Linnorm keep all genes in the results? Defualts to TRUE.
#' @param perturbation Integer >=2. To search for an optimal minimal deviation parameter (please see the article), Linnorm uses the iterated local search algorithm which perturbs away from the initial local minimum. The range of the area searched in each perturbation is exponentially increased as the area get further away from the initial local minimum, which is determined by their index. This range is calculated by 10 * (perturbation ^ index).
#' @details  If method is default, Linnorm outputs a transformed expression matrix. For users who wish to work with lambda instead, the output is a single lambda value. Please note that users with the lambda value can obtain a transformed Linnorm dataset by: log1p(lambda * datamatrix). There is no need to rerun the program if a lambda is already calculated.
#' @return This function returns either a transformed data matrix or a lambda value.
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric
#' @export
#' @examples
#' #Obtain example matrix:
#' library(seqc)
#' SampleA <- ILM_aceview_gene_BGI[,grepl("A_",colnames(ILM_aceview_gene_BGI))]
#' rownames(SampleA) <- ILM_aceview_gene_BGI[,2]
#' #Extract a portion of the matrix for an example
#' expMatrix <- SampleA[,1:3]
#' transformedExp <- Linnorm(expMatrix)
#' transformedExp <- Linnorm(expMatrix, method = "lambda")
#' @import
#' Rcpp
#' RcppArmadillo
Linnorm <- function(datamatrix, showinfo = FALSE, method="default",perturbation=10, minZeroPortion=2/3, keepAll = TRUE) {
	#data checking
	expdata <- as.matrix(datamatrix)
	if (length(expdata[1,]) < 2) {
		stop("Number of samples is less than 2.")
	}
	if (length(expdata[,1]) < 10) {
		stop("Number of features is too small.")
	}
	if (perturbation < 2) {
		stop("perturbation is too small.")
	}
	if (method != "default" && method != "lambda" && method != "internal") {
		stop("Invalid algorithm value.")
	}
	if (minZeroPortion >1 || minZeroPortion < 0) {
		stop("Invalid minZeroPortion.")
	}
	if (anyNA(expdata)) {
		stop("Dataset contains NA.")
	}
	if (sum(which(expdata < 0)) != 0) {
		stop("Dataset contains negative number.")
	}
	#Step 1: Relative Expression
	#Turn it into relative expression
	for (i in seq_along(expdata[1,])) {
		expdata[,i] <- expdata[,i]/sum(expdata[,i])
	}
	#Save the original dataset before trimming.
	datamatrix <- expdata
	expdata <- expdata[order(rowMeans(expdata)),]
	#trim outliers
	if (minZeroPortion == 0) {
		expdata <- expdata[rowSums(expdata != 0) > 0,]
	} else {
		expdata <- expdata[rowSums(expdata != 0) >= (length(expdata[1,]) * minZeroPortion),]
	}
	if (method == "default") {
		
		X <- LocateLambda(expdata,perturbation)
		if (showinfo) {
			message("Lambda is ", X,".",appendLF=TRUE)
			flush.console()
		}
		if (keepAll) {
			return( log1p(datamatrix * X) )
		} else {
			return( log1p(expdata * X) )
		}
	}
	if (method == "lambda") {
		X <- LocateLambda(expdata,perturbation)
		return( X )
	}
	if (method == "internal") {
		X <- LocateLambda(expdata,perturbation)
		if (showinfo) {
			message("Lambda is ", X,".",appendLF=TRUE)
			flush.console()
		}
		if (keepAll) {
			expdata <- datamatrix
		}
		listing <- list(expdata,X)
		return (listing)
	}
}

