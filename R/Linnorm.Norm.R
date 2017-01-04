#' Linnorm Normalization Function
#'
#' This function performs batch effect and library size difference normalization on the input dataset.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets.
#' @param showinfo	Logical. Show algorithm running information. Defaults to FALSE.
#' @param output	character. "Raw" or "XPM". Output's total count will be approximately the median of the inputs' when set to "Raw". Output CPM (if input is raw counts or CPM) or TPM (if input is RPKM FPKM or TPM) when set to "XPM". 
#' @param minZeroPortion Double >=0, <= 1. Genes not satisfying this threshold will be removed. For exmaple, if set to 0.3, genes without at least 30 percent of the samples being non-zero will be removed. Defaults to 0.5.
#' @param BE_F_p	Double >=0, <= 1. Filter genes with standard deviation and skewness less than this p value before applying Linnorm's batch effect normalization algorithm. Defaults to 0.3173.
#' @param BE_F_LC_Genes	Double >= 0.01, <= 0.95 or Character "Auto". Filter this portion of the lowest expressing genes before applying Linnorm's batch effect normalization algorithm. It can be determined automatically by setting to "Auto". Defaults to "Auto".
#' @param BE_F_HC_Genes	Double >=0, <= 1. Filter this portion of the highest expressing genes before applying Linnorm's batch effect normalization algorithm. Defaults to 0.01.
#' @param BE_strength	Double >=0, <= 1. Before Linnorm transformation, how strongly should Linnorm normalize batch effects? Defaults to 0.5.
#' @param max_F_LC	Double >=0, <= 0.95. When L_F_LC or B_F_LC is set to auto, this is the maximum threshold that Linnorm would assign. Defaults to 0.75.
#' @details  This function normalizes the input dataset using the Linnorm algorithm.
#' @return This function returns a normalized data matrix.
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization
#' @export
#' @examples
#' #Obtain example matrix:
#' data(LIHC)
#' #Transformation:
#' transformedExp <- Linnorm(LIHC)
#' @import
#' Rcpp
#' RcppArmadillo
Linnorm.Norm <- function (datamatrix, showinfo=FALSE, output="XPM", minZeroPortion = 0.5, BE_F_p = 0.3173, BE_F_LC_Genes = "Auto", BE_F_HC_Genes = 0.01, BE_strength = 0.5, max_F_LC = 0.75) {
	#data checking
	datamatrix <- as.matrix(datamatrix)
	if (length(datamatrix[1,]) < 3) {
		stop("Number of samples is less than 3.")
	}
	if (length(datamatrix[,1]) < 500) {
		stop("Number of features is too small.")
	}
	if (output != "Raw" && output != "XPM") {
		stop("Invalid output argument. It must be Raw or XPM.")
	}
	if (minZeroPortion > 1 || minZeroPortion < 0) {
		stop("Invalid minZeroPortion.")
	}
	if (BE_F_p > 1 || BE_F_p < 0) {
		stop("Invalid BE_F_p.")
	}
	if (BE_strength > 1 || BE_strength < 0) {
		stop("Invalid BE_strength.")
	}
	if (BE_F_LC_Genes > 0.75 || BE_F_LC_Genes < 0.01) {
		if (BE_F_LC_Genes != "Auto") {
			stop("Invalid BE_F_LC_Genes.")
		}
	}
	if (BE_F_HC_Genes > 0.75 || BE_F_HC_Genes < 0.01) {
		stop("Invalid BE_F_HC_Genes.")
	}
	if (anyNA(datamatrix)) {
		stop("Dataset contains NA.")
	}
	if (sum(which(datamatrix < 0)) != 0) {
		stop("Dataset contains negative number.")
	}
	if (max_F_LC > 0.95 || max_F_LC < 0) {
		stop("Invalid max_F_LC.")
	}
	
	#Step 1: Relative Expression
	#Turn it into relative expression
	#Note that expdata does not have colnames and rownames now
	RN <- rownames(datamatrix)
	CN <- colnames(datamatrix)
	multy <- 0
	if (output == "Raw") {
		multy <- median(colSums(datamatrix))
	} else {
		multy <- 1000000
	}
	datamatrix <- XPM(datamatrix) * multy
	
	Keep <- 0
	if (minZeroPortion == 0) {
		Keep <- which(rowSums(datamatrix != 0) >= ncol(datamatrix) * minZeroPortion)
	} else {
		Keep <- which(rowSums(datamatrix != 0) > ncol(datamatrix) * minZeroPortion)
	}
	if (BE_F_LC_Genes == "Auto") {
		BE_F_LC_Genes <- FindLCT(datamatrix[Keep,], multy)
		if (BE_F_LC_Genes > max_F_LC) {
			if (showinfo) {
				message(paste("Filter low count gene threshold is ", LC_Threshold, ". It is larger than max_F_LC, ", max_F_LC, ", which is now used.", sep=""))
				BE_F_LC_Genes <- max_F_LC
			}
		}
		if (showinfo) {
			message(paste("Filter low count genes threshold is set to ", LC_Threshold, sep=""),appendLF=TRUE)
			flush.console()
		}
	}

	if (BE_F_LC_Genes + BE_F_HC_Genes > 0.95){
		BE_F_HC_Genes <- 0.01
		if (showinfo) {
			message(paste("BE_F_HC_Genes Reset to ", BE_F_HC_Genes, sep=""),appendLF=TRUE)
			flush.console()
		}
	}
	
	#Normalization
	datamatrix <- BatchEffectLinnorm1(datamatrix, minZeroPortion, BE_F_LC_Genes = BE_F_LC_Genes, BE_F_HC_Genes = BE_F_HC_Genes, BE_F_p = BE_F_p, BE_strength = BE_strength)

	colnames(datamatrix) <- CN
	rownames(datamatrix) <- RN
	
	return(datamatrix)
}

