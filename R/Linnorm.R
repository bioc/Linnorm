#' Linnorm Transformation Function
#'
#' This function performs the Linear model and normality based transformation method (Linnorm) for RNA-seq expression data or large scale count data.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets.
#' @param showinfo Logical. Show lambda value calculated. Defaults to FALSE.
#' @param method	Character. "default" or "lambda" The program will output the transformed matrix if the method is "default". If the method is "lambda", the program will output a lambda value.
#' @param minZeroPortion Double >=0, <= 1. For example, setting minZeroPortion as 0.5 will remove genes with more than half data values being zero in the calculation of normalizing parameter. It is strongly suggested to change this to 0.5 for single cell RNA-seq data. Defaults to 2/3.
#' @param keepAll	Logical. After applying minZeroPortion filtering, should Linnorm keep all genes in the results? Defualts to TRUE.
#' @param perturbation Integer >=2. To search for an optimal minimal deviation parameter (please see the article), Linnorm uses the iterated local search algorithm which perturbs away from the initial local minimum. The range of the area searched in each perturbation is exponentially increased as the area get further away from the initial local minimum, which is determined by their index. This range is calculated by 10 * (perturbation ^ index).
#' @details  If method is default, Linnorm outputs a transformed expression matrix. For users who wish to work with lambda instead, the output is a single lambda value. Please note that users with the lambda value can obtain a transformed Linnorm dataset by: log1p(lambda * datamatrix). There is no need to rerun the program if a lambda is already calculated.
#' @return This function returns either a transformed data matrix or a lambda value.
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric
#' @export
#' @examples
#' #Obtain example matrix:
#' data(LIHC)
#' #Transformation:
#' transformedExp <- Linnorm(LIHC)
#' transformedExp <- Linnorm(LIHC, method = "lambda")
#' @import
#' Rcpp
#' RcppArmadillo
Linnorm <- function(datamatrix, showinfo = FALSE, method="default",perturbation=10, minZeroPortion=2/3, keepAll = TRUE) {
	if (keepAll == TRUE) {
		Filter=FALSE
	} else {
		Filter=TRUE
	}
	L_F_p = 0.3173
	L_F_LC_Genes = "Auto"
	L_F_HC_Genes = 0.01
	BE_F_p = 0.3173
	BE_F_LC_Genes = "Auto"
	BE_F_HC_Genes = 0.01
	BE_strength = 0.5
	max_F_LC=0.75

	#data checking
	RN <- rownames(datamatrix)
	CN <- colnames(datamatrix)
	datamatrix <- as.matrix(datamatrix)
	if (length(datamatrix[1,]) < 3) {
		stop("Number of samples is less than 3.")
	}
	if (length(datamatrix[,1]) < 500) {
		stop("Number of features is too small.")
	}
	if (perturbation < 2) {
		stop("perturbation is too small.")
	}
	if (minZeroPortion > 1 || minZeroPortion < 0) {
		stop("Invalid minZeroPortion.")
	}
	if (L_F_p >1 || L_F_p < 0) {
		stop("Invalid L_F_p.")
	}
	if (L_F_LC_Genes > 0.95 || L_F_LC_Genes < 0.01) {
		if (L_F_LC_Genes != "Auto") {
			stop("Invalid L_F_LC_Genes.")
		}
	}
	if (L_F_HC_Genes > 0.95 || L_F_HC_Genes < 0.01) {
		stop("Invalid L_F_HC_Genes.")
	}
	if (BE_F_p > 1 || BE_F_p < 0) {
		stop("Invalid BE_F_p.")
	}
	if (BE_F_LC_Genes > 0.95 || BE_F_LC_Genes < 0.01) {
		if (BE_F_LC_Genes != "Auto") {
			stop("Invalid BE_F_LC_Genes.")
		}
	}
	if (BE_F_HC_Genes > 0.95 || BE_F_HC_Genes < 0.01) {
		stop("Invalid BE_F_HC_Genes.")
	}
	if (BE_strength > 1 | BE_strength < 0) {
		stop("Invalid BE_strength.")
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
	datamatrix <- XPM(datamatrix)
	
	#Find maxBound
	TheMean <- NZrowMeans(datamatrix)
	MeanOrder <- order(TheMean, decreasing = FALSE)
	numZero <- sum(TheMean == 0)
	fivepercent <- floor(0.05 * nrow(datamatrix)) + 1
	nonZero <- datamatrix[MeanOrder[numZero:(numZero +fivepercent)],][which(datamatrix[MeanOrder[numZero:(numZero +fivepercent)],] != 0)]
	maxBound <- length(nonZero)/sum(nonZero)
	
	#Get filter low count genes threhsold
	Keep <- 0
	if (minZeroPortion == 0 || minZeroPortion == 1) {
		Keep <- which(rowSums(datamatrix != 0) >= ncol(datamatrix) * minZeroPortion)
	} else {
		Keep <- which(rowSums(datamatrix != 0) > ncol(datamatrix) * minZeroPortion)
	}
	if (length(Keep) < 200) {
		stop("Given the current minZeroPortion threshold, the number of remaining feature (less than 200) is too small.")
	}
	LC_Threshold <- 0
	if (BE_F_LC_Genes == "Auto" || L_F_LC_Genes == "Auto") {
		LC_Threshold <- FindLCT(datamatrix[Keep,], maxBound)
		if (LC_Threshold > max_F_LC) {
			if (showinfo) {
				message(paste("Filter low count gene threshold is ", LC_Threshold, ". It is larger than max_F_LC, ", max_F_LC, ", which is now used.", sep=""))
				LC_Threshold <- max_F_LC
			}
		}
		if (BE_F_LC_Genes == "Auto") {
			BE_F_LC_Genes <- LC_Threshold
		}
		if (L_F_LC_Genes == "Auto") {
			L_F_LC_Genes <- LC_Threshold
		}
		if (showinfo) {
			message(paste("Filter low count genes threshold is set to ", LC_Threshold, sep=""),appendLF=TRUE)
		}
	}
	if (L_F_LC_Genes + L_F_HC_Genes > 0.95){
		L_F_HC_Genes <- 0.01
		if (showinfo) {
			message(paste("L_F_HC_Genes Reset to ", L_F_HC_Genes, sep=""),appendLF=TRUE)
		}
	}
	if (BE_F_LC_Genes + BE_F_HC_Genes > 0.95){
		BE_F_HC_Genes <- 0.01
		if (showinfo) {
			message(paste("BE_F_HC_Genes Reset to ", BE_F_HC_Genes, sep=""),appendLF=TRUE)
		}
	}
	
	#Filter dataset and calculate lambda
	FilteredData <- FirstFilter(datamatrix, minZeroPortion, L_F_p = L_F_p, L_F_LC_Genes = L_F_LC_Genes, L_F_HC_Genes = L_F_HC_Genes)
	lambda <- LocateLambda(FilteredData, perturbation, maxBound)
	
	#Normalization
	if (BE_strength > 0) {
		datamatrix <- BatchEffectLinnorm1(datamatrix * lambda, minZeroPortion, BE_F_LC_Genes = BE_F_LC_Genes, BE_F_HC_Genes = BE_F_HC_Genes, BE_F_p = BE_F_p, BE_strength = BE_strength)
		colnames(datamatrix) <- CN
		rownames(datamatrix) <- RN
		datamatrix <- log1p(datamatrix)
	} else {
		colnames(datamatrix) <- CN
		rownames(datamatrix) <- RN
		datamatrix <- log1p(datamatrix * lambda)
	}
	
	DataImputation <- FALSE
	if (method=="internal") {
		Filter = TRUE
		minZeroPortion <- 0.5
		DataImputation <- TRUE
	}
	if (Filter || DataImputation) {
		if (Filter) {
			datamatrix <- datamatrix[order(NZrowMeans(datamatrix),decreasing=FALSE),]
			datamatrix <- datamatrix[rowSums(datamatrix != 0) >= ncol(datamatrix) * minZeroPortion,]
			Start <- floor(nrow(datamatrix) * LC_Threshold + 1)
			End <- nrow(datamatrix)
			Keep <- Start:End
			datamatrix <- datamatrix[Keep,]
		}
		if (DataImputation) {
			DIMZP <- 0.5
			method <- "euclidean"
			datamatrix <- Linnorm.DataImput(datamatrix,DIMZP = DIMZP ,method = method)
		}
	}
	if (showinfo) {
		message("Lambda is ", lambda,".",appendLF=TRUE)
		flush.console()
	}
	return (datamatrix)
}

