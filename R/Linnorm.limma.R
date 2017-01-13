#' Linnorm-limma pipeline for Differentially Expression Analysis
#'
#' This function first performs Linnorm transformation on the dataset. Then, it will perform limma for DEG analysis. Please cite both Linnorm and limma when you use this function for publications.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets.
#' @param design	A design matrix required for limma. Please see limma's documentation or our vignettes for more detail.
#' @param MZP Double >=0, <= 1. Minimum non-Zero Portion Threshold for this function. Genes not satisfying this threshold will be removed from HVG anlaysis. For exmaple, if set to 0.3, genes without at least 30 percent of the samples being non-zero will be removed. Defaults to 0.
#' @param output	Character. "DEResults" or "Both". Set to "DEResults" to output a matrix that contains Differential Expression Analysis Results. Set to "Both" to output a list that contains both Differential Expression Analysis Results and the transformed data matrix.
#' @param noINF	Logical. Prevent generating INF in the fold change column by using CPM+1 or TPM+1. If it is set to FALSE, INF will be generated if one of the conditions has zero expression. Defaults to TRUE.
#' @param robust Logical. In the eBayes function of Limma, run with robust setting with TRUE or FALSE. Defaults to FALSE.
#' @param ... arguments that will be passed into Linnorm's transformation function.
#' @details  This function performs both Linnorm and limma for users who are interested in differential expression analysis.
#' @return If output is set to "DEResults", this function will output a matrix with Differntial Expression Analysis Results with the following columns:
##' \itemize{
##'  \item{logFC:}{ Log 2 Fold Change}
##'  \item{XPM:}{ Average Expression. If input is raw count or CPM, this column has the CPM unit. If input is RPKM, FPKM or TPM, this column has the TPM unit.}
##'  \item{t:}{ moderated t-statistic}
##'  \item{P.Value:}{ p value}
##'  \item{adj.P.Val:}{ Adjusted p value. This is also called False Discovery Rate or q value.}
##'  \item{B:}{ log odds that the feature is differential}
##' }
#' @return If output is set to Both, this function will output a list with the following objects:
##' \itemize{
##'  \item{DEResults:}{ Differntial Expression Analysis Results as described above.}
##'  \item{Linnorm:}{ Linnorm transformed data matrix.}
##' }
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric limma
#' @export
#' @examples
#' #Obtain example matrix:
#' data(LIHC)
#' #Create limma design matrix (first 10 columns are tumor, last 10 columns are normal)
#' designmatrix <- c(rep(1,10),rep(2,10))
#' designmatrix <- model.matrix(~ 0+factor(designmatrix))
#' colnames(designmatrix) <- c("group1", "group2")
#' rownames(designmatrix) <- colnames(LIHC)
#' #DEG analysis
#' DEGResults <- Linnorm.limma(LIHC, designmatrix)

Linnorm.limma <- function(datamatrix, design=NULL, MZP = 0, output="DEResults", noINF=TRUE, robust=FALSE, ...) {
	datamatrix <- as.matrix(datamatrix)
	if (is.null(design)) {
		stop("design is null.")
	}
	if (MZP > 1 || MZP < 0) {
		stop("Invalid MZP.")
	}
	CN <- colnames(datamatrix)
	RN <- rownames(datamatrix)
	expdata <- XPM(datamatrix) * 1000000
	colnames(expdata) <- CN
	rownames(expdata) <- RN
	
	#Linnorm transformation
	datamatrix <- Linnorm(datamatrix, ...)
	
	#limma analysis:
	#limma has a lot of warnings, lets turn them off.
	options(warn=-1)
	
	#adjust design for further analysis
	CN <- c()
	for (i in seq_along(design[1,])) {
		CN <- c(CN, paste("group",i,sep=""))
	}
	colnames(design) <- CN
	fit2 <- lmFit(datamatrix,design)
	if (length(design[1,]) == 2) {
		contrast.matrix <- makeContrasts("group1-group2", levels=design)
	} else if (length(design[1,]) == 3) {
		contrast.matrix <- makeContrasts("group1-group2","group1-group3","group2-group3", levels=design)
	} else if (length(design[1,]) == 4) {
		contrast.matrix <- makeContrasts("group1-group2","group1-group3","group1-group4","group2-group3","group2-group4","group3-group4", levels=design)
	} else if (length(design[1,]) == 5) {
		contrast.matrix <- makeContrasts("group1-group2","group1-group3","group1-group4","group1-group5","group2-group3","group2-group4","group2-group5","group3-group4","group3-group5","group4-group5", levels=design)
	} else if (length(design[1,]) == 6) {
		contrast.matrix <- makeContrasts("group1-group2","group1-group3","group1-group4","group1-group5","group1-group6","group2-group3","group2-group4","group2-group5","group2-group6","group3-group4","group3-group5","group3-group6","group4-group5","group4-group6","group5-group6", levels=design)
	} else {
		stop("Error: The number of group is larger than 6. Please use limma manually.")
	}
	fit2 <- contrasts.fit(fit2, contrast.matrix)
	fit2 <- eBayes(fit2,robust=robust)
	limmaResults <- topTable(fit2, number=length(fit2$p.value), adjust.method="BH")
	datamatrix <- datamatrix[rownames(limmaResults),]
	expdata <- expdata[rownames(limmaResults),]
	if (length(design[1,]) == 2) {
		set1 <- as.numeric(which(design[,1] == 1))
		set2 <- as.numeric(which(design[,1] != 1))
		limmaResults[,2] <- rowMeans(expdata)
		if (noINF) {
			expdata <- expdata + 1
			limmaResults[,1] <- log(rowMeans(expdata[,set1])/rowMeans(expdata[,set2]),2)
		} else {
			limmaResults[,1] <- log(rowMeans(expdata[,set1])/rowMeans(expdata[,set2]),2)
		}
		colnames(limmaResults)[2] <- "XPM"
	} else {
		limmaResults[,2] <- rowMeans(expdata)
		colnames(limmaResults)[2] <- "XPM"
		limmaResults <- limmaResults[,-c(1)]
	}
	if (MZP > 0) {
		limmaResults <- limmaResults[rownames(datamatrix[rowSums(datamatrix != 0) >= MZP * ncol(datamatrix),]),]
	}
	if (output=="DEResults") {
		return (limmaResults)
	}
	if (output=="Both") {
		listing <- list(limmaResults, datamatrix)
		results <- setNames(listing, c("DEResults", "Linnorm"))
		return (results)
	}
}


