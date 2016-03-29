#' Linnorm Normalization Function
#'
#' This function performs the Linear model and normality based normalization method (Linnorm) for RNA-seq expression data or large scale count data.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Undefined values such as NA are not supported.
#' @param showinfo Logical. Show lambda value calculated. Defaults to FALSE.
#' @param method	"default" or "lambda" The program will output the transformed matrix if the method is "default". If the method is "lambda", the program will output a lambda value.
#' @param minZeroPortion Double >=0, <= 1. Featuress without at least this portion of non-zero values will not be used in the calculation of normalizing parameter. Defaults to 2/3.
#' @param perturbation Integer >=2. To search for an optimal minimal deviation parameter (please see the article), Linnorm uses the iterated local search algorithm which perturbs away from the initial local minimum. The range of the area searched in each perturbation is exponentially increased as the area get further away from the initial local minimum, which is determined by their index. This range is calculated by 10 * (perturbation ^ index).
#' @details  If method is default, Linnorm outputs a transformed expression matrix. For users who wish to work with lambda instead, the output is a single lambda value. Please note that users with the lambda value can obtain a normalized Linnorm dataset by: log1p(lambda * datamatrix). There is no need to rerun the program if a lambda is already calculated.
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
#' normalizedExp <- Linnorm(expMatrix)
#' normalizedExp <- Linnorm(expMatrix, method = "lambda")
#' @import
#' Rcpp
#' RcppArmadillo
Linnorm <- function(datamatrix, showinfo = FALSE, method="default",perturbation=10, minZeroPortion=2/3) {
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
	if (method != "default" && method != "lambda") {
		stop("Invalid algorithm value.")
	}
	if (minZeroPortion >1 || minZeroPortion < 0) {
		stop("Invalid minZeroPortion.")
	}
	#Step 1: Relative Expression
	#Turn it into relative expression
	for (i in seq_along(expdata[1,])) {
		expdata[,i] <- expdata[,i]/sum(expdata[,i])
	}
	#Save the original dataset before trimming.
	datamatrix <- expdata
	expdata <- expdata[order(rowMeans(expdata)),]	
	if (method == "default") {		
		#trim outliers
		expdata <- expdata[rowSums(expdata != 0) >= (length(expdata[1,]) * minZeroPortion),]
		
		X <- LocateLambda(expdata,perturbation)
		if (showinfo) {
			message("Lambda is ", X,".",appendLF=TRUE)
			flush.console()
		}
		return( log1p(datamatrix * X) )
	}
	if (method == "lambda") {
		#trim outliers
		expdata <- expdata[rowSums(expdata != 0) >= (length(expdata[1,]) * minZeroPortion),]
		
		X <- LocateLambda(expdata,perturbation)
		return( X )
	}
}


#' Linnorm-limma pipeline for Differentially Expression Analysis
#'
#' This function first performs Linnorm normalization on the dataset. Then, it will perform limma for DEG analysis. Finally, it will correct fold change outputs from limma results, that will be wrong otherwise. Please cite both Linnorm and limma when you use this function for publications.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Undefined values such as NA are not supported.
#' @param design	A design matrix required for limma. Please see limma's documentation or our vignettes for more detail.
#' @param output Character. "DEResults" or "Both". Set to "DEResults" to output a matrix that contains Differential Expression Analysis Results. Set to "Both" to output a list that contains both Differential Expression Analysis Results and the transformed data matrix.
#' @param noINF	Logical. Prevent generating INF in the fold change column by using Linnorm's lambda and adding one. If it is set to FALSE, INF will be generated if one of the conditions has zero expression. Defaults to TRUE. 
#' @param showinfo Logical. Show lambda value calculated. Defaults to FALSE.
#' @param minZeroPortion Double >=0, <= 1. Featuress without at least this portion of non-zero values will not be used in the calculation of normalizing parameter. Defaults to 2/3.
#' @param perturbation Integer >=2. To search for an optimal minimal deviation parameter (please see the article), Linnorm uses the iterated local search algorithm which perturbs away from the initial local minimum. The range of the area searched in each perturbation is exponentially increased as the area get further away from the initial local minimum, which is determined by their index. This range is calculated by 10 * (perturbation ^ index).
#' @param robust Logical. In the eBayes function of Limma, run with robust setting with TRUE or FALSE. Defaults to TRUE.
#' @details  This function performs both Linnorm and limma for users who are interested in differential expression analysis. Please note that if you directly use a Linnorm Nomralized dataset with limma, the output fold change and average expression with be wrong. (p values and adj.pvalues will be fine.) This is because the voom-limma pipeline assumes input to be in raw counts. This function is written to fix this problem.
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
##'  \item{TMatrix:}{ A Linnorm Normalized Expression Matrix.}
##' }
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric limma
#' @export
#' @examples
#' #Obtain example matrix:
#' library(seqc)
#' SampleA <- ILM_aceview_gene_BGI[,grepl("A_",colnames(ILM_aceview_gene_BGI))]
#' rownames(SampleA) <- ILM_aceview_gene_BGI[,2]
#' SampleB <- ILM_aceview_gene_BGI[,grepl("B_",colnames(ILM_aceview_gene_BGI))]
#' rownames(SampleB) <- ILM_aceview_gene_BGI[,2]
#' #Extract a portion of the matrix for an example
#' expMatrix <- cbind(SampleA[,1:3], SampleB[,1:3])
#' designmatrix <- c(1,1,1,2,2,2)
#' designmatrix <- model.matrix(~ 0+factor(designmatrix))
#' colnames(designmatrix) <- c("group1", "group2")
#' rownames(designmatrix) <- colnames(expMatrix)
#' 
#' #Example 1
#' DEGResults <- Linnorm.limma(expMatrix, designmatrix)
#' #Example 2
#' DEGResults <- Linnorm.limma(expMatrix, designmatrix, output="Both")
Linnorm.limma <- function(datamatrix, design=NULL, output="DEResults", noINF=TRUE, showinfo = FALSE, perturbation=10, minZeroPortion=2/3, robust=TRUE) {
	expdata <- as.matrix(datamatrix)
	
	#Linnorm normalization
	if (length(expdata[1,]) < 2) {
		stop("Number of samples is less than 2.")
	}
	if (length(expdata[,1]) < 10) {
		stop("Number of features is too small.")
	}
	if (perturbation < 2) {
		stop("perturbation is too small.")
	}
	if (minZeroPortion >1 || minZeroPortion < 0) {
		stop("Invalid minZeroPortion.")
	}
	if (is.null(design)) {
		stop("design is null.")
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
	expdata <- expdata[rowSums(expdata != 0) >= (length(expdata[1,]) * minZeroPortion),]
	
	X <- LocateLambda(expdata,perturbation)
	if (showinfo) {
		message("Lambda is ", X,".",appendLF=TRUE)
		flush.console()
	}
	
	expdata <- log1p(datamatrix * X)
	
	#limma analysis:
	
	#limma has a lot of warnings, lets turn them off.
	options(warn=-1)
	
	#adjust design for further analysis
	CN <- c()
	for (i in seq_along(design[1,])) {
		CN <- c(CN, paste("group",i,sep=""))
	}
	colnames(design) <- CN
	fit2 <- lmFit(expdata,design)
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
		stop("Error: The number of columns in design matrix is larger than 6.")
	}
	fit2 <- contrasts.fit(fit2, contrast.matrix)
	fit2 <- eBayes(fit2,robust=robust)
	limmaResults <- topTable(fit2, number=length(fit2$p.value), adjust.method="BH")
	datamatrix <- datamatrix[rownames(limmaResults),]
	if (length(design[1,]) == 2) {
		set1 <- as.numeric(which(design[,1] == 1))
		set2 <- as.numeric(which(design[,1] != 1))
		limmaResults[,2] <- unlist(apply(datamatrix,1,mean)) * 1000000
		if (noINF) {
			datamatrix <- datamatrix * X + 1
			limmaResults[,1] <- unlist(apply(datamatrix,1,function(x){return(log(mean(x[set1])/mean(x[set2] ),2) )}))
		} else {
			limmaResults[,1] <- unlist(apply(datamatrix,1,function(x){return(log(mean(x[set1])/mean(x[set2]),2) )}))
		}		
		colnames(limmaResults)[2] <- "XPM"
	} else {
		limmaResults[,2] <- unlist(apply(datamatrix,1,mean)) * 1000000
		colnames(limmaResults)[2] <- "XPM"
		limmaResults <- limmaResults[,-c(1)]
	}
	if (output=="DEResults") {
		return (limmaResults)
	}
	if (output=="Both") {
		listing <- list(limmaResults, expdata)
		results <- setNames(listing, c("DEResults", "TMatrix"))
		return (results)
	}
		
	
}


#' This function simulates a RNA-seq dataset based on a given distribution.
#' @param thisdata Matrix:	The matrix or data frame that contains your dataset. Each row is a gene and each column is a replicate. Undefined values such as NA are not supported. This program assumes that all columns are replicates of the same sample.
#' @param distribution Character: Defaults to "Poisson". This parameter controls the output distribution of the simulated RNA-seq dataset. It can be one of "Gamma" (Gamma distribution), "Poisson" (Poisson distribution), "LogNorm" (Log Normal distribution) or "NB" (Negative Binomial distribution).
#' @param NumRep Integer: The number of replicates. This is half of the number of output samples. Defaults to 3.
#' @param NumDiff Integer: The number of Differentially Changed Features. Defaults to 5000.
#' @param NumFea Integer: The number of Total Features. Defaults to 20000.
#' @param showinfo Logical: should we show data information on the console? Defaults to FALSE.
#' @param MaxLibSizelog2FC Double: The maximum library size difference from the mean that is allowed, in terms of log 2 fold change. Set to 0 to prevent program from generating library size differences. Defaults to 0.5.
#' @return This function returns a list that contains a matrix of count data in integer raw count and a vector that shows which genes are differentially expressed. In the matrix, each row is a gene and each column is a replicate. The first NumRep (see parameter) of the columns belong to sample 1, and the last NumRep (see parameter) of the columns belong to sample 2. There will be NumFea (see parameter) number of rows.
#' @keywords RNA-seq Raw Count Expression Simulation Gamma distribution Simulate Poisson "Log Normal" "Negative Binomial"
#' @export
#' @examples
#' #Obtain example matrix:
#' library(seqc)
#' SampleA <- ILM_aceview_gene_BGI[,grepl("A_",colnames(ILM_aceview_gene_BGI))]
#' rownames(SampleA) <- ILM_aceview_gene_BGI[,2]
#' #Extract a portion of the matrix for an example
#' expMatrix <- SampleA[,1:10]
#' 
#' Example for Negative Binomial distribution
#' simulateddata <- RnaXSim(expMatrix, distribution="NB", NumRep=3, NumDiff = 500, NumFea = 3000)
#' Example for Poisson distribution
#' simulateddata <- RnaXSim(expMatrix, distribution="Poisson", NumRep=3, NumDiff = 500, NumFea = 3000)
#' Example for Log Normal distribution
#' simulateddata <- RnaXSim(expMatrix, distribution="LogNorm", NumRep=3, NumDiff = 500, NumFea = 3000)
#' Example for Gamma distribution
#' simulateddata <- RnaXSim(expMatrix, distribution="Gamma", NumRep=3, NumDiff = 500, NumFea = 3000)
RnaXSim <- function(thisdata, distribution="Poisson", NumRep=3, NumDiff = 5000, NumFea = 20000, showinfo=FALSE, MaxLibSizelog2FC=0.5) {
	thisdata <- na.omit(as.matrix(thisdata))
	if (distribution != "Gamma" && distribution != "Poisson" && distribution != "LogNorm" && distribution != "NB") {
		stop("Invalid distribution.")
	}
	if (NumDiff < 0) {
		stop("Invalid NumDiff value.")
	}
	if (NumRep < 2) {
		stop("Invalid NumRep value.")
	}
	if (NumFea < 0) {
		stop("Invalid NumFea value.")
	}
	
	#Turn it into relative expression
	LibSize <- colSums(thisdata)
	for (i in seq_along(thisdata[1,])) {
		thisdata[,i] <- (thisdata[,i] * min(LibSize))/sum(thisdata[,i])
	}
	#sort and remove features with zeros and less than one.
	thisdata <- thisdata[rowMeans(thisdata) >= 1,]
	thisdata <- thisdata[order(rowMeans(thisdata)),]
	
	if (distribution == "Gamma") {		
		return (GammaSim(thisdata, NumRep=NumRep, NumDiff = NumDiff, NumFea = NumFea, showinfo=showinfo, MaxLibSizelog2FC=MaxLibSizelog2FC))
	}

	if (distribution == "Poisson") {			
		return (PoissonSim(thisdata, NumRep=NumRep, NumDiff = NumDiff, NumFea = NumFea, showinfo=showinfo, MaxLibSizelog2FC=MaxLibSizelog2FC))
	}

	if (distribution == "LogNorm") {
		return (LogNormSim(thisdata, NumRep=NumRep, NumDiff = NumDiff, NumFea = NumFea, showinfo=showinfo, MaxLibSizelog2FC=MaxLibSizelog2FC))
	}

	if (distribution == "NB") {
		return (NBSim(thisdata, NumRep=NumRep, NumDiff = NumDiff, NumFea = NumFea, showinfo=showinfo, MaxLibSizelog2FC=MaxLibSizelog2FC))
	}
	
}

