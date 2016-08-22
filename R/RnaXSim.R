#' This function simulates a RNA-seq dataset based on a given distribution.
#' @param datamatrix	Matrix. The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets.This program assumes that all columns are replicates of the same sample.
#' @param distribution	Character: Defaults to "Poisson". This parameter controls the output distribution of the simulated RNA-seq dataset. It can be one of "Gamma" (Gamma distribution), "Poisson" (Poisson distribution), "LogNorm" (Log Normal distribution) or "NB" (Negative Binomial distribution).
#' @param NumRep	Integer: The number of replicates. This is half of the number of output samples. Defaults to 3.
#' @param NumDiff	Integer: The number of Differentially Changed Features. Defaults to 2000.
#' @param NumFea	Integer: The number of Total Features. Defaults to 20000.
#' @param showinfo	Logical: should we show data information on the console? Defaults to FALSE.
#' @param DEGlog2FC	"Auto" or Double: log 2 fold change threshold that defines differentially expressed genes. If set to "Auto," DEGlog2FC is defined at the level where ANOVA can get a q value of 0.05 with the average expression, where the data values are log1p transformed. Defaults to "Auto".
#' @param MaxLibSizelog2FC	Double: The maximum library size difference from the mean that is allowed, in terms of log 2 fold change. Set to 0 to prevent program from generating library size differences. Defaults to 0.5.
#' @return This function returns a list that contains a matrix of count data in integer raw count and a vector that shows which genes are differentially expressed. In the matrix, each row is a gene and each column is a replicate. The first NumRep (see parameter) of the columns belong to sample 1, and the last NumRep (see parameter) of the columns belong to sample 2. There will be NumFea (see parameter) number of rows. The top NumCorr of genes will be positively or negatively correlated with each other (randomly); and they are evenly separated into groups. Each group is not intended to be correlated to each other, but, by chance, it can happen.
#' @keywords RNA-seq Raw Count Expression Simulation Gamma distribution Simulate Poisson "Log Normal" "Negative Binomial"
#' @export
#' @examples
#' #Obtain example matrix:
#' library(seqc)
#' SampleA <- ILM_aceview_gene_BGI[,grepl("A_",colnames(ILM_aceview_gene_BGI))]
#' rownames(SampleA) <- ILM_aceview_gene_BGI[,2]
#' #Extract a portion of the matrix for an example
#' expMatrix <- SampleA[,1:10]
#' #Example for Negative Binomial distribution
#' simulateddata <- RnaXSim(expMatrix, distribution="NB", NumRep=3, NumDiff = 200, NumFea = 2000)
#' #Example for Poisson distribution
#' simulateddata <- RnaXSim(expMatrix, distribution="Poisson", NumRep=3, NumDiff = 200, NumFea = 2000)
#' #Example for Log Normal distribution
#' simulateddata <- RnaXSim(expMatrix, distribution="LogNorm", NumRep=3, NumDiff = 200, NumFea = 2000)
#' #Example for Gamma distribution
#' simulateddata <- RnaXSim(expMatrix, distribution="Gamma", NumRep=3, NumDiff = 200, NumFea = 2000)
RnaXSim <- function(datamatrix, distribution="Poisson", NumRep=3, NumDiff = 2000, NumFea = 20000, showinfo=FALSE, DEGlog2FC="Auto", MaxLibSizelog2FC=0.5) {
	datamatrix <- na.omit(as.matrix(datamatrix))
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
	if (anyNA(datamatrix)) {
		stop("Dataset contains NA.")
	}
	if (sum(which(datamatrix < 0)) != 0) {
		stop("Dataset contains negative number.")
	}
	
	#Turn it into relative expression
	LibSize <- colSums(datamatrix)
	for (i in seq_along(datamatrix[1,])) {
		datamatrix[,i] <- (datamatrix[,i] * min(LibSize))/sum(datamatrix[,i])
	}
	#sort and remove features with zeros and less than one.
	#datamatrix <- datamatrix[rowMeans(datamatrix) >= 1,]
	datamatrix <- datamatrix[rowSums(datamatrix != 0) == length(datamatrix[1,]),]

	datamatrix <- datamatrix[order(rowMeans(datamatrix)),]
	
	if (distribution == "Gamma") {
		return (GammaSim(datamatrix, NumRep=NumRep, NumDiff = NumDiff, NumFea = NumFea, showinfo=showinfo, MaxLibSizelog2FC=MaxLibSizelog2FC, DEGlog2FC=DEGlog2FC))
	}

	if (distribution == "Poisson") {
		return (PoissonSim(datamatrix, NumRep=NumRep, NumDiff = NumDiff, NumFea = NumFea, showinfo=showinfo, MaxLibSizelog2FC=MaxLibSizelog2FC, DEGlog2FC=DEGlog2FC))
	}

	if (distribution == "LogNorm") {
		return (LogNormSim(datamatrix, NumRep=NumRep, NumDiff = NumDiff, NumFea = NumFea, showinfo=showinfo, MaxLibSizelog2FC=MaxLibSizelog2FC, DEGlog2FC=DEGlog2FC))
	}

	if (distribution == "NB") {
		return (NBSim(datamatrix, NumRep=NumRep, NumDiff = NumDiff, NumFea = NumFea, showinfo=showinfo, MaxLibSizelog2FC=MaxLibSizelog2FC, DEGlog2FC=DEGlog2FC))
	}
	
}