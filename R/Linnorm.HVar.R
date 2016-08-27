#' Linnorm-Hvar pipeline for highly variable gene discovery.
#'
#' This function first performs Linnorm transformation on the dataset. Then, it will perform highly variable gene discovery.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets. If a Linnorm transfored dataset is being used, please set the "input" argument into "Linnorm".
#' @param input	Character. "Raw" or "Linnorm". In case you have already transformed your dataset with Linnorm, set input into "Linnorm" so that you can input the Linnorm transformed dataset into the "datamatrix" argument. Defaults to "Raw".
#' @param spikein	character vector. Row names of the spike-in genes in the datamatrix. If this is provided, test of significance will be performed against the spike in genes. Defaults to NULL.
#' @param method	Character. "SE" or "SD". Use Standard Error (SE) or Standard Deviation (SD) to calculate p values. Defaults to SE.
#' @param showinfo	Logical. Show lambda value calculated. Defaults to FALSE.
#' @param minZeroPortion	Double >=0, <= 1. For example, setting minZeroPortion as 0.5 will remove genes with more than half data values being zero in the calculation of normalizing parameter. Since this test is based on variance, which requires more non-zero values, it is suggested to set it to a larger value. Defaults to 0.5.
#' @param keepAll	Logical. After applying minZeroPortion filtering, should Linnorm keep all genes in the results? Defualts to FALSE.
#' @param perturbation	Integer >=2. To search for an optimal minimal deviation parameter (please see the article), Linnorm uses the iterated local search algorithm which perturbs away from the initial local minimum. The range of the area searched in each perturbation is exponentially increased as the area get further away from the initial local minimum, which is determined by their index. This range is calculated by 10 * (perturbation ^ index).
#' @param log.p	Logical. Output p/q values in log scale. Defaults to FALSE.
#' @param sig.q	Double >0, <= 1. Significant level of q value for plotting. Defaults to 0.01.
#' @details  This function discovers highly variable gene in the dataset using Linnorm transformation.
#' @return This function will output a list with the following objects:
##' \itemize{
##'  \item{Results:}{ A matrix with the results.}
##'  \item{plot:}{ Mean vs Standard Deviation Plot which highlights significant genes.}
##'  \item{Linnorm:}{ Linnorm transformed and filtered data matrix.}
##' }
#' @return The Results matrix has the following columns:
##' \itemize{
##'  \item{XPM:}{ Average expression level in XPM. If input is raw coutns or CPM, this column is in CPM unit. If input is RPKM, FPKM or TPM, this column is in the TPM unit.}
##'  \item{XPM.SD:}{ Standard deviation of average expression.}
##'  \item{Transformed.Avg.Exp:}{ Average expression of Linnorm transformed data.}
##'  \item{Transformed.SD:}{ Standard deviation of Linnorm transformed data.}
##'  \item{Normalized.Log2.SD.Fold.Change:}{ Normalized log2 fold change of the gene's standard deviation.}
##'  \item{p.value:}{ p value of the statistical test.}
##'  \item{q.value:}{ q value/false discovery rate/adjusted p value of the statistical test.}
##' }
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric variance highly variable
#' @export
#' @examples
#' data(Islam2011)
#' results <- Linnorm.HVar(Islam2011)

Linnorm.HVar <- function(datamatrix, input="Raw", method = "SE", spikein=NULL, showinfo = FALSE, perturbation=10, minZeroPortion=0.5, keepAll=FALSE, log.p=FALSE, sig.q=0.01) {
	if (input != "Raw" && input != "Linnorm") {
		stop("input argument is not recognized.")
	}
	if (method != "SE" && method != "SD") {
		stop("method is not recognized.")
	}
	if (sig.q <= 0 || sig.q > 1) {
		stop("Invalid sig.q value.")
	}
	expdata <- 0
	XPM <- 0
	XPMSD <- 0
	filtered <- 0
	filteredXPM <- 0
	filteredXPMSD <- 0
	filtered <- 0
	if (input == "Raw") {
		#Linnorm transformation
		expdata <- Linnorm(datamatrix, showinfo = showinfo, method="internal",perturbation=perturbation, minZeroPortion = minZeroPortion, keepAll = keepAll)
		X <- expdata[[2]]
		expdata <- expdata[[1]]
		XPM <- rowMeans(expdata * 1000000) 
		XPMSD <- rowSDs(expdata * 1000000)
		expdata <- log1p(expdata * X)
		if (keepAll) {
			filteredXPM <- XPM[rowSums(expdata != 0) < ncol(expdata) * minZeroPortion]
			filteredXPMSD <- XPMSD[rowSums(expdata != 0) < ncol(expdata) * minZeroPortion]
			filtered <- expdata[rowSums(expdata != 0) < ncol(expdata) * minZeroPortion,]
		}
		XPMSD <- XPMSD[rowSums(expdata != 0) >= ncol(expdata) * minZeroPortion]
		XPM <- XPM[rowSums(expdata != 0) >= ncol(expdata) * minZeroPortion]
		expdata <- expdata[rowSums(expdata != 0) >= ncol(expdata) * minZeroPortion,]
	} 
	if (input == "Linnorm"){
		XPMdata <- exp(datamatrix) - 1
		for (i in seq_along(XPMdata[1,])) {
			XPMdata[,i] <- (XPMdata[,i] * 1000000)/sum(XPMdata[,i])
		}
		XPM <- rowMeans(XPMdata) 
		XPMSD <- rowSDs(XPMdata)
		expdata <- datamatrix
		if (keepAll) {
			filteredXPM <- XPM[rowSums(datamatrix != 0) < ncol(datamatrix) * minZeroPortion]
			filteredXPMSD <- XPMSD[rowSums(datamatrix != 0) < ncol(datamatrix) * minZeroPortion]
			filtered <- expdata[rowSums(datamatrix != 0) < ncol(datamatrix) * minZeroPortion,]
		}
		XPM <- XPM[rowSums(datamatrix != 0) >= ncol(datamatrix) * minZeroPortion]
		XPMSD <- XPMSD[rowSums(datamatrix != 0) >= ncol(datamatrix) * minZeroPortion]
		expdata <- expdata[rowSums(datamatrix != 0) >= ncol(datamatrix) * minZeroPortion,]
	}
	######First use Linnorm transformed dataset######
	datamean <- rowMeans(expdata)
	dataSD <- sqrt(rowSDs(expdata))
	
	#Logistic regression to fit technical noise.
	logitit <- loess(dataSD~datamean)
	
	#In case of negative fit, where negative stdev is impossible in our case, change them into the smallest stdev in the dataset
	logitit$fitted[which(logitit$fitted <= 0)] <- min(dataSD[which(dataSD != 0)])
	
	#Obtain Stdev ratios to adjust for technical noise.
	SDRatio <- log(dataSD/logitit$fitted, 2)

	#Calculate p values
	#if spike in list is provided, we test whether a given standard deviation is larger than the spike in.
	pvalues <- 0
	spikes <- 0
	if (is.null(spikein)) {
		if (method == "SD") {
			pvalues <- pnorm(as.numeric(SDRatio),mean(as.numeric(SDRatio),na.rm=TRUE),sd(SDRatio,na.rm=TRUE), lower.tail = FALSE, log.p=TRUE )
		}
		if (method == "SE") {
			pvalues <- pnorm(as.numeric(SDRatio),mean(as.numeric(SDRatio),na.rm=TRUE),sd(SDRatio,na.rm=TRUE)/sqrt(length(SDRatio)), lower.tail = FALSE, log.p=TRUE )
		}
	} else {
		spikes <- which(rownames(expdata) %in% spikein)
		if (method == "SD") {
			pvalues <- pnorm(as.numeric(SDRatio),mean(as.numeric(SDRatio[spikes]),na.rm=TRUE),sd(SDRatio[spikes],na.rm=TRUE), lower.tail = FALSE, log.p=TRUE )
		}
		if (method == "SE") {
			pvalues <- pnorm(as.numeric(SDRatio),mean(as.numeric(SDRatio[spikes]),na.rm=TRUE),sd(SDRatio[spikes],na.rm=TRUE)/sqrt(length(spikes)), lower.tail = FALSE, log.p=TRUE )
		}
	}
	qvalues <- p.adjust(exp(pvalues),"BH")

	results <- matrix(ncol=7, nrow=length(SDRatio))
	colnames(results) <- c("XPM", "XPM.SD", "Transformed.Avg.Exp", "Transformed.SD", "Normalized.Log2.SD.Fold.Change", "p.value", "q.value")
	rownames(results) <- rownames(expdata)
	dataSD <- dataSD^2
	results[,1] <- XPM
	results[,2] <- XPMSD
	results[,3] <- datamean
	results[,4] <- dataSD
	results[,5] <- SDRatio
	if (log.p) {
		results[,6] <- pvalues
		results[,7] <- log(qvalues)
	} else {
		results[,6] <- exp(pvalues)
		results[,7] <- qvalues
	}
		
	if (keepAll) {
		ZERO <- matrix(0, ncol=7, nrow=length(filteredXPM))
		colnames(ZERO) <- c("XPM", "XPM.SD", "Transformed.Avg.Exp", "Transformed.SD", "Normalized.Log2.SD.Fold.Change", "p.value", "q.value")
		rownames(ZERO) <- rownames(filtered)
		ZERO[,1] <- filteredXPM
		ZERO[,2] <- filteredXPMSD
		ZERO[,3] <- rowMeans(filtered)
		ZERO[,4] <- rowSDs(filtered)
		ZERO[,5] <- NA
		ZERO[,6] <- NA
		ZERO[,7] <- NA
		results <- rbind(results, ZERO)
		expdata <- rbind(expdata, filtered)
	}

	groups <- rep("non-sig", length(SDRatio))
	groups[which(qvalues <= sig.q)] <- "Significant"
	groups[spikes] <- "Spike in"
	plotdata <- data.frame(mean=datamean,SD=dataSD,group=groups)
	render_plot <- ggplot_build(ggplot(plotdata, aes(x=mean, y=SD, color=group)) + geom_point(aes(shape=group), size = 1) + scale_x_continuous("Transformed Mean") + scale_y_continuous("Transformed Standard Deviation") + ggtitle("Mean vs SD plot of highly variable genes") + theme(aspect.ratio=3/4))
	listing <- list(results, render_plot, expdata)
	result <- setNames(listing, c("Results", "plot", "Linnorm"))
	return (result)
}


