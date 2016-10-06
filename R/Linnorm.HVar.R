#' Linnorm-Hvar pipeline for highly variable gene discovery.
#'
#' This function first performs Linnorm transformation on the dataset. Then, it will perform highly variable gene discovery.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets. If a Linnorm transfored dataset is being used, please set the "input" argument into "Linnorm".
#' @param input	Character. "Raw" or "Linnorm". In case you have already transformed your dataset with Linnorm, set input into "Linnorm" so that you can input the Linnorm transformed dataset into the "datamatrix" argument. Defaults to "Raw".
#' @param spikein	character vector. Row names of the spike-in genes in the datamatrix. If this is provided, test of significance will be performed against the spike in genes. Defaults to NULL.
#' @param method	Character. "SE" or "SD". Use Standard Error (SE) or Standard Deviation (SD) to calculate p values. Defaults to SD.
#' @param showinfo	Logical. Show lambda value calculated. Defaults to FALSE.
#' @param minZeroPortion	Double >=0, <= 1. For example, setting minZeroPortion as 0.5 will remove genes with more than half data values being zero in the calculation of normalizing parameter. Since this test is based on variance, which requires more non-zero values, it is suggested to set it to a larger value. Defaults to 2/3.
#' @param keepAll	Logical. After applying minZeroPortion filtering, should Linnorm keep all genes in the results? Defualts to FALSE.
#' @param perturbation	Integer >=2. To search for an optimal minimal deviation parameter (please see the article), Linnorm uses the iterated local search algorithm which perturbs away from the initial local minimum. The range of the area searched in each perturbation is exponentially increased as the area get further away from the initial local minimum, which is determined by their index. This range is calculated by 10 * (perturbation ^ index).
#' @param log.p	Logical. Output p/q values in log scale. Defaults to FALSE.
#' @param sig.value	Character. "p" or "q". Use p or q value for highlighting significant genes. Defaults to "p".
#' @param sig	Double >0, <= 1. Significant level of p or q value for plotting. Defaults to 0.05.
#' @details  This function discovers highly variable gene in the dataset using Linnorm transformation.
#' @return This function will output a list with the following objects:
##' \itemize{
##'  \item{Results:}{ A matrix with the results.}
##'  \item{plot:}{ Mean vs Standard Deviation Plot which highlights significant genes.}
##'  \item{Linnorm:}{ Linnorm transformed and filtered data matrix.}
##' }
#' @return The Results matrix has the following columns:
##' \itemize{
##'  \item{XPM:}{ Average non-zero expression level in XPM. If input is raw coutns or CPM, this column is in CPM unit. If input is RPKM, FPKM or TPM, this column is in the TPM unit.}
##'  \item{XPM.SD:}{ Standard deviation of average non-zero expression.}
##'  \item{Transformed.Avg.Exp:}{ Average expression of non-zero Linnorm transformed data.}
##'  \item{Transformed.SD:}{ Standard deviation of non-zero Linnorm transformed data.}
##'  \item{Normalized.Log2.SD.Fold.Change:}{ Normalized log2 fold change of the gene's standard deviation.}
##'  \item{p.value:}{ p value of the statistical test.}
##'  \item{q.value:}{ q value/false discovery rate/adjusted p value of the statistical test.}
##' }
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric variance highly variable
#' @export
#' @examples
#' data(Islam2011)
#' results <- Linnorm.HVar(Islam2011)

Linnorm.HVar <- function(datamatrix, input="Raw", method = "SD", spikein=NULL, showinfo = FALSE, perturbation=10, minZeroPortion=2/3, keepAll=FALSE, log.p=FALSE, sig.value="p", sig=0.05) {
	if (input != "Raw" && input != "Linnorm") {
		stop("input argument is not recognized.")
	}
	if (method != "SE" && method != "SD") {
		stop("method is not recognized.")
	}
	if (sig <= 0 || sig > 1) {
		stop("Invalid sig value.")
	}
	if ( sig.value != "p" &&  sig.value != "q") {
		stop("Invalid sig.value.")
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
		XPM <- NZrowMeans(expdata * 1000000) 
		XPMSD <- NZrowSDs(expdata * 1000000)
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
		XPM <- NZrowMeans(XPMdata) 
		XPMSD <- NZrowSDs(XPMdata)
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
	datamean <- NZrowMeans(expdata)
	dataSD <- sqrt(NZrowSDs(expdata))
	
	#Logistic regression to fit technical noise.
	logitit <- loess(dataSD~datamean, weights=datamean)
	
	#In case of negative fit, where negative stdev is impossible in our case, change them into the smallest stdev in the dataset
	logitit$fitted[which(logitit$fitted <= 0)] <- min(dataSD[which(dataSD != 0)])
	
	#Obtain Stdev ratios to adjust for technical noise.
	SDRatio <- dataSD/logitit$fitted

	#Calculate p values
	#if spike in list is provided, we test whether a given standard deviation is larger than the spike in.
	pvalues <- 0
	spikes <- 0
	
	
	if (length(spikein) < 2) {
		if (showinfo) {
			message("Length of spikein <= 2. Not used.",appendLF=TRUE)
			flush.console()
		}
		#Remove outlier
		SDRatio2 <- SDRatio[!SDRatio %in% boxplot.stats(SDRatio)$out]
		if (method == "SD") {
			pvalues <- pnorm(as.numeric(SDRatio),mean(as.numeric(SDRatio2),na.rm=TRUE),sd(SDRatio2,na.rm=TRUE), lower.tail = FALSE, log.p=TRUE )
		}
		if (method == "SE") {
			pvalues <- pnorm(as.numeric(SDRatio),mean(as.numeric(SDRatio2),na.rm=TRUE),sd(SDRatio2,na.rm=TRUE)/sqrt(length(SDRatio2)), lower.tail = FALSE, log.p=TRUE )
		}
	} else {
		spikes <- which(rownames(expdata) %in% spikein)
		SDRatio2 <- SDRatio[spikes]
		SDRatio2 <- SDRatio2[!SDRatio2 %in% boxplot.stats(SDRatio2)$out]
		if (method == "SD") {
			pvalues <- pnorm(as.numeric(SDRatio),mean(as.numeric(SDRatio2),na.rm=TRUE),sd(SDRatio2,na.rm=TRUE), lower.tail = FALSE, log.p=TRUE )
		}
		if (method == "SE") {
			pvalues <- pnorm(as.numeric(SDRatio),mean(as.numeric(SDRatio2),na.rm=TRUE),sd(SDRatio2,na.rm=TRUE)/sqrt(length(spikes)), lower.tail = FALSE, log.p=TRUE )
		}
	}
	epvalues <- exp(pvalues)
	qvalues <- p.adjust(epvalues,"BH")

	results <- matrix(ncol=7, nrow=length(SDRatio))
	colnames(results) <- c("XPM", "XPM.SD", "Transformed.Avg.Exp", "Transformed.SD", "Normalized.Log2.SD.Fold.Change", "p.value", "q.value")
	rownames(results) <- rownames(expdata)
	dataSD <- dataSD^2
	results[,1] <- XPM
	results[,2] <- XPMSD
	results[,3] <- datamean
	results[,4] <- dataSD
	results[,5] <- log(SDRatio,2)
	if (log.p) {
		results[,6] <- pvalues
		results[,7] <- log(qvalues)
	} else {
		results[,6] <- epvalues
		results[,7] <- qvalues
	}
		
	if (keepAll) {
		ZERO <- matrix(0, ncol=7, nrow=length(filteredXPM))
		colnames(ZERO) <- c("XPM", "XPM.SD", "Transformed.Avg.Exp", "Transformed.SD", "Normalized.Log2.SD.Fold.Change", "p.value", "q.value")
		rownames(ZERO) <- rownames(filtered)
		ZERO[,1] <- filteredXPM
		ZERO[,2] <- filteredXPMSD
		ZERO[,3] <- NZrowMeans(filtered)
		ZERO[,4] <- NZrowSDs(filtered)
		ZERO[,5] <- NA
		ZERO[,6] <- NA
		ZERO[,7] <- NA
		results <- rbind(results, ZERO)
		expdata <- rbind(expdata, filtered)
	}

	groups <- rep("non-sig", length(SDRatio))
	if (sig.value == "p") {
		groups[which(epvalues <= sig)] <- "Significant"
	}
	if (sig.value == "q") {
		groups[which(qvalues <= sig)] <- "Significant"
	}	
	groups[spikes] <- "Spike in"
	plotdata <- data.frame(mean=datamean,SD=dataSD,group=groups)
	render_plot <- ggplot_build(ggplot(plotdata, aes(x=mean, y=SD, color=group)) + geom_point(aes(shape=group), size = 1) + scale_x_continuous("Transformed Mean") + scale_y_continuous("Transformed Standard Deviation") + ggtitle("Mean vs SD plot of highly variable genes") + theme(aspect.ratio=3/4))
	listing <- list(results, render_plot, expdata)
	result <- setNames(listing, c("Results", "plot", "Linnorm"))
	return (result)
}


