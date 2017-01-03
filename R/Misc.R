#Convert Dataset into XPM
XPM <- function(x) {
    .Call(XPMCpp, x)
}

FindLCT <- function(datamatrix, Multy,showinfo) {
	MeanSD <- rowLog1pMeanSD(datamatrix,Multy)
	LC_Threshold <- 0
	MeanSD <- MeanSD[,which(!is.nan(MeanSD[2,]))]
	Portion <- 1/3
	MeanOrder <- order(MeanSD[1,])
	Slope <- 1
	while(Slope > 0 && LC_Threshold < 1) {
		LC_Threshold <- LC_Threshold + 0.01
		Range <- floor(ncol(MeanSD) * LC_Threshold + 1):ncol(MeanSD)
		Range <- Range[1:floor(length(Range) * Portion + 1)]
		Slope <- getSlope(MeanSD[1,MeanOrder[Range]],MeanSD[2,MeanOrder[Range]])
	}
	return (round(LC_Threshold,2))
}

#Filter dataset
FirstFilter <- function(x, minZeroPortion, L_F_p = 0.25, L_F_LC_Genes = 0.01, L_F_HC_Genes = 0.01, spikein = NULL) {
	MeanSDSkew <- NZrowLogMeanSDSkew(x)
	a <- which(!is.nan(MeanSDSkew[3,]))
	MeanSDSkew <- MeanSDSkew[,a]
	x <- x[a,]
	
	#Sort data and filter LowGeneFil and minZeroPortion
	if (minZeroPortion == 0) {
		Keep <- which(rowSums(x != 0) >= ncol(x) * minZeroPortion)
		MeanSDSkew <- MeanSDSkew[,Keep]
		x <- x[Keep,]
	} else {
		Keep <- which(rowSums(x != 0) > ncol(x) * minZeroPortion)
		MeanSDSkew <- MeanSDSkew[,Keep]
		x <- x[Keep,]
	}
	
	MeanOrder <- order(MeanSDSkew[1,], decreasing = FALSE)
	x <- x[MeanOrder,]
	MeanSDSkew <- MeanSDSkew[,MeanOrder]
	
	Start <- floor(nrow(x) * L_F_LC_Genes + 1)
	End <- nrow(x) - floor(nrow(x) * L_F_HC_Genes)
	Keep <- Start:End
	x <- x[Keep,]
	MeanSDSkew <- MeanSDSkew[,Keep]
	spikein <- spikein[which(spikein %in% rownames(x))]
	
	
	allStableGenes <- 0
	logitit <- loessFit(MeanSDSkew[2,],MeanSDSkew[1,], weights=1/MeanSDSkew[2,]^2)
	Keep <- which(logitit$fitted > 0)
	LogFit <- logitit$fitted[Keep]
	x <- x[Keep,]
	MeanSDSkew <- MeanSDSkew[,Keep]
	SDRatio <- log(as.numeric(MeanSDSkew[2,]/LogFit))
	
	#normalize SDRatio
	LR <- LinearRegression(MeanSDSkew[1,],SDRatio)
	Residual <- (SDRatio - (LR$coefficients[[2]] * MeanSDSkew[1,] + LR$coefficients[[1]]))
	LR2 <- LinearRegression(MeanSDSkew[1,],abs(Residual))
	SDRatio <- SDRatio * (LR2$coefficients[[2]] * MeanSDSkew[1,1] + LR2$coefficients[[1]])/(LR2$coefficients[[2]] * MeanSDSkew[1,] + LR2$coefficients[[1]])
	
	
	pvalueMatrix <- matrix(nrow=ncol(MeanSDSkew), ncol=2)
	
	#Column 1, SD p values. Column 2, Skew p values
	if (length(spikein) < 10) {
		warning("Not enough sufficiently expressed spike-in genes (less than 10), they will be ignored.")
	
		SDnoOutlier <- SDRatio[!SDRatio %in% boxplot.stats(SDRatio)$out]
		Themean <- mean(SDnoOutlier,na.rm=TRUE)
		TheSD <- sd(SDnoOutlier,na.rm=TRUE)
		Bigger <- which(SDRatio >= Themean)
		smaller <- which(SDRatio < Themean)
		pvalues <- pnorm(SDRatio,Themean,TheSD, lower.tail = FALSE)
		pvalueMatrix[Bigger,1] <- 2 * pvalues[Bigger]
		pvalueMatrix[smaller,1] <- 2 * (1 - pvalues[smaller])
		#pvalueMatrix[,1] <- pnorm(SDRatio,mean(SDRatio,na.rm=TRUE),sd(SDRatio,na.rm=TRUE), lower.tail = FALSE )
		
		SkewLR <- LinearRegression(MeanSDSkew[1,],MeanSDSkew[3,])
		SkewResidual <- MeanSDSkew[3,] - (SkewLR$coefficients[[2]] * MeanSDSkew[3,] - SkewLR$coefficients[[1]])
		SkewnoOutlier <- SkewResidual[!SkewResidual %in% boxplot.stats(SkewResidual)$out]
		
		Themean <- mean(SkewnoOutlier,na.rm=TRUE)
		TheSD <- sd(SkewnoOutlier,na.rm=TRUE)
		Bigger <- which(SkewResidual >= Themean)
		smaller <- which(SkewResidual < Themean)
		pvalues <- pnorm(SkewResidual,Themean,TheSD, lower.tail = FALSE)
		pvalueMatrix[Bigger,2] <- 2 * pvalues[Bigger]
		pvalueMatrix[smaller,2] <- 2 * (1 - pvalues[smaller])
	} else {
		spikes <- which(rownames(x) %in% spikein)
		SDnoOutlier <- SDRatio[spikes]
		Themean <- mean(SDnoOutlier,na.rm=TRUE)
		TheSD <- sd(SDnoOutlier,na.rm=TRUE)
		Bigger <- which(SDRatio >= Themean)
		smaller <- which(SDRatio < Themean)
		pvalues <- pnorm(SDRatio,Themean,TheSD, lower.tail = FALSE)
		pvalueMatrix[Bigger,1] <- 2 * pvalues[Bigger]
		pvalueMatrix[smaller,1] <- 2 * (1 - pvalues[smaller])
		#pvalueMatrix[,1] <- pnorm(SDRatio,mean(SDRatio,na.rm=TRUE),sd(SDRatio,na.rm=TRUE), lower.tail = FALSE )
		
		SkewLR <- LinearRegression(MeanSDSkew[1,],MeanSDSkew[3,])
		SkewResidual <- MeanSDSkew[3,] - (SkewLR$coefficients[[2]] * MeanSDSkew[3,] - SkewLR$coefficients[[1]])
		SkewnoOutlier <- SkewResidual[spikes]
		
		Themean <- mean(SkewnoOutlier,na.rm=TRUE)
		TheSD <- sd(SkewnoOutlier,na.rm=TRUE)
		Bigger <- which(SkewResidual >= Themean)
		smaller <- which(SkewResidual < Themean)
		pvalues <- pnorm(SkewResidual,Themean,TheSD, lower.tail = FALSE)
		pvalueMatrix[Bigger,2] <- 2 * pvalues[Bigger]
		pvalueMatrix[smaller,2] <- 2 * (1 - pvalues[smaller])
	}
	
	#combinedPvalues <- apply(pvalueMatrix,1,FisherMethod)
	combinedPvalues <- empiricalBrownsMethod(MeanSDSkew[2:3,],pvalueMatrix)
	allStableGenes <- which(combinedPvalues > L_F_p)
	
	#allStableGenes <- which(pvalueMatrix[,1] < L_F_p)
	#allStableGenes <- which(allStableGenes %in% which(pvalueMatrix[,2] < L_F_p))
	#allStableGenes <- which(!((1:nrow(pvalueMatrix)) %in% allStableGenes))
	
	#Safety, need 100 genes at least
	while (length(allStableGenes) < 100) {
		porder <- order(combinedPvalues, decreasing=TRUE)
		allStableGenes <- porder[1:100]
	}
	#Slope of th elowest 25% of the genes
	return (x[allStableGenes,])
}

#Normalization for batch effect.
BatchEffectLinnorm1 <- function(x, minZeroPortion, BE_F_LC_Genes = 0.25,BE_F_HC_Genes = 0.05, BE_F_p = 0.5, BE_strength = 0.25, spikein = NULL) {
	x2 <- x
	MeanSDSkew2 <- NZrowLogMeanSDSkew(x)
	MeanSDSkew <- MeanSDSkew2
	
	MeanOrder <- order(MeanSDSkew[1,], decreasing = FALSE)
	x <- x[MeanOrder,]
	MeanSDSkew <- MeanSDSkew[,MeanOrder]
	MeanSDSkew2 <- MeanSDSkew2[,MeanOrder]
	
	
	a <- which(!is.nan(MeanSDSkew[3,]))
	x <- x[a,]
	MeanSDSkew <- MeanSDSkew[,a]
	
	#Sort data and filter LowGeneFil and minZeroPortion
	if (minZeroPortion == 0) {
		Keep <- which(rowSums(x != 0) >= ncol(x) * minZeroPortion)
		MeanSDSkew <- MeanSDSkew[,Keep]
		x <- x[Keep,]
	} else {
		Keep <- which(rowSums(x != 0) > ncol(x) * minZeroPortion)
		MeanSDSkew <- MeanSDSkew[,Keep]
		x <- x[Keep,]
	}
	
	Start <- floor(nrow(x) * BE_F_LC_Genes + 1)
	End <- nrow(x) - floor(nrow(x) * BE_F_HC_Genes)
	
	Keep <- Start:End
	x <- x[Keep,]
	MeanSDSkew <- MeanSDSkew[,Keep]
	
	
	
	allStableGenes <- 0
	logitit <- loessFit(MeanSDSkew[2,],MeanSDSkew[1,], weights=1/MeanSDSkew[2,]^2)
	Keep <- which(logitit$fitted > 0)
	LogFit <- logitit$fitted[Keep]
	x <- x[Keep,]
	MeanSDSkew <- MeanSDSkew[,Keep]
	
	SDRatio <- log(as.numeric(MeanSDSkew[2,]/LogFit))
	#normalize SDRatio
	LR <- LinearRegression(MeanSDSkew[1,],SDRatio)
	Residual <- SDRatio - (LR$coefficients[[2]] * MeanSDSkew[1,] + LR$coefficients[[1]])
	LR2 <- LinearRegression(MeanSDSkew[1,],abs(Residual))
	SDRatio <- SDRatio * (LR2$coefficients[[2]] * MeanSDSkew[1,1] + LR2$coefficients[[1]])/(LR2$coefficients[[2]] * MeanSDSkew[1,] + LR2$coefficients[[1]])
	
	
	
	pvalueMatrix <- matrix(nrow=ncol(MeanSDSkew), ncol=2)
	
	if (length(spikein) < 10) {
		SDnoOutlier <- SDRatio[!SDRatio %in% boxplot.stats(SDRatio)$out]
		Themean <- mean(SDnoOutlier,na.rm=TRUE)
		TheSD <- sd(SDnoOutlier,na.rm=TRUE)
		Bigger <- which(SDRatio >= Themean)
		smaller <- which(SDRatio < Themean)
		pvalues <- pnorm(SDRatio,Themean,TheSD, lower.tail = FALSE)
		pvalueMatrix[Bigger,1] <- 2 * pvalues[Bigger]
		pvalueMatrix[smaller,1] <- 2 * (1 - pvalues[smaller])
		#pvalueMatrix[,1] <- pnorm(SDRatio,mean(SDRatio,na.rm=TRUE),sd(SDRatio,na.rm=TRUE), lower.tail = FALSE )
		
		SkewLR <- LinearRegression(MeanSDSkew[1,],MeanSDSkew[3,])
		SkewResidual <- MeanSDSkew[3,] - (SkewLR$coefficients[[2]] * MeanSDSkew[3,] - SkewLR$coefficients[[1]])
		SkewnoOutlier <- SkewResidual[!SkewResidual %in% boxplot.stats(SkewResidual)$out]
		
		Themean <- mean(SkewnoOutlier,na.rm=TRUE)
		TheSD <- sd(SkewnoOutlier,na.rm=TRUE)
		Bigger <- which(SkewResidual >= Themean)
		smaller <- which(SkewResidual < Themean)
		pvalues <- pnorm(SkewResidual,Themean,TheSD, lower.tail = FALSE)
		pvalueMatrix[Bigger,2] <- 2 * pvalues[Bigger]
		pvalueMatrix[smaller,2] <- 2 * (1 - pvalues[smaller])
	} else {
		spikes <- which(rownames(x) %in% spikein)
		SDnoOutlier <- SDRatio[spikes]
		Themean <- mean(SDnoOutlier,na.rm=TRUE)
		TheSD <- sd(SDnoOutlier,na.rm=TRUE)
		Bigger <- which(SDRatio >= Themean)
		smaller <- which(SDRatio < Themean)
		pvalues <- pnorm(SDRatio,Themean,TheSD, lower.tail = FALSE)
		pvalueMatrix[Bigger,1] <- 2 * pvalues[Bigger]
		pvalueMatrix[smaller,1] <- 2 * (1 - pvalues[smaller])
		#pvalueMatrix[,1] <- pnorm(SDRatio,mean(SDRatio,na.rm=TRUE),sd(SDRatio,na.rm=TRUE), lower.tail = FALSE )
		
		SkewLR <- LinearRegression(MeanSDSkew[1,],MeanSDSkew[3,])
		SkewResidual <- MeanSDSkew[3,] - (SkewLR$coefficients[[2]] * MeanSDSkew[3,] - SkewLR$coefficients[[1]])
		SkewnoOutlier <- SkewResidual[spikes]
		
		Themean <- mean(SkewnoOutlier,na.rm=TRUE)
		TheSD <- sd(SkewnoOutlier,na.rm=TRUE)
		Bigger <- which(SkewResidual >= Themean)
		smaller <- which(SkewResidual < Themean)
		pvalues <- pnorm(SkewResidual,Themean,TheSD, lower.tail = FALSE)
		pvalueMatrix[Bigger,2] <- 2 * pvalues[Bigger]
		pvalueMatrix[smaller,2] <- 2 * (1 - pvalues[smaller])
	}
	
	#combinedPvalues <- apply(pvalueMatrix,1,FisherMethod)
	combinedPvalues <- empiricalBrownsMethod(MeanSDSkew[2:3,],pvalueMatrix)

	combinedPvalues[is.na(combinedPvalues)] <- 0

	allStableGenes <- which(combinedPvalues > BE_F_p)

	#Safety, need 100 genes at least
	while (length(allStableGenes) < 100) {
		porder <- order(combinedPvalues, decreasing=TRUE)
		allStableGenes <- porder[1:100]
	}
	
	Results <- BatchEffect2(x[allStableGenes,], x2, MeanSDSkew[1,allStableGenes], BE_strength)
	colnames(Results) <- colnames(x2)
	rownames(Results) <- rownames(x2)
	#x2[is.infinite(x2)] <- 0
	return (Results)
}
#Batch effect results
BatchEffect2 <- function(x,y,z,z2) {
	.Call(BatchEffectCpp, x,y,z,z2)
}

#Linnorm's main funciton. Find optimal Lambda. This is implemented in C++.
LocateLambda <- function(x,y,z) {
    .Call(LocateLambdaCpp, x,y,z)
}
SkewVar <- function(x,y) {
    .Call(SkewVarCpp, x,y)
}
SkewAVar <- function(x,y) {
    .Call(SkewAVarCpp, x,y)
}

#Get Slope from x and y vectors
getSlope <- function(x,y) {
	.Call(getSlopeCpp, x,y)
}

#Create index for parsing correlation matrix on the "upper triangle" vector.
createUpperIndex <- function(colLength,TotalLength) {
	index <- matrix(0,nrow=TotalLength,ncol=2)
	TL <- 1
	for (i in 2:colLength) {
		newmatrix <- matrix(i, ncol=2, nrow=(i-1))
		newmatrix[,1] <- seq(1,(i-1),1)
		index[TL:(TL+i-2),] <- newmatrix
		TL <- TL + i - 1
	}
	return (index)
}
#for "lower triangle"
createLowerIndex <- function(rowLength,TotalLength) {
	index <- matrix(0,nrow=TotalLength,ncol=2)
	TL <- 1
	for (i in 2:rowLength) {
		newmatrix <- matrix(i, ncol=2, nrow=(i-1))
		newmatrix[,2] <- seq(1,(i-1),1)
		index[TL:(TL+i-2),] <- newmatrix
		TL <- TL + i - 1
	}
	return (index)
}
#Convert  "upper triangle" vector to matrix.
UpperToMatrix <- function(datavalues,UpperIndex) {
	theMatrix <- matrix(1, ncol=max(UpperIndex[,2]), nrow=max(UpperIndex[,2]))
	upperrow <- order(UpperIndex[,1])
	lowerrow <- order(UpperIndex[,2])
	upperrowi <- 1
	upperrowj <- 1
	lowerrowi <- 1
	lowerrowj <- 1
	for (i in 1:ncol(theMatrix)) {
		if (upperrowj < length(upperrow)) {
			while (UpperIndex[upperrow[upperrowj],1] == i) {
				upperrowj <- upperrowj + 1
				if (upperrowj == length(upperrow)) {
					upperrowj <- upperrowj + 1
					break
				}
			}
			theMatrix[i,UpperIndex[upperrow[upperrowi:(upperrowj-1)],2]] <- datavalues[upperrow[upperrowi:(upperrowj-1)]]
		}
		if (lowerrowj < length(lowerrow)) {
			while (UpperIndex[lowerrow[lowerrowj],2] == i) {
				lowerrowj <- lowerrowj + 1
				if (lowerrowj == length(lowerrow)) {
					lowerrowj <- lowerrowj + 1
					break
				}
			}
			theMatrix[i,UpperIndex[lowerrow[lowerrowi:(lowerrowj-1)],1]] <- datavalues[lowerrow[lowerrowi:(lowerrowj-1)]]
		}
		upperrowi <- upperrowj
		lowerrowi <- lowerrowj
	}
	return (theMatrix)
}

#Check if input are colors
areColors <- function(x) {
	sapply(x, function(X) {
		tryCatch(is.matrix(col2rgb(X)), 
		error = function(e) FALSE)
	})
}

#RangeRatio for HVG discovery
RangeRatio <- function(x) {
	Keep <- sort(x[x!=0])
	if (length(Keep) <= 2) {
		return(NA)
	}
	twentypercent <- floor(length(Keep) * 0.15) + 1
	Top5 <- mean(Keep[(length(Keep)-twentypercent -1):length(Keep)])
	Bottom5 <- mean(Keep[1:twentypercent])
	answer <- Top5 - Bottom5
	return(answer) 
}
