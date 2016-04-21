rowVars <- function(x) {
  return ((rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1) )  )
}
rowSDs <- function(x) {
  return ((rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1) ) ^ 0.5 )
}


gammaShape <- function(x) {
    .Call(gammaShapeCpp, x)
}

LocateLambda <- function(x,y) {
    .Call(LocateLambdaCpp, x,y)
}


GammaSim <- function(thisdata, NumRep=3, NumDiff = 2000, NumFea = 20000, showinfo=FALSE, MaxLibSizelog2FC=0.5, DEGlog2FC="Auto") {
	thisdata_ori <- thisdata
	thisdata <- thisdata[rowSums(thisdata != 0) == length(thisdata[1,]),]

	#Capture distribution from the dataset for simulation
	#K is the shape parameter from Gamma distribution
	Fit <- lm((log(rowVars(thisdata)))~(log(rowMeans(thisdata))))
	FindVar <- function(inputmean) {
		return (exp(log(inputmean) * Fit$coefficients[[2]] + Fit$coefficients[[1]] ))
	}
	
	gammamatrix <- matrix(0, ncol=(2 * NumRep), nrow=NumFea)
	RN <- vector(mode="character", NumFea)
	for (i in 1:NumFea) {
		RN[i] <- paste("Gene",i,sep="")
	}
	rownames(gammamatrix) <- RN	
	
	
	#Obtain Significant fold change threshold
	Proportion <- NumDiff/NumFea
	if (Proportion > 0.9){
		stop("Error: NumDiff is too large. Proportion of Differential Feature is larger than 90%.")
	}
	if (Proportion <= 0){
		stop("Error: NumDiff is too small. Proportion of Differential Feature is smaller than 5%.")
	}
	#Minimum Fold change (FC) for pvalue to reach 0.05
	pvalue <- 1
	minBound <- 1
	#FC of 5 is sure to be significant, so there is no need to search for boundary. If it ever happens to be not so, it also serves as safety for the program, since FC > 5  not being significant doesn't make sense.
	maxBound <- 5
	midBound <- round((minBound + maxBound)/2,4)
	s <- median(rowMeans(thisdata_ori))
	theta1 <- FindVar(s)/s
	theK <-s/theta1
	
	NR2 <- NumRep
	design <- matrix(nrow=(NR2 * 2), ncol=2)
	colnames(design) <- c("SampleInfo", "Gam")
	col1 <- c()
	for (i in 1:NR2) {
		col1 <- c(col1, "Normal")
	}
	for (i in 1:NR2) {
		col1 <- c(col1, "Tumor")
	}
	design[,1] <- col1
	design <- as.data.frame(design)
	pvstore <- vector(mode="numeric",100)
	if (DEGlog2FC != "Auto") {
		SigFC <- log(2^DEGlog2FC)
		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find its stddev.
		minBound <- 0
		#FC of 100 is sure to be significant, so there is no need to search for boundary
		maxBound <- 100
		midBound <- (minBound + maxBound)/2
		Probability <- vector(mode="double",NumFea)
		while (midBound != minBound && midBound != maxBound) {
			normmodel <- rnorm(NumFea,0,midBound)
			for (i in seq_along(normmodel)) {
				Probability[i] <- 1 - (2 * (pnorm(abs(normmodel[i]),0,SigFC/2, lower.tail = TRUE) - 0.5 ) )
			}
			Probability[is.na(Probability)] <- 1
			answer1 <- Proportion - length(Probability[Probability < 0.05])/NumFea
			if (answer1 < 0) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		FCstddev <- midBound
		if (showinfo) {
			message("Significant Log 2 Fold Change for p value is ", log(exp(SigFC),2),".",appendLF=TRUE)
			message("Final Standard Deviation of log 2 fold change in the dataset is ", log(exp(FCstddev),2),".",appendLF=TRUE)
			flush.console()
		}
	} else { 
		while (midBound != minBound || midBound != maxBound) {
			s2 <- s * midBound
			theta2 <- FindVar(s2)/s2
			theK2 <- s2/theta2
			#Do the test 100 times for an average p value.
			for (i in 1:100) {
				a <- as.numeric(rgamma(NR2,shape=theK,scale=theta1))
				b <- as.numeric(rgamma(NR2,shape=theK2,scale=theta2))
				#WilTest <- wilcox.test(a,b,alternative = c("two.sided"))
				#pvstore[i] <- as.numeric(WilTest[3])
				expressionData <- c(log1p(a),log1p(b))
				design[,2] <- as.numeric(expressionData)
				AnovaSum <- summary(aov(Gam~SampleInfo,data=design))
				pvstore[i] <- AnovaSum[[1]][,5][1]
				if (is.na(pvstore[i])) {
					pvstore[i] <- 1
				}
			}
			if (mean(pvstore, na.rm=TRUE) < 0.05) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		SigFC <- log(midBound)
		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find its stddev.
		minBound <- 0
		#FC of 100 is sure to be significant, so there is no need to search for boundary
		maxBound <- 100
		midBound <- (minBound + maxBound)/2
		Probability <- vector(mode="double",NumFea)
		while (midBound != minBound && midBound != maxBound) {
			normmodel <- rnorm(NumFea,0,midBound)
			for (i in seq_along(normmodel)) {
				Probability[i] <- 1 - (2 * (pnorm(abs(normmodel[i]),0,SigFC/2, lower.tail = TRUE) - 0.5 ) )
			}
			Probability[is.na(Probability)] <- 1
			Probability <- as.numeric(p.adjust(Probability,method ="BH"))
			answer1 <- Proportion - length(Probability[Probability < 0.05])/NumFea
			if (answer1 < 0) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		FCstddev <- midBound
		if (showinfo) {
			message("Significant Log 2 Fold Change for p value is ", log(exp(SigFC),2),".",appendLF=TRUE)
			message("Final Standard Deviation of log 2 fold change in the dataset is ", log(exp(FCstddev),2),".",appendLF=TRUE)
			flush.console()
		}
	}
	#Define Fold Change of Genes.
	#Non differential features will have mean of 0 and stddev of SigFC
	changes <- sort(rnorm(NumFea,0,FCstddev))
	begin <- round(NumDiff/2,0)
	ending <- round(NumFea - NumDiff/2,0)
	unchanged <- changes[begin:ending]
	changing <- c(changes[1:begin],changes[ending:length(changes)])

	
	unchanged <- exp(unchanged)
	changing <- exp(changing)
	tobechanged <- sample(1:NumFea,NumDiff)
	for (i in 1:NumFea) {
		#sample from thisdata
		TDsample <- sample(1:length(thisdata_ori[,1]),1)
		thisrow <- thisdata_ori[TDsample,]
		dmean <- mean(thisrow)
		theta <- FindVar(dmean)/dmean
		theK <- dmean/theta
		if (dmean == 0) {
			gammamatrix[i,] <- rep(0,length(gammamatrix[1,]))
		} else {
			if ( i %in% tobechanged) {
				if (sample(1:2,1) == 1) {
					gammamatrix[i,1:NumRep] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
					TBC <- sample(changing,1)
					dmean <- dmean * TBC
					theta <- FindVar(dmean)/dmean
					theK <- dmean/theta
					gammamatrix[i,(NumRep + 1):length(gammamatrix[1,])] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
				} else {
					gammamatrix[i,(NumRep + 1):length(gammamatrix[1,])] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
					TBC <- sample(changing,1)
					dmean <- dmean * TBC
					theta <- FindVar(dmean)/dmean
					theK <- dmean/theta
					gammamatrix[i,1:NumRep] <- sample(as.vector(rgamma(NumRep * 2, shape=theK,scale=theta)),NumRep)
				}
			} else {
				if (sample(1:2,1) == 1) {
					gammamatrix[i,1:NumRep] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
					TBC <- sample(unchanged,1)
					dmean <- dmean * TBC
					theta <- FindVar(dmean)/dmean
					theK <- dmean/theta
					gammamatrix[i,(NumRep + 1):length(gammamatrix[1,])] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
				} else {
					gammamatrix[i,(NumRep + 1):length(gammamatrix[1,])] <- sample(as.vector(rgamma(NumRep* 2, shape=theK,scale=theta)),NumRep)
					TBC <- sample(unchanged,1)
					dmean <- dmean * TBC
					theta <- FindVar(dmean)/dmean
					theK <- dmean/theta
					gammamatrix[i,1:NumRep] <- sample(as.vector(rgamma(NumRep * 2,shape=theK,scale=theta)),NumRep)
				}
			}
		}
	}
	gammamatrix[is.na(gammamatrix)] = 0
	gammamatrix[is.infinite(gammamatrix)] = 0
	gammamatrix[gammamatrix < 0] = 0
	changes <- 2^(seq(-MaxLibSizelog2FC,0,0.001) )
	thesechanges <- sample(changes,length(gammamatrix[1,]),replace=TRUE)

	for (i in seq_along(gammamatrix[1,])) {
		gammamatrix[,i] <- gammamatrix[,i] * thesechanges[i]
	}

	gammamatrix <-  as.matrix(floor(gammamatrix))
	listing <- list(gammamatrix,as.numeric(sort(tobechanged)))
	results <- setNames(listing, c("data", "DiffList"))
	return (results)
}

PoissonSim <- function(thisdata, NumRep=3, NumDiff = 2000, NumFea = 20000, showinfo=FALSE, MaxLibSizelog2FC=0.5, DEGlog2FC="Auto") {
	thisdata_ori <- thisdata
	thisdata <- thisdata[rowSums(thisdata != 0) == length(thisdata[1,]),]
	poismatrix <- matrix(0, ncol=(2 * NumRep), nrow=NumFea)
	RN <- vector(mode="character", NumFea)
	for (i in 1:NumFea) {
		RN[i] <- paste("Gene",i,sep="")
	}
	rownames(poismatrix) <- RN	
	
	#Obtain Significant fold change threshold
	Proportion <- NumDiff/NumFea
	if (Proportion > 0.9){
		stop("Error: NumDiff is too large. Proportion of Differential Feature is larger than 90%.")
	}
	if (Proportion <= 0){
		stop("Error: NumDiff is too small. Proportion of Differential Feature is smaller than 5%.")
	}
	#Minimum FC for pvalue to reach 0.05
	pvalue <- 1
	minBound <- 1
	#FC of 5 is sure to be significant, so there is no need to search for boundary. If it ever happens to be not so, it also serves as safety for the program, since FC > 5  not being significant doesn't make sense.
	maxBound <- 5
	midBound <- round((minBound + maxBound)/2,4)
	s <- median(rowMeans(thisdata_ori))
	NR2 <- NumRep
	design <- matrix(nrow=(NR2 * 2), ncol=2)
	colnames(design) <- c("SampleInfo", "POIS")
	col1 <- c()
	for (i in 1:NR2) {
		col1 <- c(col1, "Normal")
	}
	for (i in 1:NR2) {
		col1 <- c(col1, "Tumor")
	}
	#####################################
	
	design[,1] <- col1
	design <- as.data.frame(design)
	pvstore <- vector(mode="numeric",100)
	
	if (DEGlog2FC != "Auto") {
		SigFC <- log(2^DEGlog2FC)
		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find its stddev.
		minBound <- 0
		#FC of 100 is sure to be significant, so there is no need to search for boundary
		maxBound <- 100
		midBound <- (minBound + maxBound)/2
		Probability <- vector(mode="double",NumFea)
		while (midBound != minBound && midBound != maxBound) {
			normmodel <- rnorm(NumFea,0,midBound)
			for (i in seq_along(normmodel)) {
				Probability[i] <- 1 - (2 * (pnorm(abs(normmodel[i]),0,SigFC/2, lower.tail = TRUE) - 0.5 ) )
			}
			Probability[is.na(Probability)] <- 1
			answer1 <- Proportion - length(Probability[Probability < 0.05])/NumFea
			if (answer1 < 0) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		FCstddev <- midBound
		if (showinfo) {
			message("Significant Log 2 Fold Change for p value is ", log(exp(SigFC),2),".",appendLF=TRUE)
			message("Final Standard Deviation of log 2 fold change in the dataset is ", log(exp(FCstddev),2),".",appendLF=TRUE)
			flush.console()
		}
	} else { 
		while (midBound != minBound || midBound != maxBound) {
			s2 <- (s*midBound)
			#Do the test 100 times for an average p value.
			for (i in 1:100) {
				expressionData <- c(log1p(as.numeric(rpois(NR2,s))),log1p(as.numeric(rpois(NR2,s2))))
				design[,2] <- as.numeric(expressionData)
				AnovaSum <- summary(aov(POIS~SampleInfo,data=design))
				pvstore[i] <- AnovaSum[[1]][,5][1]
				#WilTest <- wilcox.test(as.numeric(rpois(NR2,s)),as.numeric(rpois(NR2,s2)),alternative = c("two.sided"))
				#pvstore[i] <- as.numeric(WilTest[3])
				if (is.na(pvstore[i])) {
					pvstore[i] <- 1
				}
			}
			if (mean(pvstore) < 0.05) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		SigFC <- log(midBound)

		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find the stddev.
		minBound <- 0
		#FC of 100 is sure to be significant, so there is no need to search for boundary
		maxBound <- 100
		midBound <- (minBound + maxBound)/2
		Probability <- vector(mode="double",NumFea)
		while (midBound != minBound && midBound != maxBound) {
			normmodel <- rnorm(NumFea,0,midBound)
			for (i in seq_along(normmodel)) {
				Probability[i] <- 1 - (2 * (pnorm(abs(normmodel[i]),0,SigFC/2, lower.tail = TRUE) - 0.5 ) )
			}
			Probability[is.na(Probability)] <- 1
			Probability <- as.numeric(p.adjust(Probability,method ="BH"))
			answer1 <- Proportion - length(Probability[Probability < 0.05])/NumFea
			if (answer1 < 0) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		FCstddev <- midBound
		if (showinfo) {
			message("Significant Log 2 Fold Change for p value is ", log(exp(SigFC),2),".",appendLF=TRUE)
			message("Final Standard Deviation of log 2 fold change in the dataset is ", log(exp(FCstddev),2),".",appendLF=TRUE)
			flush.console()
		}
	}
	#Non differential features will have mean of 0 and stddev of SigFC
	changes <- sort(rnorm(NumFea,0,FCstddev))
	begin <- round(NumDiff/2,0)
	ending <- round(NumFea - NumDiff/2,0)
	unchanged <- changes[begin:ending]
	changing <- c(changes[1:begin],changes[ending:length(changes)])

	
	unchanged <- exp(unchanged)
	changing <- exp(changing)
	tobechanged <- sample(1:NumFea,NumDiff)
	for (i in 1:NumFea) {
		#sample from thisdata
		TDsample <- sample(1:length(thisdata_ori[,1]),1)
		dmean <- mean(as.numeric(thisdata_ori[TDsample,]))
		if (dmean == 0) {
			poismatrix[i,] <- rep(0,length(poismatrix[1,]))
		} else {
			if ( i %in% tobechanged) {
				if (sample(1:2,1) == 1) {
					poismatrix[i,1:NumRep] <- sample(as.vector(rpois(NumRep* 2, dmean)),NumRep)
					TBC <- sample(changing,1)
					dmean <- dmean * TBC
					poismatrix[i,(NumRep + 1):length(poismatrix[1,])] <- sample(as.vector(rpois(NumRep* 2, dmean)),NumRep)
				} else {
					poismatrix[i,(NumRep + 1):length(poismatrix[1,])] <- sample(as.vector(rpois(NumRep* 2, dmean)),NumRep)
					TBC <- sample(changing,1)
					dmean <- dmean * TBC
					poismatrix[i,1:NumRep] <- sample(as.vector(rpois(NumRep * 2, dmean)),NumRep)
				}
			} else {
				if (sample(1:2,1) == 1) {
					poismatrix[i,1:NumRep] <- sample(as.vector(rpois(NumRep* 2, dmean)),NumRep)
					TBC <- sample(unchanged,1)
					dmean <- dmean * TBC
					poismatrix[i,(NumRep + 1):length(poismatrix[1,])] <- sample(as.vector(rpois(NumRep* 2, dmean)),NumRep)
				} else {
					poismatrix[i,(NumRep + 1):length(poismatrix[1,])] <- sample(as.vector(rpois(NumRep* 2, dmean)),NumRep)
					TBC <- sample(unchanged,1)
					dmean <- dmean * TBC
					poismatrix[i,1:NumRep] <- sample(as.vector(rpois(NumRep * 2, dmean)),NumRep)
				}
			}
		}
	}
	poismatrix[is.na(poismatrix)] = 0
	poismatrix[is.infinite(poismatrix)] = 0
	poismatrix[poismatrix < 0] = 0
	changes <- 2^(seq(-MaxLibSizelog2FC,0,0.001) )
	thesechanges <- sample(changes,length(poismatrix[1,]),replace=TRUE)
	
	for (i in seq_along(poismatrix[1,])) {
		poismatrix[,i] <- poismatrix[,i] * thesechanges[i]
	}

	poismatrix <-  as.matrix(floor(poismatrix))
	listing <- list(poismatrix,as.numeric(sort(tobechanged)))
	results <- setNames(listing, c("data", "DiffList"))

	return (results)
}

LogNormSim <- function(thisdata, NumRep=3, NumDiff = 2000, NumFea = 20000, showinfo=FALSE, MaxLibSizelog2FC=0.5, DEGlog2FC="Auto") {
	thisdata_ori <- thisdata
	thisdata <- thisdata[rowSums(thisdata != 0) == length(thisdata[1,]),]
	#Capture distribution from the dataset for simulation
	Fit <- lm((log(rowSDs(thisdata)))~(log(rowMeans(thisdata))))
	FindSD <- function(inputmean) {
		return (exp(log(inputmean) * Fit$coefficients[[2]] + Fit$coefficients[[1]] ))
	}
	lnormmatrix <- matrix(0, ncol=(2 * NumRep), nrow=NumFea)
	RN <- vector(mode="character", NumFea)
	for (i in 1:NumFea) {
		RN[i] <- paste("Gene",i,sep="")
	}
	rownames(lnormmatrix) <- RN	
	
	
	#Obtain Significant fold change threshold
	Proportion <- NumDiff/NumFea
	if (Proportion > 0.9){
		stop("Error: NumDiff is too large. Proportion of Differential Feature is larger than 90%.")
	}
	if (Proportion <= 0){
		stop("Error: NumDiff is too small. Proportion of Differential Feature is smaller than 5%.")
	}
	#Minimum FC for pvalue to reach 0.05
	pvalue <- 1
	minBound <- 1
	#FC of 5 is sure to be significant, so there is no need to search for boundary. If it ever happens to be not so, it also serves as safety for the program, since FC > 5  not being significant doesn't make sense.
	maxBound <- 5
	midBound <- round((minBound + maxBound)/2,4)
	theMean <- median(rowMeans(thisdata_ori))
	theSD <- FindSD(theMean)
	LNMean <- log(theMean) - 0.5 * log((theSD/theMean)^2 + 1)
	LNSD <- (log((theSD/theMean)^2 + 1) )^0.5
	NR2 <- NumRep
	design <- matrix(nrow=(NR2 * 2), ncol=2)
	colnames(design) <- c("SampleInfo", "lnorm")
	col1 <- c()
	for (i in 1:NR2) {
		col1 <- c(col1, "Normal")
	}
	for (i in 1:NR2) {
		col1 <- c(col1, "Tumor")
	}
	design[,1] <- col1
	design <- as.data.frame(design)
	pvstore <- vector(mode="numeric",100)
	if (DEGlog2FC != "Auto") {
		SigFC <- log(2^DEGlog2FC)
		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find its stddev.
		minBound <- 0
		#FC of 100 is sure to be significant, so there is no need to search for boundary
		maxBound <- 100
		midBound <- (minBound + maxBound)/2
		Probability <- vector(mode="double",NumFea)
		while (midBound != minBound && midBound != maxBound) {
			normmodel <- rnorm(NumFea,0,midBound)
			for (i in seq_along(normmodel)) {
				Probability[i] <- 1 - (2 * (pnorm(abs(normmodel[i]),0,SigFC/2, lower.tail = TRUE) - 0.5 ) )
			}
			Probability[is.na(Probability)] <- 1
			answer1 <- Proportion - length(Probability[Probability < 0.05])/NumFea
			if (answer1 < 0) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		FCstddev <- midBound
		if (showinfo) {
			message("Significant Log 2 Fold Change for p value is ", log(exp(SigFC),2),".",appendLF=TRUE)
			message("Final Standard Deviation of log 2 fold change in the dataset is ", log(exp(FCstddev),2),".",appendLF=TRUE)
			flush.console()
		}
	} else { 
		while (midBound != minBound || midBound != maxBound) {
			theMean2 <- theMean * midBound
			theSD2 <- FindSD(theMean2)
			LNMean2 <- log(theMean2) - 0.5 * log((theSD2/theMean2)^2 + 1)
			LNSD2 <- (log((theSD2/theMean2)^2 + 1) )^0.5
			#Do the test 100 times for an average p value.
			for (i in 1:100) {
				a <- as.numeric(rlnorm(NR2,meanlog = LNMean, sdlog = LNSD))
				b <- as.numeric(rlnorm(NR2,meanlog = LNMean2, sdlog = LNSD2))
				a[a < 0] <- 0
				b[b < 0] <- 0
				#WilTest <- wilcox.test(a,b,alternative = c("two.sided"))
				#pvstore[i] <- as.numeric(WilTest[3])
				expressionData <- c(log1p(a),log1p(b))
				design[,2] <- as.numeric(expressionData)
				AnovaSum <- summary(aov(lnorm~SampleInfo,data=design))
				pvstore[i] <- AnovaSum[[1]][,5][1]
				if (is.na(pvstore[i])) {
					pvstore[i] <- 1
				}
			}
			if (mean(pvstore) < 0.05) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		SigFC <- log(midBound)

		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find the stddev.
		minBound <- 0
		#FC of 100 is sure to be significant, so there is no need to search for boundary
		maxBound <- 100
		midBound <- (minBound + maxBound)/2
		Probability <- vector(mode="double",NumFea)
		while (midBound != minBound && midBound != maxBound) {
			normmodel <- rnorm(NumFea,0,midBound)
			for (i in seq_along(normmodel)) {
				Probability[i] <- 1 - (2 * (pnorm(abs(normmodel[i]),0,SigFC/2, lower.tail = TRUE) - 0.5 ) )
			}
			Probability[is.na(Probability)] <- 1
			Probability <- as.numeric(p.adjust(Probability,method ="BH"))
			answer1 <- Proportion - length(Probability[Probability < 0.05])/NumFea
			if (answer1 < 0) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		FCstddev <- midBound
		if (showinfo) {
			message("Significant Log 2 Fold Change for p value is ", log(exp(SigFC),2),".",appendLF=TRUE)
			message("Final Standard Deviation of log 2 fold change in the dataset is ", log(exp(FCstddev),2),".",appendLF=TRUE)
			flush.console()
		}
	}
	#Non differential features will have mean of 0 and stddev of SigFC
	changes <- sort(rnorm(NumFea,0,FCstddev))
	begin <- round(NumDiff/2,0)
	ending <- round(NumFea - NumDiff/2,0)
	unchanged <- changes[begin:ending]
	changing <- c(changes[1:begin],changes[ending:length(changes)])

	unchanged <- exp(unchanged)
	changing <- exp(changing)
	tobechanged <- sample(1:NumFea,NumDiff)
	for (i in 1:NumFea) {
		#sample from thisdata
		TDsample <- sample(1:length(thisdata_ori[,1]),1)
		thisrow <- thisdata_ori[TDsample,]
		theMean <- mean(thisrow)
		theSD <- FindSD(theMean)
		LNMean <- log(theMean) - 0.5 * log((theSD/theMean)^2 + 1)
		LNSD <- (log((theSD/theMean)^2 + 1) )^0.5
		if (theMean == 0) {
			lnormmatrix[i,] <- rep(0,length(lnormmatrix[1,]))
		} else {
			if ( i %in% tobechanged) {
				if (sample(1:2,1) == 1) {
					lnormmatrix[i,1:NumRep] <- (sample(as.vector(rlnorm(NumRep* 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
					TBC <- sample(changing,1)
					theMean <- theMean * TBC
					theSD <- FindSD(theMean)
					LNMean <- log(theMean) - 0.5 * log((theSD/theMean)^2 + 1)
					LNSD <- (log((theSD/theMean)^2 + 1) )^0.5
					lnormmatrix[i,(NumRep + 1):length(lnormmatrix[1,])] <- (sample(as.vector(rlnorm(NumRep* 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
				} else {
					lnormmatrix[i,(NumRep + 1):length(lnormmatrix[1,])] <- (sample(as.vector(rlnorm(NumRep* 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
					TBC <- sample(changing,1)
					theMean <- theMean * TBC
					theSD <- FindSD(theMean)
					LNMean <- log(theMean) - 0.5 * log((theSD/theMean)^2 + 1)
					LNSD <- (log((theSD/theMean)^2 + 1) )^0.5
					lnormmatrix[i,1:NumRep] <- (sample(as.vector(rlnorm(NumRep * 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
				}
			} else {
				if (sample(1:2,1) == 1) {
					lnormmatrix[i,1:NumRep] <- (sample(as.vector(rlnorm(NumRep* 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
					TBC <- sample(unchanged,1)
					theMean <- theMean * TBC
					theSD <- FindSD(theMean)
					LNMean <- log(theMean) - 0.5 * log((theSD/theMean)^2 + 1)
					LNSD <- (log((theSD/theMean)^2 + 1) )^0.5
					lnormmatrix[i,(NumRep + 1):length(lnormmatrix[1,])] <- (sample(as.vector(rlnorm(NumRep* 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
				} else {
					lnormmatrix[i,(NumRep + 1):length(lnormmatrix[1,])] <- (sample(as.vector(rlnorm(NumRep* 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
					TBC <- sample(unchanged,1)
					theMean <- theMean * TBC
					theSD <- FindSD(theMean)
					LNMean <- log(theMean) - 0.5 * log((theSD/theMean)^2 + 1)
					LNSD <- (log((theSD/theMean)^2 + 1) )^0.5
					lnormmatrix[i,1:NumRep] <- (sample(as.vector(rlnorm(NumRep * 2,meanlog = LNMean, sdlog = LNSD)),NumRep))
				}
			}
		}
	}
	
	lnormmatrix[is.na(lnormmatrix)] = 0
	lnormmatrix[is.infinite(lnormmatrix)] = 0
	lnormmatrix[lnormmatrix < 0] = 0
	changes <- 2^(seq(-MaxLibSizelog2FC,0,0.001) )
	thesechanges <- sample(changes,length(lnormmatrix[1,]),replace=TRUE)
	for (i in seq_along(lnormmatrix[1,])) {
		lnormmatrix[,i] <- lnormmatrix[,i] * thesechanges[i]
	}
	lnormmatrix <-  as.matrix(floor(lnormmatrix))
	listing <- list(lnormmatrix,as.numeric(sort(tobechanged)))
	results <- setNames(listing, c("data", "DiffList"))
	return (results)
}

NBSim <- function(thisdata, NumRep=3, NumDiff = 2000, NumFea = 20000, showinfo=FALSE, MaxLibSizelog2FC=0.5, DEGlog2FC="Auto") {
	thisdata_ori <- thisdata
	thisdata <- thisdata[rowSums(thisdata != 0) == length(thisdata[1,]),]
	#Capture distribution from the dataset for simulation
	MeanList <- vector(mode="numeric",length(thisdata[,1]))
	d <- vector(mode="numeric",length(thisdata[,1]))
	
	for (i in seq_along(thisdata[,1])){
		x <- thisdata[i,]
		MeanList[i] <- mean(x)
		d[i] <- as.numeric(unlist((glm.nb(x~1))[[24]]))
	}
	MeanList <- MeanList[!is.na(d)]
	d <- d[!is.na(d)]

	Fit <- lm((log(d))~(log(MeanList)))
	FindDispersion <- function(inputmean) {
		return (exp(log(inputmean) * Fit$coefficients[[2]] + Fit$coefficients[[1]] ))
	}
	
	nbmatrix <- matrix(0, ncol=(2 * NumRep), nrow=NumFea)
	RN <- vector(mode="character", NumFea)
	for (i in 1:NumFea) {
		RN[i] <- paste("Gene",i,sep="")
	}
	rownames(nbmatrix) <- RN	
	
	
	#Obtain Significant fold change threshold
	Proportion <- NumDiff/NumFea
	if (Proportion > 0.9){
		stop("Error: NumDiff is too large. Proportion of Differential Feature is larger than 90%.")
	}
	if (Proportion <= 0){
		stop("Error: NumDiff is too small. Proportion of Differential Feature is smaller than 5%.")
	}
	#Minimum FC for pvalue to reach 0.05
	pvalue <- 1
	minBound <- 1
	#FC of 5 is sure to be significant, so there is no need to search for boundary. If it ever happens to be not so, it also serves as safety for the program, since FC > 5  not being significant doesn't make sense.
	maxBound <- 5
	midBound <- round((minBound + maxBound)/2,4)
	themean <- median(rowMeans(thisdata_ori))
	theDis <- FindDispersion(themean)
	
	NR2 <- NumRep
	design <- matrix(nrow=(NR2 * 2), ncol=2)
	colnames(design) <- c("SampleInfo", "Gam")
	col1 <- c()
	for (i in 1:NR2) {
		col1 <- c(col1, "Normal")
	}
	for (i in 1:NR2) {
		col1 <- c(col1, "Tumor")
	}
	design[,1] <- col1
	design <- as.data.frame(design)
	pvstore <- vector(mode="numeric",100)
	if (DEGlog2FC != "Auto") {
		SigFC <- log(2^DEGlog2FC)
		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find its stddev.
		minBound <- 0
		#FC of 100 is sure to be significant, so there is no need to search for boundary
		maxBound <- 100
		midBound <- (minBound + maxBound)/2
		Probability <- vector(mode="double",NumFea)
		while (midBound != minBound && midBound != maxBound) {
			normmodel <- rnorm(NumFea,0,midBound)
			for (i in seq_along(normmodel)) {
				Probability[i] <- 1 - (2 * (pnorm(abs(normmodel[i]),0,SigFC/2, lower.tail = TRUE) - 0.5 ) )
			}
			Probability[is.na(Probability)] <- 1
			answer1 <- Proportion - length(Probability[Probability < 0.05])/NumFea
			if (answer1 < 0) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		FCstddev <- midBound
		if (showinfo) {
			message("Significant Log 2 Fold Change for p value is ", log(exp(SigFC),2),".",appendLF=TRUE)
			message("Final Standard Deviation of log 2 fold change in the dataset is ", log(exp(FCstddev),2),".",appendLF=TRUE)
			flush.console()
		}
	} else { 
		while (midBound != minBound || midBound != maxBound) {
			themean2 <- themean * midBound
			theDis2 <- FindDispersion(themean2)
			
			#Do the test 100 times for an average p value.
			for (i in 1:100) {
				a <- as.numeric(rnbinom(NR2,size=theDis,mu=themean))
				b <- as.numeric(rnbinom(NR2,size=theDis2,mu=themean2))
				#WilTest <- wilcox.test(a,b,alternative = c("two.sided"))
				#pvstore[i] <- as.numeric(WilTest[3])
				expressionData <- c(log1p(a),log1p(b))
				design[,2] <- as.numeric(expressionData)
				AnovaSum <- summary(aov(Gam~SampleInfo,data=design))
				pvstore[i] <- AnovaSum[[1]][,5][1]
				if (is.na(pvstore[i])) {
					pvstore[i] <- 1
				}
			}
			if (mean(pvstore, na.rm=TRUE) < 0.05) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		SigFC <- log(midBound)

		#To create Normal Distributed of FCs, where portion of Differential features are satisfied, we need to find the stddev.
		minBound <- 0
		#FC of 100 is sure to be significant, so there is no need to search for boundary
		maxBound <- 100
		midBound <- (minBound + maxBound)/2
		Probability <- vector(mode="double",NumFea)
		while (midBound != minBound && midBound != maxBound) {
			normmodel <- rnorm(NumFea,0,midBound)
			for (i in seq_along(normmodel)) {
				Probability[i] <- 1 - (2 * (pnorm(abs(normmodel[i]),0,SigFC/2, lower.tail = TRUE) - 0.5 ) )
			}
			Probability[is.na(Probability)] <- 1
			Probability <- as.numeric(p.adjust(Probability,method ="BH"))
			answer1 <- Proportion - length(Probability[Probability < 0.05])/NumFea
			if (answer1 < 0) {
				maxBound <- midBound
			} else {
				minBound <- midBound
			}
			midBound <- (minBound + maxBound)/2
		}
		FCstddev <- midBound
		if (showinfo) {
			message("Significant Log 2 Fold Change for p value is ", log(exp(SigFC),2),".",appendLF=TRUE)
			message("Final Standard Deviation of log 2 fold change in the dataset is ", log(exp(FCstddev),2),".",appendLF=TRUE)
			flush.console()
		}
	}
	#Non differential features will have mean of 0 and stddev of SigFC
	changes <- sort(rnorm(NumFea,0,FCstddev))
	begin <- round(NumDiff/2,0)
	ending <- round(NumFea - NumDiff/2,0)
	unchanged <- changes[begin:ending]
	changing <- c(changes[1:begin],changes[ending:length(changes)])

	
	unchanged <- exp(unchanged)
	changing <- exp(changing)
	tobechanged <- sample(1:NumFea,NumDiff)
	for (i in 1:NumFea) {
		#sample from thisdata
		TDsample <- sample(1:length(thisdata_ori[,1]),1)
		thisrow <- thisdata_ori[TDsample,]
		themean <- mean(thisrow)
		theDis <- FindDispersion(themean)
	
		
		if (themean == 0) {
			nbmatrix[i,] <- rep(0,length(nbmatrix[1,]))
		} else {
			if ( i %in% tobechanged) {
				if (sample(1:2,1) == 1) {
					nbmatrix[i,1:NumRep] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
					TBC <- sample(changing,1)
					themean <- themean * TBC
					theDis <- FindDispersion(themean)
					nbmatrix[i,(NumRep + 1):length(nbmatrix[1,])] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
				} else {
					nbmatrix[i,(NumRep + 1):length(nbmatrix[1,])] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
					TBC <- sample(changing,1)
					themean <- themean * TBC
					theDis <- FindDispersion(themean)
					nbmatrix[i,1:NumRep] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
				}
			} else {
				if (sample(1:2,1) == 1) {
					nbmatrix[i,1:NumRep] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
					TBC <- sample(unchanged,1)
					themean <- themean * TBC
					theDis <- FindDispersion(themean)
					nbmatrix[i,(NumRep + 1):length(nbmatrix[1,])] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
				} else {
					nbmatrix[i,(NumRep + 1):length(nbmatrix[1,])] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
					TBC <- sample(unchanged,1)
					themean <- themean * TBC
					theDis <- FindDispersion(themean)
					nbmatrix[i,1:NumRep] <- sample(as.vector(rnbinom(NumRep* 2,size=theDis,mu=themean)),NumRep)
				}
			}
		}
	}
	nbmatrix[is.na(nbmatrix)] = 0
	nbmatrix[is.infinite(nbmatrix)] = 0
	nbmatrix[nbmatrix < 0] = 0
	changes <- 2^(seq(-MaxLibSizelog2FC,0,0.001) )
	thesechanges <- sample(changes,length(nbmatrix[1,]),replace=TRUE)

	for (i in seq_along(nbmatrix[1,])) {
		nbmatrix[,i] <- nbmatrix[,i] * thesechanges[i]
	}

	nbmatrix <-  as.matrix(floor(nbmatrix))
	listing <- list(nbmatrix,as.numeric(sort(tobechanged)))
	results <- setNames(listing, c("data", "DiffList"))

	return (results)
}


