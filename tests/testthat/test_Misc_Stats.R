library(Linnorm)
library(moments)
library(matrixStats)

context("Moment calculation check")

#Initialize datasets
data(LIHC)
LIHC <- LIHC/1000000
LIHC <- LIHC[rowSums(LIHC != 0) > 3, ]
LIHC <- LIHC[order(rowMeans(LIHC)),]

data(SEQC)
SEQC <- XPM(as.matrix(SEQC))
SEQC <- SEQC[rowSums(SEQC != 0) > 3, ]
SEQC <- SEQC[order(rowMeans(SEQC)),]

#rowSDs
matrixStatsAnswerLIHC <- rowSds(LIHC)
matrixStatsAnswerSEQC <- rowSds(SEQC)

LinnormAnswerLIHC <- rowSDs(LIHC)
LinnormAnswerSEQC <- rowSDs(SEQC)

test_that("rowSDs is accurate", {
	expect_equal(matrixStatsAnswerLIHC[14069],LinnormAnswerLIHC[14069])
	expect_equal(matrixStatsAnswerLIHC[15810],LinnormAnswerLIHC[15810])
	expect_equal(matrixStatsAnswerLIHC[12988],LinnormAnswerLIHC[12988])
	
	expect_equal(matrixStatsAnswerSEQC[13061],LinnormAnswerSEQC[13061])
	expect_equal(matrixStatsAnswerSEQC[16106],LinnormAnswerSEQC[16106])
	expect_equal(matrixStatsAnswerSEQC[10132],LinnormAnswerSEQC[10132])
})

#NZcolMeans
SEQC <- t(SEQC)
LIHC <- t(LIHC)
NZcolMeans2 <- function(x) {
	answer <- rep(0, ncol(x))
	for (i in 1:ncol(x)) {
		answer[i] <- mean(x[x[,i] != 0,i])
	}
	return(answer)
}
AnswerLIHC <- NZcolMeans2(LIHC)
AnswerSEQC <- NZcolMeans2(SEQC)

LinnormAnswerLIHC <- NZcolMeans(LIHC)
LinnormAnswerSEQC <- NZcolMeans(SEQC)

test_that("NZcolMeans is accurate", {
	expect_equal(AnswerLIHC[14069],LinnormAnswerLIHC[14069])
	expect_equal(AnswerLIHC[15810],LinnormAnswerLIHC[15810])
	expect_equal(AnswerLIHC[12988],LinnormAnswerLIHC[12988])
	
	expect_equal(AnswerSEQC[13061],LinnormAnswerSEQC[13061])
	expect_equal(AnswerSEQC[16106],LinnormAnswerSEQC[16106])
	expect_equal(AnswerSEQC[10132],LinnormAnswerSEQC[10132])
})

#NZcolLogMeanSDSkew
NZcolLogMeanSDSkew2 <- function(x) {
	answer <- matrix(nrow=3, ncol=ncol(x))
	for (i in 1:ncol(x)) {
		thisdata <- log(x[x[,i] != 0,i])
		answer[1,i] <- mean(thisdata)
		answer[2,i] <- sd(thisdata)
		answer[3,i] <- skewness(thisdata)
	}
	return(answer)
}

AnswerLIHC <- NZcolLogMeanSDSkew2(LIHC)
AnswerSEQC <- NZcolLogMeanSDSkew2(SEQC)

LinnormAnswerLIHC <- NZcolLogMeanSDSkew(LIHC)
LinnormAnswerSEQC <- NZcolLogMeanSDSkew(SEQC)

test_that("NZcolLogMeanSDSkew mean is accurate", {
	expect_equal(AnswerLIHC[1,14069],LinnormAnswerLIHC[1,14069])
	expect_equal(AnswerLIHC[1,15810],LinnormAnswerLIHC[1,15810])
	expect_equal(AnswerLIHC[1,12988],LinnormAnswerLIHC[1,12988])
	
	expect_equal(AnswerSEQC[1,13061],LinnormAnswerSEQC[1,13061])
	expect_equal(AnswerSEQC[1,16106],LinnormAnswerSEQC[1,16106])
	expect_equal(AnswerSEQC[1,10132],LinnormAnswerSEQC[1,10132])
})

test_that("NZcolLogMeanSDSkew sd is accurate", {
	expect_equal(AnswerLIHC[2,14069],LinnormAnswerLIHC[2,14069])
	expect_equal(AnswerLIHC[2,15810],LinnormAnswerLIHC[2,15810])
	expect_equal(AnswerLIHC[2,12988],LinnormAnswerLIHC[2,12988])
	
	expect_equal(AnswerSEQC[2,13061],LinnormAnswerSEQC[2,13061])
	expect_equal(AnswerSEQC[2,16106],LinnormAnswerSEQC[2,16106])
	expect_equal(AnswerSEQC[2,10132],LinnormAnswerSEQC[2,10132])
})

test_that("NZcolLogMeanSDSkew skew is accurate", {
	expect_equal(AnswerLIHC[3,14069],LinnormAnswerLIHC[3,14069])
	expect_equal(AnswerLIHC[3,15810],LinnormAnswerLIHC[3,15810])
	expect_equal(AnswerLIHC[3,12988],LinnormAnswerLIHC[3,12988])
	
	expect_equal(AnswerSEQC[3,13061],LinnormAnswerSEQC[3,13061])
	expect_equal(AnswerSEQC[3,16106],LinnormAnswerSEQC[3,16106])
	expect_equal(AnswerSEQC[3,10132],LinnormAnswerSEQC[3,10132])
})