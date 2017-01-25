library(Linnorm)

context("Locate Lambda")
data(SEQC)
SEQC <- tXPM(SEQC)
SEQC <- SEQC[rowSums(SEQC != 0) > 19, ]
SEQC <- SEQC[order(rowMeans(SEQC)),]

LocateLambdaSEQC <- LocateLambda(SEQC, 100, 60000000)
LocateLambda_legacySEQC <- LocateLambda_legacy(SEQC, 100, 60000000)

test_that("LocateLambda is accurate", {
	expect_equal(LocateLambdaSEQC, LocateLambda_legacySEQC)
})
