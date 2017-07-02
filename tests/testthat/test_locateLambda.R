library(Linnorm)

context("Locate Lambda")
data(SEQC)
SEQC <- tXPM(SEQC)
SEQC <- SEQC[rowSums(SEQC != 0) > 4, ]
SEQC <- SEQC[order(rowMeans(SEQC)),]

LocateLambdaSEQC <- LocateLambda(SEQC, 100, 6000000)
LocateLambda_legacySEQC <- LocateLambda_legacy(SEQC, 100, 6000000)

PDSEQC <- abs((LocateLambdaSEQC - LocateLambda_legacySEQC)/mean(c(LocateLambda_legacySEQC,LocateLambdaSEQC))) * 100
AcceptableSEQC <- PDSEQC < 5

test_that("LocateLambda SEQC", {
	expect_true(AcceptableSEQC)
})
