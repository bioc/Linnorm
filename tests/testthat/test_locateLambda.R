library(Linnorm)

context("Locate Lambda")
data(SEQC)
SEQC <- tXPM(SEQC)
SEQC <- SEQC[rowSums(SEQC != 0) > 4, ]
SEQC <- SEQC[order(rowMeans(SEQC)),]

LocateLambdaSEQC <- LocateLambda(SEQC, 100, 60000000)
LocateLambda_legacySEQC <- LocateLambda_legacy(SEQC, 100, 60000000)

PDSEQC <- abs((LocateLambdaSEQC - LocateLambda_legacySEQC)/mean(c(LocateLambda_legacySEQC,LocateLambdaSEQC))) * 100
AcceptableSEQC <- PDSEQC < 5

test_that("LocateLambda SEQC", {
	expect_true(AcceptableSEQC)
})

data(LIHC)
LIHC <- LIHC/1000000
LIHC <- LIHC[rowSums(LIHC != 0) > 4, ]
LIHC <- LIHC[order(rowMeans(LIHC)),]

LocateLambdaLIHC <- LocateLambda(LIHC, 100, 60000000)
LocateLambda_legacyLIHC <- LocateLambda_legacy(LIHC, 100, 60000000)

PDLIHC <- abs((LocateLambdaLIHC - LocateLambda_legacyLIHC)/mean(c(LocateLambda_legacyLIHC,LocateLambdaLIHC))) * 100
AcceptableLIHC <- PDLIHC < 5

test_that("LocateLambda LIHC", {
	expect_true(AcceptableLIHC)
})
