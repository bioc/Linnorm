#Calculate Moments:

NZrowLog1pMeanSD <- function(x,y) {
    .Call(NZrowLog1pMeanSDCpp, x,y)
}
HZrowLog1pMeanSD <- function(x,y) {
    .Call(HZrowLog1pMeanSDCpp, x,y)
}

NZrowLog1pMeanSqrtSD <- function(x,y) {
    .Call(NZrowLog1pMeanSqrtSDCpp, x,y)
}
NZrowLogMeanSDSkew <- function(x) {
    .Call(NZrowLogMeanSDSkewCpp, x)
}
NZrowMeanSDSkew <- function(x) {
    .Call(NZrowMeanSDSkewCpp, x)
}
NZrowMeanSD <- function(x) {
    .Call(NZrowMeanSDCpp, x)
}

rowLog1pMeanSD <- function(x,y) {
    .Call(rowLog1pMeanSDCpp, x,y)
}

rowLogMeanSDSkew <- function(x) {
    .Call(rowLogMeanSDSkewCpp, x)
}
rowMeanSDSkew <- function(x) {
    .Call(rowMeanSDSkewCpp, x)
}
rowMeanSD <- function(x) {
    .Call(rowMeanSDCpp, x)
}

NZrowMeans <- function(x) {
    .Call(NZrowMeansCpp, x)
}
rowVars <- function(x) {
    .Call(rowVarsCpp, x)
}
rowSDs <- function(x) {
    .Call(rowSDsCpp, x)
}
colSDs <- function(x) {
    .Call(colSDsCpp, x)
}

#Weighted rowMeans
WrowMeans <- function(x,y) {
	.Call(WrowMeansCpp, x,y)
}
WNZrowMeans <- function(x,y) {
	.Call(WNZrowMeansCpp, x,y)
}


#Fisher's method of combining p values.
FishersMethod <- function(p) {
	pchisq( -2*sum(log(p)), df=length(p)*2, lower.tail=FALSE)
}

#Significance of Pearson correlation coefficient.
r.sig <- function(r,n) {
	tvalue <- abs(r) * sqrt((n - 2)/(1 - r^2))
	return(2*pt(tvalue, n, lower.tail =FALSE))
}
