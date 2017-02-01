#Calculate Moments:

rowVars <- function(x) {
    .Call(colVarsCpp, t(x))
}
rowSDs <- function(x) {
    .Call(colSDsCpp, t(x))
}
colVars <- function(x) {
    .Call(colVarsCpp, x)
}
colSDs <- function(x) {
    .Call(colSDsCpp, x)
}
colMeanSD <- function(x) {
    .Call(colMeanSDCpp, x)
}
rowMeanSD <- function(x) {
    .Call(colMeanSDCpp, t(x))
}

NZrowMeans <- function(x) {
    .Call(NZcolMeansCpp, t(x))
}

NZcolMeans <- function(x) {
    .Call(NZcolMeansCpp, x)
}
NZcolMeanSD <- function(x) {
    .Call(NZcolMeanSDCpp, x)
}

colLog1pMeanSD <- function(x,y) {
    .Call(colLog1pMeanSDCpp, x,y)
}
rowLog1pMeanSD <- function(x,y) {
    .Call(colLog1pMeanSDCpp, t(x),y)
}

NZcolLogMeanSDSkew <- function(x) {
    .Call(NZcolLogMeanSDSkewCpp, x)
}
NZrowLogMeanSDSkew  <- function(x) {
    .Call(NZcolLogMeanSDSkewCpp, t(x))
}

#Weighted rowMeans
WNZcolMeans <- function(x,y) {
	.Call(WNZcolMeansCpp, x,y)
}
WcolMeans <- function(x,y) {
	.Call(WcolMeansCpp, x,y)
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
