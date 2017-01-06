/*
The MIT License (MIT) 
Copyright (c) <2016> <Shun Hang Yip>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/
#include <R_ext/Rdynload.h>
#include <R.h>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <RcppArmadillo.h>
#include <vector>
#include <limits>
using namespace std;
using namespace Rcpp;

//Convert dataset into XPM, each row is a gene(feature) each column is a sample.
SEXP XPMCpp (SEXP xSEXP){
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//Two pass
	//arma::vec ColumnSums = sum(GeneExp); This comes with error. commented.
	arma::vec ColumnSums(GeneExp.n_cols);
	ColumnSums.fill(0);
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		for (int j = 0; j < int(GeneExp.n_cols); j ++) {
			ColumnSums.at(j) += GeneExp.at(i,j);
		}
	}
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		for (int j = 0; j < int(GeneExp.n_cols); j ++) {
			GeneExp.at(i,j) = GeneExp.at(i,j)/ColumnSums.at(j);
		}
	}
	return Rcpp::wrap(GeneExp);
}
//The following three functions calculate the optimal lambda for the dataset.
//This function calculates F(lambda) (see article) based on the expression matrix and lambda.
double VarOnly (const arma::mat& GeneExp, const double& lambda2) {
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	//arma::vec skewvec(GeneExp.n_rows);
	arma::vec SDevvec(GeneExp.n_rows);
	//arma::vec kurtvec(GeneExp.n_rows);
	arma::vec meanvec(GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	//double SumKurt = 0;
	//double SumSkew = 0;
	double SumSD = 0;
	double SumMean = 0;
	//double SumKurtMean = 0;
	//double SumSkewMean = 0;
	double SumSDMean = 0;
	double SumMeanSq = 0;
	
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		//double M3 = 0;
		//double M4 = 0;
		//double delta_n2;
		double delta, delta_n, term1;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (GeneExp.at(i,n) == 0) {
				continue;
			}
			delta = log1p(GeneExp.at(i,n) * lambda2) - mean;
			delta_n = delta / (numData + 1);
			//delta_n2 = delta_n * delta_n;
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			//M4 = M4 + term1 * delta_n2 * (pow(numData + 1,2) - 3 * (numData + 1) + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
			//M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			numData++;
		}
		if (numData < 3) {
			continue;
		}
		//Here, calculate kurtosis, skewness and SD using 4th, 3rd and 2nd moments.
		//kurtvec.at(i) =  (GeneExp.n_cols * M4)/pow(M2,2) - 3;
		//skewvec.at(i) =  (sqrt(GeneExp.n_cols) * M3) / pow(M2,1.5);
		SDevvec.at(i) =  sqrt(M2/(GeneExp.n_cols - 1));
		
		//Linear regression
		//meanvec.at(i) = 1 + (mean - minMean)/range;
		//meanvec.at(i) = 1 + i;
		meanvec.at(i) = mean;
		//SumKurt += kurtvec.at(i);
		//SumSkew += skewvec.at(i);
		SumSD += SDevvec.at(i);
		SumMean += meanvec.at(i);
		//SumKurtMean += kurtvec.at(i) * meanvec.at(i);
		//SumSkewMean += skewvec.at(i) * meanvec.at(i);
		SumSDMean += SDevvec.at(i) * meanvec.at(i);
		SumMeanSq += pow(meanvec.at(i),2);
	}
	//y = Mx + C, here are the M and C results of the linear regression
	//double kurtM = (GeneExp.n_rows * SumKurtMean - SumKurt * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) ) ;
	//double kurtC = (SumKurt - kurtM * SumMean)/GeneExp.n_rows;
	//double skewM = (GeneExp.n_rows * SumSkewMean - SumSkew * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) ) ;
	//double skewC = (SumSkew - skewM * SumMean)/GeneExp.n_rows;
	double SDevM = (GeneExp.n_rows * SumSDMean - SumSD * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) );
	/*
	//integral of the linear equation of skewness
	double Skewzerointercept = -skewC/skewM;
	double Skewintegral = 0;
	skewM = skewM/2;
	if (Skewzerointercept > meanvec.at(0) && Skewzerointercept < meanvec.at(skewvec.n_elem-1)) {
		Skewintegral = (abs(skewM * (meanvec.at(skewvec.n_elem-1) + Skewzerointercept) + skewC) * (meanvec.at(skewvec.n_elem-1) - Skewzerointercept) + abs(skewM * (meanvec.at(0) + Skewzerointercept) + skewC)* (Skewzerointercept - meanvec.at(0)))/(meanvec.at(skewvec.n_elem-1) - meanvec.at(0));
	} else {
		Skewintegral = abs(skewM * (meanvec.at(skewvec.n_elem-1) + meanvec.at(0)) + skewC);
	}
	
	//integral of the linear equation of kurtosis
	double Kurtzerointercept = -kurtC/kurtM;
	double Kurtintegral = 0;
	kurtM = kurtM/2;
	if (Kurtzerointercept > meanvec.at(0) && Kurtzerointercept < meanvec.at(kurtvec.n_elem-1)) {
		Kurtintegral = (abs(kurtM * (meanvec.at(kurtvec.n_elem-1) + Kurtzerointercept) + kurtC) * (meanvec.at(kurtvec.n_elem-1) - Kurtzerointercept) + abs(kurtM * (meanvec.at(0) + Kurtzerointercept) + kurtC)* (Kurtzerointercept - meanvec.at(0)) )/(meanvec.at(kurtvec.n_elem-1) - meanvec.at(0));
	} else {
		Kurtintegral = abs(kurtM * (meanvec.at(kurtvec.n_elem-1) + meanvec.at(0)) + kurtC);
	}
	*/
	//return 2 * abs(SDevM) + Skewintegral + Kurtintegral;
	return abs(SDevM);
}
double SkewOnly (const arma::mat& GeneExp, const double& lambda2) {
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::vec skewvec(GeneExp.n_rows);
	//arma::vec SDevvec(GeneExp.n_rows);
	//arma::vec kurtvec(GeneExp.n_rows);
	arma::vec meanvec(GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	//double SumKurt = 0;
	double SumSkew = 0;
	//double SumSD = 0;
	double SumMean = 0;
	//double SumKurtMean = 0;
	double SumSkewMean = 0;
	//double SumSDMean = 0;
	double SumMeanSq = 0;
	
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		double M3 = 0;
		//double M4 = 0;
		//double delta_n2;
		double delta, delta_n, term1;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (GeneExp.at(i,n) == 0) {
				continue;
			}
			delta = log1p(GeneExp.at(i,n) * lambda2) - mean;
			delta_n = delta / (numData + 1);
			//delta_n2 = delta_n * delta_n;
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			//M4 = M4 + term1 * delta_n2 * (pow(numData + 1,2) - 3 * (numData + 1) + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
			M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			numData++;
		}
		if (numData < 3) {
			continue;
		}
		//Here, calculate kurtosis, skewness and SD using 4th, 3rd and 2nd moments.
		//kurtvec.at(i) =  (GeneExp.n_cols * M4)/pow(M2,2) - 3;
		skewvec.at(i) =  (sqrt(GeneExp.n_cols) * M3) / pow(M2,1.5);
		//SDevvec.at(i) =  sqrt(M2/(GeneExp.n_cols - 1));
		
		//Linear regression
		//meanvec.at(i) = 1 + (mean - minMean)/range;
		//meanvec.at(i) = 1 + i;
		meanvec.at(i) = mean;
		//SumKurt += kurtvec.at(i);
		SumSkew += skewvec.at(i);
		//SumSD += SDevvec.at(i);
		SumMean += meanvec.at(i);
		//SumKurtMean += kurtvec.at(i) * meanvec.at(i);
		SumSkewMean += skewvec.at(i) * meanvec.at(i);
		//SumSDMean += SDevvec.at(i) * meanvec.at(i);
		SumMeanSq += pow(meanvec.at(i),2);
	}
	//y = Mx + C, here are the M and C results of the linear regression
	//double kurtM = (GeneExp.n_rows * SumKurtMean - SumKurt * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) ) ;
	//double kurtC = (SumKurt - kurtM * SumMean)/GeneExp.n_rows;
	double skewM = (GeneExp.n_rows * SumSkewMean - SumSkew * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) ) ;
	double skewC = (SumSkew - skewM * SumMean)/GeneExp.n_rows;
	//double SDevM = (GeneExp.n_rows * SumSDMean - SumSD * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) );
	
	//integral of the linear equation of skewness
	double Skewzerointercept = -skewC/skewM;
	double Skewintegral = 0;
	skewM = skewM/2;
	if (Skewzerointercept > meanvec.at(0) && Skewzerointercept < meanvec.at(skewvec.n_elem-1)) {
		Skewintegral = (abs(skewM * (meanvec.at(skewvec.n_elem-1) + Skewzerointercept) + skewC) * (meanvec.at(skewvec.n_elem-1) - Skewzerointercept) + abs(skewM * (meanvec.at(0) + Skewzerointercept) + skewC)* (Skewzerointercept - meanvec.at(0)))/(meanvec.at(skewvec.n_elem-1) - meanvec.at(0));
	} else {
		Skewintegral = abs(skewM * (meanvec.at(skewvec.n_elem-1) + meanvec.at(0)) + skewC);
	}
	/*
	//integral of the linear equation of kurtosis
	double Kurtzerointercept = -kurtC/kurtM;
	double Kurtintegral = 0;
	kurtM = kurtM/2;
	if (Kurtzerointercept > meanvec.at(0) && Kurtzerointercept < meanvec.at(kurtvec.n_elem-1)) {
		Kurtintegral = (abs(kurtM * (meanvec.at(kurtvec.n_elem-1) + Kurtzerointercept) + kurtC) * (meanvec.at(kurtvec.n_elem-1) - Kurtzerointercept) + abs(kurtM * (meanvec.at(0) + Kurtzerointercept) + kurtC)* (Kurtzerointercept - meanvec.at(0)) )/(meanvec.at(kurtvec.n_elem-1) - meanvec.at(0));
	} else {
		Kurtintegral = abs(kurtM * (meanvec.at(kurtvec.n_elem-1) + meanvec.at(0)) + kurtC);
	}
	*/
	//return 2 * abs(SDevM) + Skewintegral + Kurtintegral;
	return abs(Skewintegral);
}

double SkewVar (const arma::mat& GeneExp, const double& lambda2) {
	//vectors to store skewness, standard deviation and mean of each gene/feature.
	arma::vec skewvec(GeneExp.n_rows);
	arma::vec SDevvec(GeneExp.n_rows);
	arma::vec meanvec(GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	double SumSkew = 0;
	double SumSD = 0;
	double SumMean = 0;
	double SumSkewMean = 0;
	double SumSDMean = 0;
	double SumMeanSq = 0;
	
	//One pass linear regression with one pass variance, skewness
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		double M3 = 0;
		double delta, delta_n, term1;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			delta = log1p(GeneExp.at(i,n) * lambda2) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			numData++;
		}
		if (numData < 3) {
			continue;
		}
		//Here, calculate skewness and SD using 3rd and 2nd moments.
		skewvec.at(i) =  (sqrt(GeneExp.n_cols) * M3) / pow(M2,1.5);
		SDevvec.at(i) =  sqrt(M2/(GeneExp.n_cols - 1));
		
		//Linear regression
		meanvec.at(i) = mean;
		SumSkew += skewvec.at(i);
		SumSD += SDevvec.at(i);
		SumMean += meanvec.at(i);
		SumSkewMean += skewvec.at(i) * meanvec.at(i);
		SumSDMean += SDevvec.at(i) * meanvec.at(i);
		SumMeanSq += pow(meanvec.at(i),2);
	}
	//y = Mx + C, here are the M and C results of the linear regression
	double skewM = (GeneExp.n_rows * SumSkewMean - SumSkew * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) ) ;
	double skewC = (SumSkew - skewM * SumMean)/GeneExp.n_rows;
	double SDevM = (GeneExp.n_rows * SumSDMean - SumSD * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) );
	
	//integral of the linear equation of skewness
	double Skewzerointercept = -skewC/skewM;
	double Skewintegral = 0;
	skewM = skewM/2;
	if (Skewzerointercept > meanvec.at(0) && Skewzerointercept < meanvec.at(skewvec.n_elem-1)) {
		Skewintegral = (abs(skewM * (meanvec.at(skewvec.n_elem-1) + Skewzerointercept) + skewC) * (meanvec.at(skewvec.n_elem-1) - Skewzerointercept) + abs(skewM * (meanvec.at(0) + Skewzerointercept) + skewC)* (Skewzerointercept - meanvec.at(0)))/(meanvec.at(skewvec.n_elem-1) - meanvec.at(0));
	} else {
		Skewintegral = abs(skewM * (meanvec.at(skewvec.n_elem-1) + meanvec.at(0)) + skewC);
	}

	return pow(log1p(abs(SDevM))+1,2) + pow(log1p(Skewintegral)+1,2);
}
SEXP SkewVarCpp(SEXP xSEXP,SEXP ySEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	double lambda2 = Rcpp::as<double>(ySEXP);
	return Rcpp::wrap(SkewVar(GeneExp, lambda2));
}
arma::vec SkewAVar (const arma::mat& GeneExp, const double& lambda2) {
	//vectors to store skewness, standard deviation and mean of each gene/feature.
	arma::vec skewvec(GeneExp.n_rows);
	arma::vec SDevvec(GeneExp.n_rows);
	arma::vec meanvec(GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	double SumSkew = 0;
	double SumSD = 0;
	double SumMean = 0;
	double SumSkewMean = 0;
	double SumSDMean = 0;
	double SumMeanSq = 0;
	
	//One pass linear regression with one pass variance, skewness
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		double M3 = 0;
		double delta, delta_n, term1;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			delta = log1p(GeneExp.at(i,n) * lambda2) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			numData++;
		}
		if (numData < 3) {
			continue;
		}
		//Here, calculate skewness and SD using 3rd and 2nd moments.
		skewvec.at(i) =  (sqrt(GeneExp.n_cols) * M3) / pow(M2,1.5);
		SDevvec.at(i) =  sqrt(M2/(GeneExp.n_cols - 1));
		
		//Linear regression
		meanvec.at(i) = mean;
		SumSkew += skewvec.at(i);
		SumSD += SDevvec.at(i);
		SumMean += meanvec.at(i);
		SumSkewMean += skewvec.at(i) * meanvec.at(i);
		SumSDMean += SDevvec.at(i) * meanvec.at(i);
		SumMeanSq += pow(meanvec.at(i),2);
	}
	//y = Mx + C, here are the M and C results of the linear regression
	double skewM = (GeneExp.n_rows * SumSkewMean - SumSkew * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) ) ;
	double skewC = (SumSkew - skewM * SumMean)/GeneExp.n_rows;
	double SDevM = (GeneExp.n_rows * SumSDMean - SumSD * SumMean)/(GeneExp.n_rows * SumMeanSq - pow(SumMean,2) );
	
	//integral of the linear equation of skewness
	double Skewzerointercept = -skewC/skewM;
	double Skewintegral = 0;
	skewM = skewM/2;
	if (Skewzerointercept > meanvec.at(0) && Skewzerointercept < meanvec.at(skewvec.n_elem-1)) {
		Skewintegral = (abs(skewM * (meanvec.at(skewvec.n_elem-1) + Skewzerointercept) + skewC) * (meanvec.at(skewvec.n_elem-1) - Skewzerointercept) + abs(skewM * (meanvec.at(0) + Skewzerointercept) + skewC)* (Skewzerointercept - meanvec.at(0)))/(meanvec.at(skewvec.n_elem-1) - meanvec.at(0));
	} else {
		Skewintegral = abs(skewM * (meanvec.at(skewvec.n_elem-1) + meanvec.at(0)) + skewC);
	}
	arma::vec Answer(2);
	Answer.at(0) = SDevM;
	Answer.at(1) = Skewintegral;
	return Answer;
}
SEXP SkewAVarCpp(SEXP xSEXP,SEXP ySEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	double lambda2 = Rcpp::as<double>(ySEXP);
	return Rcpp::wrap(SkewAVar(GeneExp, lambda2));
}


//Given a range of lambda, this function finds lambda that minimizes F(lambda) (see article) based on the expression matrix by using binary search.
double LocalSearch(const arma::mat& GeneExp, double minBound, double maxBound, double& smallest, double search_exponent) {
	minBound = round(minBound);
	maxBound = round(maxBound);
	//cout << "here1 " << minBound<<   " " << maxBound << endl;
	double GminBound = minBound, GmaxBound = maxBound;
	double midBound = round((minBound + maxBound)/2);
	double om = SkewVar(GeneExp, midBound);
	double om2 = SkewVar(GeneExp,(midBound + 1));
	//Initialize smallestBound. This object is for the program to remember the smallest F(lambda) ever calculated. If the local minimal found is larger that this smallestBound, the boundaries will be reset and binary search will be rerun using a smaller boundary, with the smallestBound as center.
	double smallestBound;
	if (om > om2) {
		minBound = midBound + 1;
		smallest = om2;
		smallestBound = midBound + 1;
	} else {
		maxBound = midBound;
		smallest = om;
		smallestBound = midBound;
	}
	midBound = round((minBound + maxBound)/2);
	int index = 0;
	while (smallestBound != midBound) {
		//cout << "here1 " << smallestBound<<   " " << midBound << endl;
		if (index > 0) {
			//Reset boundary to center at smallestBound
			minBound = round(smallestBound - (smallestBound - GminBound)/pow(search_exponent,index));
			maxBound = round(smallestBound + (GmaxBound - smallestBound)/pow(search_exponent,index));
			midBound = round((minBound + maxBound)/2);
			//cout << minBound << " " << maxBound << endl;
		}
		index++;
		//Binary search for smallest F(lambda)
		while (minBound != midBound && maxBound != midBound) {
			om = SkewVar(GeneExp, midBound);
			om2 = SkewVar(GeneExp,midBound + 1);
			if (om > om2) {
				minBound = midBound;
				if (om2 < smallest) {
					smallest = om2;
					smallestBound = midBound + 1;
				}
			} else {
				maxBound = midBound;
				if (om < smallest) {
					smallest = om;
					smallestBound = midBound;
				}
			}
			midBound = round((minBound + maxBound)/2);
		}
	}
	return midBound;
}
double LocalSearchVarOnly(const arma::mat& GeneExp, double minBound, double maxBound, double& smallest, double search_exponent) {
	minBound = round(minBound);
	maxBound = round(maxBound);
	//cout << "here1 " << minBound<<   " " << maxBound << endl;
	double GminBound = minBound, GmaxBound = maxBound;
	double midBound = round((minBound + maxBound)/2);
	double om = VarOnly(GeneExp, midBound);
	double om2 = VarOnly(GeneExp,(midBound + 1));
	//Initialize smallestBound. This object is for the program to remember the smallest F(lambda) ever calculated. If the local minimal found is larger that this smallestBound, the boundaries will be reset and binary search will be rerun using a smaller boundary, with the smallestBound as center.
	double smallestBound;
	if (om > om2) {
		minBound = midBound + 1;
		smallest = om2;
		smallestBound = midBound + 1;
	} else {
		maxBound = midBound;
		smallest = om;
		smallestBound = midBound;
	}
	midBound = round((minBound + maxBound)/2);
	int index = 0;
	while (smallestBound != midBound) {
		//cout << "here1 " << smallestBound<<   " " << midBound << endl;
		if (index > 0) {
			//Reset boundary to center at smallestBound
			minBound = round(smallestBound - (smallestBound - GminBound)/pow(search_exponent,index));
			maxBound = round(smallestBound + (GmaxBound - smallestBound)/pow(search_exponent,index));
			midBound = round((minBound + maxBound)/2);
			//cout << minBound << " " << maxBound << endl;
		}
		index++;
		//Binary search for smallest F(lambda)
		while (minBound != midBound && maxBound != midBound) {
			om = VarOnly(GeneExp, midBound);
			om2 = VarOnly(GeneExp,midBound + 1);
			if (om > om2) {
				minBound = midBound;
				if (om2 < smallest) {
					smallest = om2;
					smallestBound = midBound + 1;
				}
			} else {
				maxBound = midBound;
				if (om < smallest) {
					smallest = om;
					smallestBound = midBound;
				}
			}
			midBound = round((minBound + maxBound)/2);
		}
	}
	return midBound;
}
double LocalSearchSkewOnly(const arma::mat& GeneExp, double minBound, double maxBound, double& smallest, double search_exponent) {
	minBound = round(minBound);
	maxBound = round(maxBound);
	//cout << "here1 " << minBound<<   " " << maxBound << endl;
	double GminBound = minBound, GmaxBound = maxBound;
	double midBound = round((minBound + maxBound)/2);
	double om = SkewOnly(GeneExp, midBound);
	double om2 = SkewOnly(GeneExp,(midBound + 1));
	//Initialize smallestBound. This object is for the program to remember the smallest F(lambda) ever calculated. If the local minimal found is larger that this smallestBound, the boundaries will be reset and binary search will be rerun using a smaller boundary, with the smallestBound as center.
	double smallestBound;
	if (om > om2) {
		minBound = midBound + 1;
		smallest = om2;
		smallestBound = midBound + 1;
	} else {
		maxBound = midBound;
		smallest = om;
		smallestBound = midBound;
	}
	midBound = round((minBound + maxBound)/2);
	int index = 0;
	while (smallestBound != midBound) {
		//cout << "here1 " << smallestBound<<   " " << midBound << endl;
		if (index > 0) {
			//Reset boundary to center at smallestBound
			minBound = round(smallestBound - (smallestBound - GminBound)/pow(search_exponent,index));
			maxBound = round(smallestBound + (GmaxBound - smallestBound)/pow(search_exponent,index));
			midBound = round((minBound + maxBound)/2);
			//cout << minBound << " " << maxBound << endl;
		}
		index++;
		//Binary search for smallest F(lambda)
		while (minBound != midBound && maxBound != midBound) {
			om = SkewOnly(GeneExp, midBound);
			om2 = SkewOnly(GeneExp,midBound + 1);
			if (om > om2) {
				minBound = midBound;
				if (om2 < smallest) {
					smallest = om2;
					smallestBound = midBound + 1;
				}
			} else {
				maxBound = midBound;
				if (om < smallest) {
					smallest = om;
					smallestBound = midBound;
				}
			}
			midBound = round((minBound + maxBound)/2);
		}
	}
	return midBound;
}

//Using the LocalSearch function, performs iterated local search to find minimal lambda.
SEXP LocateLambdaCpp(SEXP xSEXP,SEXP ySEXP, SEXP zSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	double search_exponent = Rcpp::as<double>(ySEXP);
	double maxBound = Rcpp::as<double>(zSEXP);
   //Boundary of ILS are minBound and maxBound
	//Here, we define the range of lambda based on the dataset.
	double OriginalmaxBound = maxBound;
	double minBound = 1;
	//Find local minima
	double localminIntegral;
	//Intital minimal
	double localmin = LocalSearch(GeneExp, minBound, maxBound,localminIntegral,search_exponent);
	if (minBound + 10 >= localmin || maxBound - 10 <= localmin) {
		localmin = round((minBound + maxBound)/4);
		localminIntegral = SkewVar(GeneExp, localmin);
		//cout << "starting at middle" << endl;
	}
	
	//First, iterated local search starting on the left hand side of the local minima.  
	//LHS search
	int searchIndex = 1;
	//1.perturbing the current local minimum;
	double newminBound = localmin - 10 * search_exponent;
	double lastnewminBound = localmin;
	double newmin, newminIntegral;
	double finallocalmin = localmin + 1;
	int numMaxBoundIncrease = 0;
	while ( (finallocalmin != localmin || maxBound - 10 <= localmin || minBound + 10 >= localmin) && numMaxBoundIncrease <= 5) {
		if (minBound + 10 >= localmin || maxBound - 10 <= localmin) {
			//enlarge maxBound if global minimum is not found yet.
			minBound = maxBound/2;
			maxBound = maxBound * 10;
			localminIntegral = SkewVar(GeneExp, localmin);
			numMaxBoundIncrease++;
			//cout << "maxBound enlarged to " << maxBound << endl;
		}
		finallocalmin = localmin;
		//Note: instead of randomly choosing a new range to search, we perturb from the local minimal by increasing the range exponentially. We first search through left hand side(LHS), then search through right hand side(RHS) of the local minimal. The stop condition is that no new local minimal is found after a full LHS and RHS search. This allow us to keep the running time around mlog(n). We use a history recording binary search algorithm to find local minimum within a range.
		while (newminBound >= minBound) {
			//cout << "LHS " << searchIndex <<  " " << newmin <<endl;
			//2.applying local search after starting from the modified solution.
			newmin = LocalSearch(GeneExp, newminBound, lastnewminBound, newminIntegral,search_exponent);
			//If new minimum is smaller than the local minimal, reset local minimal and searchIndex.
			if (newminIntegral < localminIntegral) {
				localmin = newmin;
				localminIntegral = newminIntegral;
				lastnewminBound = (newminBound + lastnewminBound)/2;
				if (lastnewminBound + 10 >= localmin) {
					searchIndex++;
				} else {
					searchIndex = 2;
				}
				
				//1.perturbing the current local minimum;
				if (newminBound != minBound) {
					newminBound = localmin -  10 * pow(search_exponent,searchIndex);
					if (newminBound < minBound) {
						newminBound = minBound;
					}
				} else {
					newminBound = minBound - 1;
				}
			} else {
				searchIndex++;
				lastnewminBound = (newminBound + lastnewminBound)/2;
				//1.perturbing the current local minimum;
				if (newminBound != minBound) {
					newminBound = localmin -  10 * pow(search_exponent,searchIndex);
					if (newminBound < minBound) {
						newminBound = minBound;
					}
				} else {
					newminBound = minBound - 1;
				}
			}
		}

		
		//Lastly, iterated local search on the right hand side of the local minima.  
		//RHS search
		searchIndex = 1;
		//1.perturbing the current local maximum;
		double newmaxBound = localmin + 10 * search_exponent;
		double lastnewmaxBound = localmin;
		double newmax, newmaxIntegral;
		while (newmaxBound <= maxBound) {
			//cout << "RHS " << searchIndex <<  " " << newmax <<endl;
			//2.applying local search after starting from the modified solution.
			newmax = LocalSearch(GeneExp, lastnewmaxBound, newmaxBound, newmaxIntegral,search_exponent);
			//If new maximum is smaller than the local maximal, reset local maximal and searchIndex.
			if (newmaxIntegral < localminIntegral) {
				localmin = newmax;
				localminIntegral = newmaxIntegral;
				lastnewmaxBound = (newmaxBound + lastnewmaxBound)/2;
				if (lastnewmaxBound - 10 <= localmin) {
					searchIndex++;
				} else {
					searchIndex = 2;
				}
				//1.perturbing the current local maximum;
				if (newmaxBound != maxBound) {
					newmaxBound = localmin +  10 * pow(search_exponent,searchIndex);
					if (newmaxBound > maxBound) {
						
						newmaxBound = maxBound;
					}
				} else {
					newmaxBound = maxBound + 1;
				}
			} else {
				searchIndex++;
				lastnewmaxBound = (newmaxBound + lastnewmaxBound)/2;
				//1.perturbing the current local maximum;
				if (newmaxBound != maxBound) {
					newmaxBound = localmin +  10 * pow(search_exponent,searchIndex);
					if (newmaxBound > maxBound) {
						newmaxBound = maxBound;
					}
				} else {
					newmaxBound = maxBound + 1;
				}
			}
		}
	}
	if (numMaxBoundIncrease == 6 && (maxBound - 10 <= localmin || minBound + 10 >= localmin)) {
		//Safety: This should never happen on count datasets, but we will never know.
		//Failed to find global minimum. If both skew and SD got minimum, use their mean. Otherwise, use only one of their minimum. If none have minimum, use the original maxBound.
		double tryVar = LocalSearchVarOnly(GeneExp, 1, OriginalmaxBound, localminIntegral, search_exponent);
		double trySkew = LocalSearchSkewOnly(GeneExp, 1, OriginalmaxBound, localminIntegral, search_exponent);
		if ( (OriginalmaxBound - 10 <= tryVar || minBound + 10 >= tryVar) && (OriginalmaxBound - 10 <= trySkew || minBound + 10 >= trySkew)) {
			localmin = OriginalmaxBound;
		} else if ( (OriginalmaxBound - 10 <= tryVar || minBound + 10 >= tryVar) || (OriginalmaxBound - 10 <= trySkew || minBound + 10 >= trySkew)) {
			if (OriginalmaxBound - 10 <= tryVar || minBound + 10 >= tryVar) {
				localmin = trySkew;
			} else {
				localmin = tryVar;
			}
		} else {
			localmin = (trySkew + tryVar)/2;
		}
	}
	return Rcpp::wrap(localmin);
}

//Get Slope from x and y vectors
SEXP getSlopeCpp(SEXP xSEXP, SEXP ySEXP) {
	arma::vec xvec = Rcpp::as<arma::vec>(xSEXP);
	arma::vec yvec = Rcpp::as<arma::vec>(ySEXP);
	
	double SumX = 0;
	double SumXSq = 0;
	double SumXY = 0;
	double SumY = 0;
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(xvec.n_elem); i ++) {
		SumX += xvec.at(i);
		SumXSq += pow(xvec.at(i),2);
		SumXY += xvec.at(i) * yvec.at(i);
		SumY += yvec.at(i);
	}
	//y = Mx + C, here are the M and C results of the linear regression
	double Slope = (xvec.n_elem * SumXY - SumY * SumX)/(xvec.n_elem * SumXSq - pow(SumX,2) ) ;
	//double constant = (SumY.at(n) - Slope * SumX)/GeneExp.n_rows;
	return Rcpp::wrap(Slope);
}
SEXP LinearRegressionCpp(SEXP xSEXP, SEXP ySEXP) {
	arma::vec xvec = Rcpp::as<arma::vec>(xSEXP);
	arma::vec yvec = Rcpp::as<arma::vec>(ySEXP);
	
	double SumX = 0;
	double SumXSq = 0;
	double SumXY = 0;
	double SumY = 0;
	
	//One pass linear regression
	for (int i = 0; i < int(xvec.n_elem); i++) {
		SumX += xvec.at(i);
		SumXSq += pow(xvec.at(i),2);
		SumXY += xvec.at(i) * yvec.at(i);
		SumY += yvec.at(i);
	}
	//y = Mx + C, here are the M and C results of the linear regression
	arma::vec Answer(2);
	Answer.at(0) = (xvec.n_elem * SumXY - SumY * SumX)/(xvec.n_elem * SumXSq - pow(SumX,2) ) ;
	Answer.at(1) = (SumY - Answer.at(0) * SumX)/xvec.n_elem;
	return Rcpp::wrap(Answer);
}
SEXP LinearRegressionZeroCpp(SEXP xSEXP, SEXP ySEXP) {
	arma::vec xvec = Rcpp::as<arma::vec>(xSEXP);
	arma::vec yvec = Rcpp::as<arma::vec>(ySEXP);
	
	double SumXSq = 0;
	double SumXY = 0;
	
	//One pass linear regression
	for (int i = 0; i < int(xvec.n_elem); i ++) {
		SumXY += xvec.at(i) *  yvec.at(i);
		SumXSq += pow(xvec.at(i),2);
	}
	//y = Mx + C, here are the M and C results of the linear regression
	arma::vec Answer(2);
	Answer.at(0) = 0;
	Answer.at(1) = SumXY/SumXSq;
	return Rcpp::wrap(Answer);
}
SEXP LinearRegressionFPCpp(SEXP xSEXP, SEXP ySEXP,SEXP x1SEXP, SEXP y1SEXP) {
	arma::vec xvec = Rcpp::as<arma::vec>(xSEXP);
	arma::vec yvec = Rcpp::as<arma::vec>(ySEXP);
	
	double x1 = Rcpp::as<double>(x1SEXP);
	double y1 = Rcpp::as<double>(y1SEXP);
	
	double SumX = 0;
	double SumXSq = 0;
	double SumXY = 0;
	double SumY = 0;

	//One pass linear regression
	for (int i = 0; i < int(xvec.n_elem); i ++) {
		SumX += xvec.at(i);
		SumXSq += pow(xvec.at(i),2);
		SumXY += xvec.at(i) * yvec.at(i);
		SumY += yvec.at(i);
	}
	//y = Mx + C, here are the M and C results of the linear regression
	arma::vec Answer(2);
	Answer.at(1) = (x1*y1*SumX-y1*SumXSq-pow(x1,2)*SumY+x1*SumXY)/(2*x1*SumX-SumXSq-pow(x1,2)*xvec.n_elem);
	Answer.at(0) = (y1 - Answer.at(1))/x1;
	return Rcpp::wrap(Answer);
}

//Algorithm to normalize batch effect
SEXP BatchEffectCpp(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP z2SEXP){
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	arma::mat Output = Rcpp::as<arma::mat>(ySEXP);
	arma::vec meanvec = Rcpp::as<arma::vec>(zSEXP);
	double BE_strength = Rcpp::as<double>(z2SEXP);
	
	//Objects to store the sums of variables that will be needed to perform linear regression.

	arma::vec SumX(GeneExp.n_cols);
	SumX.fill(0);
	arma::vec SumXSq(GeneExp.n_cols);
	SumXSq.fill(0);
	arma::vec SumXY(GeneExp.n_cols);
	SumXY.fill(0);
	arma::vec SumY(GeneExp.n_cols);
	SumY.fill(0);
	arma::vec nRows(GeneExp.n_cols);
	nRows.fill(0);
	
	//cout << SumMean << endl;
	double store;
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (std::isnan(GeneExp.at(i,n)) || GeneExp.at(i,n) == 0) {
				continue;
			}
			store = log(GeneExp.at(i,n));
			SumXY.at(n) += store *  meanvec.at(i);
			SumX.at(n) += store;
			nRows.at(n)++;
			SumXSq.at(n) += pow( store,2);
			SumY.at(n) += meanvec.at(i);
		}
	}
	//y = Mx + C, here are the M and C results of the linear regression
	double add = (1/BE_strength) - 1;
	for (int n = 0; n < int(Output.n_cols); n++) {
		double Slope = ((nRows.at(n) * SumXY.at(n) - SumY.at(n) * SumX.at(n))/(nRows.at(n) * SumXSq.at(n) - pow(SumX.at(n),2) )+ add) * BE_strength;
		double constant = ((SumY.at(n) - Slope * SumX.at(n))/nRows.at(n)) * BE_strength;
		//cout << Slope << " " << constant << endl;
		for (int i = 0; i < int(Output.n_rows); i ++) {
			if (std::isnan(Output.at(i,n)) || Output.at(i,n) == 0) {
				continue;
			}
			Output.at(i,n) = exp(Slope * log(Output.at(i,n)) + constant);
		}
	}
	
	return Rcpp::wrap(Output);
}

//Misc R functions for mean, sd and skewness calculations. Different ones are used in different situations. They are all implemented to improve running time.
SEXP NZrowLogMeanSDSkewCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::mat Answer(3,GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		double M3 = 0;
		double delta, delta_n, term1;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (GeneExp.at(i,n) == 0) {
				continue;
			}
			delta = log(GeneExp.at(i,n)) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			numData++;
		}
		if (numData == 0) {
			Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = 0;
		} else if (numData == 1) {
			Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = mean;
		} else if (numData == 2) {
			Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(1,i) = sqrt(M2/(GeneExp.n_cols - 1));
			Answer.at(0,i) = mean;
		} else {
			Answer.at(2,i) = (sqrt(GeneExp.n_cols) * M3) / pow(M2,1.5);
			Answer.at(1,i) = sqrt(M2/(GeneExp.n_cols - 1));
			Answer.at(0,i) = mean;
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP rowLogMeanSDSkewCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::mat Answer(3,GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		double M3 = 0;
		double delta, delta_n, term1;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (std::isnan(GeneExp.at(i,n)) || GeneExp.at(i,n) == 0) {
				continue;
			}
			delta = log(GeneExp.at(i,n)) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			numData++;
		}
		if (numData == 0) {
			Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = 0;
		} else if (numData == 1) {
			Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = mean;
		} else if (numData == 2) {
			Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(1,i) = sqrt(M2/(GeneExp.n_cols - 1));
			Answer.at(0,i) = mean;
		} else {
			Answer.at(2,i) = (sqrt(GeneExp.n_cols) * M3) / pow(M2,1.5);
			Answer.at(1,i) = sqrt(M2/(GeneExp.n_cols - 1));
			Answer.at(0,i) = mean;
		}
	}
	return Rcpp::wrap(Answer);
}

SEXP NZrowLog1pMeanSDSkewCpp (SEXP xSEXP, SEXP ySEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	double lambda = Rcpp::as<double>(ySEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::mat Answer(3,GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		double M3 = 0;
		double delta, delta_n, term1;
		int numData = 0;

		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (GeneExp.at(i,n) == 0) {
				continue;
			}
			delta = log1p(GeneExp.at(i,n) * lambda) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			numData++;
		}
		if (numData == 0) {
			Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = 0;
		} else if (numData == 1) {
			Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = mean;
		} else if (numData == 2) {
			Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(1,i) = sqrt(M2/(GeneExp.n_cols - 1));
			Answer.at(0,i) = mean;
		} else {
			Answer.at(2,i) = (sqrt(GeneExp.n_cols) * M3) / pow(M2,1.5);
			Answer.at(1,i) = sqrt(M2/(GeneExp.n_cols - 1));
			Answer.at(0,i) = mean;
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP NZrowLog1pMeanSDCpp (SEXP xSEXP, SEXP ySEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	double lambda = Rcpp::as<double>(ySEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::mat Answer(2,GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		//double M3 = 0;
		double delta, delta_n, term1;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (GeneExp.at(i,n) == 0) {
				continue;
			}
			delta = log1p(GeneExp.at(i,n) * lambda) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			//M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			numData++;
		}
		if (numData == 0) {
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = 0;
		} else if (numData == 1) {
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = mean;
		} else {
			Answer.at(1,i) = sqrt(M2/(GeneExp.n_cols - 1));
			Answer.at(0,i) = mean;
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP NZrowLog1pMeanSqrtSDCpp (SEXP xSEXP, SEXP ySEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	double lambda = Rcpp::as<double>(ySEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::mat Answer(2,GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		//double M3 = 0;
		double delta, delta_n, term1;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (GeneExp.at(i,n) == 0) {
				continue;
			}
			delta = log1p(GeneExp.at(i,n) * lambda) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			//M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			numData++;
		}
		if (numData == 0) {
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = 0;
		} else if (numData == 1) {
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = mean;
		} else {
			Answer.at(1,i) = pow(M2/(GeneExp.n_cols - 1),0.25);
			Answer.at(0,i) = mean;
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP rowLog1pMeanSDCpp (SEXP xSEXP, SEXP ySEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	double lambda = Rcpp::as<double>(ySEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::mat Answer(2,GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		//double M3 = 0;
		double delta, delta_n, term1;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (std::isnan(GeneExp.at(i,n))) {
				continue;
			}
			delta = log1p(GeneExp.at(i,n) * lambda) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			//M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			numData++;
		}
		if (numData == 0) {
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = 0;
		} else if (numData == 1) {
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = mean;
		} else {
			Answer.at(1,i) = sqrt(M2/(GeneExp.n_cols - 1));
			Answer.at(0,i) = mean;
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP HZrowLog1pMeanSDCpp (SEXP xSEXP, SEXP ySEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	double lambda = Rcpp::as<double>(ySEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::mat Answer(2,GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		//double M3 = 0;
		double delta, delta_n, term1;
		int numData = 0;
		bool includeZero=false;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (GeneExp.at(i,n) == 0) {
				if (includeZero) {
					includeZero = true;
					continue;
				} else {
					includeZero = true;
				}
			}
			delta = log1p(GeneExp.at(i,n) * lambda) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			//M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			numData++;
		}
		if (numData == 0) {
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = 0;
		} else if (numData == 1) {
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = mean;
		} else {
			Answer.at(1,i) = sqrt(M2/(GeneExp.n_cols - 1));
			Answer.at(0,i) = mean;
		}
	}
	return Rcpp::wrap(Answer);
}

SEXP NZrowMeanSDSkewCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::mat Answer(3,GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		double M3 = 0;
		double delta, delta_n, term1;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (GeneExp.at(i,n) == 0) {
				continue;
			}
			delta = GeneExp.at(i,n) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			numData++;
		}
		if (numData == 0) {
			Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = 0;
		} else if (numData == 1) {
			Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = mean;
		} else if (numData == 2) {
			Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(1,i) = sqrt(M2/(GeneExp.n_cols - 1));
			Answer.at(0,i) = mean;
		} else {
			Answer.at(2,i) = (sqrt(GeneExp.n_cols) * M3) / pow(M2,1.5);
			Answer.at(1,i) = sqrt(M2/(GeneExp.n_cols - 1));
			Answer.at(0,i) = mean;
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP NZrowMeanSDCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::mat Answer(2,GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		double delta, delta_n, term1;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (GeneExp.at(i,n) == 0) {
				continue;
			}
			delta = GeneExp.at(i,n) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M2 = M2 + term1;
			numData++;
		}
		if (numData == 0) {
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = 0;
		} else if (numData == 1) {
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = mean;
		} else {
			Answer.at(1,i) = sqrt(M2/(GeneExp.n_cols - 1));
			Answer.at(0,i) = mean;
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP rowMeanSDSkewCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::mat Answer(3,GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		double M3 = 0;
		double delta, delta_n, term1;
		int numData = 0;
		int numNonZero = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (std::isnan(GeneExp.at(i,n))) {
				continue;
			} else if (GeneExp.at(i,n) != 0) {
				numNonZero++;
			}
			delta = GeneExp.at(i,n) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M3 = M3 + term1 * delta_n * (numData - 1) - 3 * delta_n * M2;
			M2 = M2 + term1;
			numData++;
		}
		if (numNonZero <= 2) {
			Answer.at(2,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(1,i) = std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = std::numeric_limits<double>::quiet_NaN();
		} else {
			Answer.at(2,i) = (sqrt(GeneExp.n_cols) * M3) / pow(M2,1.5);
			Answer.at(1,i) = sqrt(M2/(GeneExp.n_cols - 1));
			Answer.at(0,i) = mean;
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP rowMeanSDCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::mat Answer(2,GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		double delta, delta_n, term1;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (std::isnan(GeneExp.at(i,n))) {
				continue;
			}
			delta = GeneExp.at(i,n) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M2 = M2 + term1;
			numData++;
		}
		if (numData <= 1) {
			Answer.at(1,i) =  std::numeric_limits<double>::quiet_NaN();
			Answer.at(0,i) = std::numeric_limits<double>::quiet_NaN();
		} else {
			Answer.at(1,i) =  sqrt(M2/(GeneExp.n_cols - 1));
			Answer.at(0,i) = mean;
		}
	}
	return Rcpp::wrap(Answer);
}


SEXP rowVarsCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::vec Answer(GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		double delta, delta_n, term1;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (std::isnan(GeneExp.at(i,n))) {
				continue;
			}
			delta = GeneExp.at(i,n) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M2 = M2 + term1;
			numData++;
		}
		if (numData <= 1) {
			Answer.at(i) = std::numeric_limits<double>::quiet_NaN();
		} else {
			Answer.at(i) =  M2/(GeneExp.n_cols - 1);
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP rowSDsCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::vec Answer(GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		double M2 = 0;
		double delta, delta_n, term1;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (std::isnan(GeneExp.at(i,n))) {
				continue;
			}
			delta = GeneExp.at(i,n) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M2 = M2 + term1;
			numData++;
		}
		if (numData <= 1) {
			Answer.at(i) = std::numeric_limits<double>::quiet_NaN();
		} else {
			Answer.at(i) =  sqrt(M2/(GeneExp.n_cols - 1));
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP colSDsCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::vec Answer(GeneExp.n_cols);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_cols); i ++) {
		double mean = 0;
		double M2 = 0;
		double delta, delta_n, term1;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_rows); n++) {
			if (std::isnan(GeneExp.at(n,i))) {
				continue;
			}
			delta = GeneExp.at(n,i) - mean;
			delta_n = delta / (numData + 1);
			term1 = delta * delta_n * numData;
			mean = mean + delta_n;
			M2 = M2 + term1;
			numData++;
		}
		if (numData <= 1) {
			Answer.at(i) = std::numeric_limits<double>::quiet_NaN();
		} else {
			Answer.at(i) =  sqrt(M2/(GeneExp.n_rows - 1));
		}
	}
	return Rcpp::wrap(Answer);
}
SEXP NZrowMeansCpp(SEXP xSEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::vec Answer(GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	//One pass linear regression with one pass variance, skewness and kurtosis
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		//double M2 = 0;
		double delta, delta_n;
		int numData = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (GeneExp.at(i,n) == 0) {
				continue;
			}
			delta = GeneExp.at(i,n) - mean;
			delta_n = delta / (numData + 1);
			mean = mean + delta_n;
			numData++;
		}
		Answer.at(i) = mean;
	}
	return Rcpp::wrap(Answer);
}
SEXP WNZrowMeansCpp(SEXP xSEXP, SEXP ySEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	arma::vec weights = Rcpp::as<arma::vec>(ySEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::vec Answer(GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	//At least two pass, because we need to check for zero and renew weights
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double sumNewWeights = 0;
		int numNZData = 0;
		arma::vec NZIndex(GeneExp.n_cols);
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			if (std::isnan(GeneExp.at(i,n)) || GeneExp.at(i,n) == 0) {
				continue;
			}
			sumNewWeights += weights.at(n);
			NZIndex.at(numNZData) = n;
			numNZData++;
		}
		if (numNZData == 0) {
			Answer.at(i) = 0;
			continue;
		}
		double mean = 0;
		for (int n = 0; n < numNZData; n++) {
			mean += GeneExp.at(i,NZIndex.at(n)) * weights.at(NZIndex.at(n)) / sumNewWeights;
		}
		Answer.at(i) = mean;
	}
	return Rcpp::wrap(Answer);
}
SEXP WrowMeansCpp(SEXP xSEXP, SEXP ySEXP) {
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	arma::vec weights = Rcpp::as<arma::vec>(ySEXP);
	//vectors to store skewness, standard deviation, kurtosis and mean of each gene/feature.
	arma::vec Answer(GeneExp.n_rows);
	//Objects to store the sums of variables that will be needed to perform linear regression.
	
	for (int i = 0; i < int(GeneExp.n_rows); i ++) {
		double mean = 0;
		for (int n = 0; n < int(GeneExp.n_cols); n++) {
			mean += weights.at(n) * GeneExp.at(i,n);
		}
		Answer.at(i) = mean;
	}
	return Rcpp::wrap(Answer);
}

static const R_CallMethodDef callMethods[] = {
	{"LocateLambdaCpp", (DL_FUNC) &LocateLambdaCpp, 3},
	{"SkewVarCpp", (DL_FUNC) &SkewVarCpp, 2},
	{"SkewAVarCpp", (DL_FUNC) &SkewAVarCpp, 2},
	{"NZrowLogMeanSDSkewCpp", (DL_FUNC) &NZrowLogMeanSDSkewCpp, 1},
	{"rowLogMeanSDSkewCpp", (DL_FUNC) &rowLogMeanSDSkewCpp, 1},
	{"NZrowLog1pMeanSDCpp", (DL_FUNC) &NZrowLog1pMeanSDCpp, 2},
	{"NZrowLog1pMeanSqrtSDCpp", (DL_FUNC) &NZrowLog1pMeanSqrtSDCpp, 2},
	{"rowLog1pMeanSDCpp", (DL_FUNC) &rowLog1pMeanSDCpp, 2},
	{"NZrowMeanSDSkewCpp", (DL_FUNC) &NZrowMeanSDSkewCpp, 1},
	{"NZrowMeanSDCpp", (DL_FUNC) &NZrowMeanSDCpp, 1},
	{"rowMeanSDSkewCpp", (DL_FUNC) &rowMeanSDSkewCpp, 1},
	{"rowMeanSDCpp", (DL_FUNC) &rowMeanSDCpp, 1},
	{"HZrowLog1pMeanSDCpp", (DL_FUNC) &HZrowLog1pMeanSDCpp, 2},
	{"rowVarsCpp", (DL_FUNC) &rowVarsCpp, 1},
	{"rowSDsCpp", (DL_FUNC) &rowSDsCpp, 1},
	{"colSDsCpp", (DL_FUNC) &colSDsCpp, 1},
	{"NZrowMeansCpp", (DL_FUNC) &NZrowMeansCpp, 1},
	{"WNZrowMeansCpp", (DL_FUNC) &WNZrowMeansCpp, 2},
	{"WrowMeansCpp", (DL_FUNC) &WrowMeansCpp, 2},
	
	{"XPMCpp", (DL_FUNC) &XPMCpp, 1},
	{"BatchEffectCpp", (DL_FUNC) &BatchEffectCpp, 4},
	{"getSlopeCpp", (DL_FUNC) &getSlopeCpp, 2},
	{"LinearRegressionCpp", (DL_FUNC) &LinearRegressionCpp, 2},
	{"LinearRegressionZeroCpp", (DL_FUNC) &LinearRegressionZeroCpp, 2},
	{"LinearRegressionFPCpp", (DL_FUNC) &LinearRegressionFPCpp, 4},
	{NULL, NULL, 0}
};

extern "C" {
	void R_init_Linnorm(DllInfo *info) { 
		R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	}
}
