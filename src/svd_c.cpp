// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h> 
using namespace Rcpp;

// [[Rcpp::export]]
List svd_c(arma::mat X) {

	arma::mat U;
	arma::vec s;
	arma::mat V;

	arma::svd(U,s,V,X);
	
	return List::create(Named("U") = U, Named("s") = s, Named("V") = V);
	
}







