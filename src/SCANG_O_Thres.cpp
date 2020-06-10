// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_64BIT_WORD 1
#include <STAAR.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
using namespace STAAR;
using namespace Rcpp;

// declare Liumod
double Liumod(arma::mat Cov, double q);
// declare maxO
arma::mat maxO(int p, int Lmax, int Lmin, int steplength, arma::mat x, arma::mat weights_B, arma::mat weights_S, arma::mat Cov, double filter, int times);


// [[Rcpp::export]]
arma::mat SCANG_O_Thres(arma::sp_mat G, arma::mat X, arma::vec working, double sigma, int fam, int times, int Lmax, int Lmin, int steplength, arma::mat weights_B, arma::mat weights_S, double filter) {

	// variants number
	int p = G.n_cols;
	// sample size
	int n = G.n_rows;
	// covariates number
	int q = X.n_cols;

	// pesudo-residual
	arma::mat y;
	y = arma::randn<arma::mat>(times,n);

	// pesudo-score
	arma::mat x;
	x.zeros(times,p);

	// t(X)*G
	arma::mat tX_G;
	tX_G.zeros(q,p);

	// Cov
	arma::mat Cov;
	Cov.zeros(p,p);

	if(fam == 0)
	{
		tX_G = trans(X)*G;
		Cov = trans(G)*G - trans(tX_G)*inv(trans(X)*X)*tX_G;
		x = (y*G - y*X*inv(trans(X)*X)*tX_G)*sigma;
	}else
	{
		tX_G = trans(X)*(arma::diagmat(working))*G;
		Cov = trans(G)*arma::diagmat(working)*G - trans(tX_G)*inv(trans(X)*arma::diagmat(working)*X)*tX_G;
		y = y*arma::diagmat(sqrt(working));
		x = (y*G - y*X*inv(trans(X)*arma::diagmat(working)*X)*tX_G)*sigma;

	}

	Cov = Cov*pow(sigma,2);

	return maxO(p, Lmax, Lmin, steplength, trans(x), weights_B, weights_S, Cov, filter, times);

}


