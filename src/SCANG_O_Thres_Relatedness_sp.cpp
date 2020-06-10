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
arma::mat SCANG_O_Thres_Relatedness_sp(arma::sp_mat G, arma::sp_mat Sigma_i, arma::mat Sigma_iX, arma::mat cov, arma::mat pesudo_residuals, int times, int Lmax, int Lmin, int steplength, arma::mat weights_B, arma::mat weights_S, double filter) {

	// variants number
	int p = G.n_cols;
	// covariates number
	int q = Sigma_iX.n_cols;

	// pesudo-score
	arma::mat x;
	x.zeros(times,p);

	// t(X)*G
	arma::mat tSigma_iX_G;
	tSigma_iX_G.zeros(q,p);

	// Cov
	arma::mat Cov;
	Cov.zeros(p,p);

	tSigma_iX_G = trans(Sigma_iX)*G;
	Cov = trans(trans(Sigma_i*G)*G) - trans(tSigma_iX_G)*cov*tSigma_iX_G;

	x = pesudo_residuals*G;

	return maxO(p, Lmax, Lmin, steplength, trans(x), weights_B, weights_S, Cov, filter, times);
}





