// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;
// declare K
double K(double x, arma::vec egvalues);
// declare K1 (first derivative)
double K1(double x, arma::vec egvalues, double q);
// declare K2 (second derivative)
double K2(double x, arma::vec egvalues);
// declare bisection
double Bisection(arma::vec egvalues, double q, double xmin, double xmax);
// declare saddlepoint
double Saddle(double q, arma::vec egvalues);
// declare maxO
arma::mat maxO(int p, int Lmax, int Lmin, int steplength, arma::mat x, arma::mat weights_B, arma::mat weights_S, arma::mat Cov, double filter, int times);
// declare Liumod
double Liumod(arma::mat Cov, double q);
// declare CCT_pval
double CCT_pval(arma::vec x, arma::vec weights);

// [[Rcpp::export]]
arma::mat SCANG_O_Thres_Relatedness(arma::sp_mat G, arma::mat P, arma::mat pesudo_residuals, int times, int Lmax, int Lmin, int steplength, arma::mat weights_B, arma::mat weights_S, double filter) {

	// variants number
	int p = G.n_cols;

	// pesudo-score
	arma::mat x;
	x.zeros(times,p);


	// Cov
	arma::mat Cov;
	Cov.zeros(p,p);

	Cov = trans(trans(P*G)*G);

	x = pesudo_residuals*G;

	return maxO(p, Lmax, Lmin, steplength, trans(x), weights_B, weights_S, Cov, filter, times);

}



