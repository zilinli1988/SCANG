// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


// [[Rcpp::export]]
double Liumod(arma::mat Cov, double q)
{
	double c2 = 0.0;
	double c4 = 0.0;
	int n = Cov.n_cols;
	int i;
	double l = 0.0;
	double res = 0.0;

	bool lower = false;
	bool logp = false;

	Cov = Cov*Cov;
	for(i = 0; i < n; i++)
	{
		c2 = c2 + Cov(i,i);
	}
	Cov = Cov*Cov;
	for(i = 0; i < n; i++)
	{
		c4 = c4 + Cov(i,i);
	}
	l = pow(c2,2)/c4;
	res = R::pchisq(q*sqrt(2*l)+l,l,lower,logp);

	return res;

}




