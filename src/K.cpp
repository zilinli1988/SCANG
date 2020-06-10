// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


// [[Rcpp::export]]
double K(double x, arma::vec egvalues)
{
    double res = 0.0;
    const int en = egvalues.size();

    for(int i = 0; i < en; i++)
    {
        res = res + log(1-2*egvalues(i)*x);
    }

    res = res*(-0.5);

    return res;
}



