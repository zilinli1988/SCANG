// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat regionfilter(arma::mat candidate, const double f)
{
	int ii,jj,kk;
	int num =  candidate.n_rows;

	arma::uvec indices = sort_index(-candidate.col(0));
	candidate = candidate.rows(indices);

	double loc_left = 0;
	double loc_right = 0;

	for (ii = 0; ii < (num-1); ii++)
	{
		if (candidate(ii,3) < 1)
		{
			for (jj = ii + 1; jj < num; jj++)
			{
				if(candidate(ii, 1) < candidate(jj, 1))
				{
					loc_left = candidate(jj, 1);
				}else
				{
					loc_left = candidate(ii, 1);
				}
				if(candidate(ii, 2) < candidate(jj, 2))
				{
					loc_right = candidate(ii, 2);
				}else
				{
					loc_right = candidate(jj, 2);
				}

				if (loc_right > loc_left - 1)
				{
					if((loc_right-loc_left + 1)/(candidate(jj,2) - candidate(jj,1) + 1) > f)
					{
						candidate(jj, 3) = 1;
					}


				}
			}
		}
	}

	int num1 = 0;
	arma::mat res(1, 4);
	res.zeros();
	for (kk = 0; kk<num; kk++)
	{
		if (candidate(kk, 3)<1)
		{
			num1 = num1 + 1;
			res.resize(num1, 4);
			res(num1 - 1, 0) = candidate(kk, 0);
			res(num1 - 1, 1) = candidate(kk, 1);
			res(num1 - 1, 2) = candidate(kk, 2);
		}
	}
	return res;

}







