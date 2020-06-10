// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// declare CCT_pval
double CCT_pval(arma::vec x, arma::vec weights);


// [[Rcpp::export]]
List SCANG_B_Search_Relatedness_sp(arma::sp_mat G, arma::sp_mat Sigma_i, arma::mat Sigma_iX, arma::mat cov, arma::vec residuals, const double threshold, const int Lmax, const int Lmin, int steplength, arma::mat weights, const int begid, const double filter, const double f)
{
	int ij, i, j, k, ii, jj, kk, ss;

	// Number of weights
	int w_num = weights.n_cols;

	int Li = 0;

	// Intermediate parameter for Burden
	arma::vec w;
	w.zeros(w_num);

	// Sum of U-Scores
	arma::vec sum0_b;
	sum0_b.zeros(w_num);

	// Q_burden
	arma::vec sumx_b;
	sumx_b.zeros(w_num);

	// p_burden
	arma::vec sump_b;
	sump_b.ones(w_num);

	// Intermediate parameter for Omnibus test
	//p_O
	double sump_o = 1.0;

	// Weights for ACAT
	arma::vec CCT_weights;
	CCT_weights.ones(w_num);

	// p-value
	arma::vec CCT_p;
	CCT_p.ones(w_num);


	// Searching Algorithm
	int num = 0;
	arma::mat candidate(1, 4);
	candidate.zeros();

	double summax = -100000.0;
	arma::mat candidatemax(1,4);
	candidatemax.zeros();

	int p = G.n_cols;
	int q = Sigma_iX.n_cols;

	// U-scores
	arma::rowvec x = trans(residuals)*G;

	// t(X)*G
	arma::mat tSigma_iX_G;
	tSigma_iX_G.zeros(q,p);

	// Weights Matrix
	arma::mat W;
	W.zeros(p,p);

	// Cov
	arma::mat Cov;
	Cov.zeros(p,p);

	// Cov matrix with weights
	arma::mat Covw;
	Covw.zeros(w_num*p,p);

	tSigma_iX_G = trans(Sigma_iX)*G;
	Cov = trans(trans(Sigma_i*G)*G) - trans(tSigma_iX_G)*cov*tSigma_iX_G;

	// weights_vec
	arma::vec weights_vec;
	weights_vec.zeros(p);

	for(k = 0; k < w_num; k++)
	{
		for(kk = 0; kk < p; kk++)
		{
			weights_vec[kk] = weights(kk,k);
		}

		W.each_col() = weights_vec;
		Covw(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1)) = W%Cov;
		W.each_row() = trans(weights_vec);
		Covw(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1)) = Covw(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1))%W;
	}


	int lengthnum = (Lmax-Lmin)/steplength + 1;

	double filtervalue_b = R::qchisq(filter,1,false,false);


	for (ij = 0; ij < lengthnum; ij++)
	{
		i = Lmin + ij*steplength;

		// initial values
		for(ss = 0; ss < w_num; ss++)
		{
			sum0_b(ss) = 0;
			sumx_b(ss) = 0;
			sump_b(ss) = 1;
		}

		for(ss = 0; ss < w_num; ss++)
		{
			w[ss] = 0;
		}

		// Q_burden
		for(ss = 0; ss < w_num; ss++)
		{
			for (k = 0; k < i; k++)
			{
				sum0_b(ss) = sum0_b(ss) + x(k)*weights(k,ss);
			}
		}

		for(ss = 0; ss < w_num; ss++)
		{
		  if(i > 1)
		  {
		    w[ss] = arma::accu(Covw(arma::span(ss*p, ss*p + i - 1), arma::span(0, i - 1)));
		  }
		  if(i == 1)
		  {
		    w[ss] = Covw(ss*p + i - 1, i - 1);
		  }

		}

		for(ss = 0; ss < w_num; ss++)
		{
			// Burden
			sumx_b(ss) = pow(sum0_b(ss),2)/w[ss];
		}

		Li = 0;

		for(ss = 0; ss < w_num; ss++)
		{
			if(sumx_b(ss) > filtervalue_b)
			{
				Li = 1;
			}
		}

		if(Li > 0)
		{
			for(ss = 0; ss < w_num; ss++)
			{
				sump_b(ss) = R::pchisq(sumx_b(ss),1,false,false);
			}

			for(ss = 0; ss < w_num; ss++)
			{
				CCT_p[ss] = sump_b(ss);
			}

			sump_o = CCT_pval(CCT_p,CCT_weights);
			sump_o = -log(sump_o);

			if(sump_o > threshold)
			{
				num = num + 1;
				candidate.resize(num, 4);
				candidate(num - 1, 0) = sump_o;
				candidate(num - 1, 1) = 1;
				candidate(num - 1, 2) = i;
			}

			if(sump_o > summax)
			{
				summax = sump_o;
				candidatemax(0,0) = sump_o;
				candidatemax(0,1) = 1 + begid - 1;
				candidatemax(0,2) = i + begid - 1;
			}
		}


		for (j = 1; j < (p - i + 1); j++)
		{


			for(ss = 0; ss < w_num; ss++)
			{
				if (i > 1)
				{
					w[ss] = w[ss] - arma::accu(Covw(arma::span(ss*p + j - 1, ss*p + j - 1), arma::span(j, j + i - 2))) - arma::accu(Covw(arma::span(ss*p +j, ss*p + j + i - 2), arma::span(j - 1, j - 1))) - Covw(ss*p +j - 1, j - 1);
					w[ss] = w[ss] + arma::accu(Covw(arma::span(ss*p + j + i - 1, ss*p + j + i - 1), arma::span(j, j + i - 2))) + arma::accu(Covw(arma::span(ss*p + j, ss*p + j + i - 2), arma::span(j + i - 1, j + i - 1))) + Covw(ss*p +j + i - 1, j + i - 1);
				}
				if (i == 1)
				{
					w[ss] = w[ss] - Covw(ss*p + j - 1, j - 1) + Covw(ss*p + j + i - 1, j + i - 1);
				}
			}

			for(ss = 0; ss < w_num; ss++)
			{
				sump_b(ss) = 1;

				// Burden
				sum0_b(ss) = sum0_b(ss) - x(j - 1)*weights(j - 1,ss) + x(j + i - 1)*weights(j + i - 1,ss);
				sumx_b(ss) = pow(sum0_b(ss),2)/w[ss];
			}


			Li = 0;

			for(ss = 0; ss < w_num; ss++)
			{
				if(sumx_b(ss) > filtervalue_b)
				{
					Li = 1;
				}
			}

			if(Li > 0)
			{
				for(ss = 0; ss < w_num; ss++)
				{
					// Burden
					sump_b(ss) = R::pchisq(sumx_b(ss),1,false,false);
				}

				for(ss = 0; ss < w_num; ss++)
				{
					CCT_p[ss] = sump_b(ss);
				}

				sump_o = CCT_pval(CCT_p,CCT_weights);
				sump_o = -log(sump_o);

				if(sump_o > threshold)
				{
					num = num + 1;
					candidate.resize(num, 4);
					candidate(num - 1, 0) = sump_o;
					candidate(num - 1, 1) = j + 1;
					candidate(num - 1, 2) = j + i;
				}

				if(sump_o > summax)
				{
					summax = sump_o;
					candidatemax(0,0) = sump_o;
					candidatemax(0,1) = j + 1 + begid - 1;
					candidatemax(0,2) = j + i + begid - 1;
				}
			}

		}

	}


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
			res(num1 - 1, 1) = candidate(kk, 1) + begid - 1;
			res(num1 - 1, 2) = candidate(kk, 2) + begid - 1;
		}
	}
	return List::create(Named("res") = res, Named("resmost") = candidatemax);

}




