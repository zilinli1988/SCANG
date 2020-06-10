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
// declare Liumod
double Liumod(arma::mat Cov, double q);
// declare CCT_pval
double CCT_pval(arma::vec x, arma::vec weights);


// [[Rcpp::export]]
List SCANG_S_Search_Relatedness_sp(arma::sp_mat G, arma::sp_mat Sigma_i, arma::mat Sigma_iX, arma::mat cov, arma::vec residuals, const double threshold, const int Lmax, const int Lmin, int steplength, arma::mat weights, const int begid, const double filter, const double f)
{
	int ij, i, j, k, ii, jj, kk, s, ss, rr;

	// Number of weights
	int w_num = weights.n_cols;

	int Li = 0;

	// Intermediate parameter for SKAT
	arma::vec w1;
	w1.zeros(w_num);

	arma::vec w2;
	w2.zeros(w_num);

	// Moments
	arma::vec c1;
	c1.zeros(w_num);

	arma::vec c2;
	c2.zeros(w_num);

	// eigen dicomposition indicator
	arma::vec testcal;
	testcal.zeros(w_num);

	// Q_skat
	arma::vec sum0_s;
	sum0_s.zeros(w_num);

	// p_skat
	arma::vec sump_s;
	sump_s.ones(w_num);

	// tstar used for filtering
	arma::vec tstar;
	tstar.zeros(w_num);

	// eigen values
	arma::mat eigenvals;
	eigenvals.zeros(Lmax,w_num);

	arma::vec eigenvals_vec;
	eigenvals_vec.zeros(Lmax);


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

	double filtervalue = 0.0;

	for (ij = 0; ij < lengthnum; ij++)
	{
		i = Lmin + ij*steplength;

		filtervalue = (R::qchisq(filter,i,false,false) - i)/sqrt(2*i);

		// initial values
		for(ss = 0; ss < w_num; ss++)
		{
			sum0_s(ss) = 0;
			sump_s(ss) = 1;
		}

		for(s = 0; s < Lmax; s++)
		{
			for(ss = 0; ss < w_num; ss++)
			{
				eigenvals(s,ss) = 0;
			}
		}

		for(ss = 0; ss < w_num; ss++)
		{
			testcal[ss] = 0;
			w1[ss] = 0;
			w2[ss] = 0;
			c1[ss] = 0;
			c2[ss] = 0;
		}

		// Q_skat
		for(ss = 0; ss < w_num; ss++)
		{
			for (k = 0; k < i; k++)
			{
				sum0_s(ss) = sum0_s(ss) + pow(x(k),2)*pow(weights(k,ss),2);
			}
		}

		// Moments
		for(k = 0; k < i; k++)
		{
			for(ss = 0; ss < w_num; ss++)
			{
				w1(ss) = w1(ss) + Covw(ss*p + k, k);
			}
		}



		for(ss = 0; ss < w_num; ss++)
		{
			for (ii = 0; ii < (i - 1); ii++)
			{
				for (jj = (ii + 1); jj < i; jj++)
				{
					w2(ss) = w2(ss) + Covw(ss*p + ii, ii)*Covw(ss*p + jj, jj) - Covw(ss*p + ii, jj)*Covw(ss*p + jj, ii);
				}
			}
		}


		for(ss = 0; ss < w_num; ss++)
		{
			c1[ss] = w1[ss];
			c2[ss] = pow(w1[ss],2) - 2*w2[ss];

		}

		for(ss = 0; ss < w_num; ss++)
		{
			// SKAT
			tstar(ss) = (sum0_s(ss) - c1(ss))/sqrt(2*c2(ss));
			if(tstar(ss) > filtervalue)
			{
				if(testcal(ss) < 1)
				{
					eigenvals(arma::span(0,i - 1),ss) = arma::eig_sym(Covw(arma::span(ss*p + 0,ss*p + i - 1),arma::span(0,i - 1)));
					testcal(ss) = 1;
					for(rr = 0; rr < i; rr++)
					{
						if(eigenvals(rr,ss) < 1e-8)
						{
							eigenvals(rr,ss) = 0;
						}
					}

				}

				for(rr = 0; rr < i; rr++)
				{
					eigenvals_vec[rr] = eigenvals(rr,ss);
				}

				sump_s(ss) = Saddle(sum0_s(ss),eigenvals_vec(arma::span(0,i - 1)));
				if(sump_s(ss) == 2)
				{
					sump_s(ss) = Liumod(Covw(arma::span(ss*p + 0,ss*p + i - 1),arma::span(0,i-1)),tstar(ss));
				}
			}
		}

		Li = 0;

		for(ss = 0; ss < w_num; ss++)
		{
			if(sump_s(ss) < filter)
			{
				Li = 1;
			}
		}

		if(Li > 0)
		{
			for(ss = 0; ss < w_num; ss++)
			{
				if(tstar(ss) <= filtervalue)
				{
					if(testcal(ss) < 1)
					{
						eigenvals(arma::span(0,i - 1),ss) = arma::eig_sym(Covw(arma::span(ss*p + 0,ss*p + i - 1),arma::span(0,i - 1)));
						testcal[ss] = 1;
						for(rr = 0; rr < i; rr++)
						{
							if(eigenvals(rr,ss) < 1e-8)
							{
								eigenvals(rr,ss) = 0;
							}
						}

					}

					for(rr = 0; rr < i; rr++)
					{
						eigenvals_vec[rr] = eigenvals(rr,ss);
					}

					sump_s(ss) = Saddle(sum0_s(ss),eigenvals_vec(arma::span(0,i - 1)));
					if(sump_s(ss) == 2)
					{
						sump_s(ss) = Liumod(Covw(arma::span(ss*p + 0,ss*p + i - 1),arma::span(0,i-1)),tstar(ss));
					}
				}
			}

			for(ss = 0; ss < w_num; ss++)
			{
				CCT_p[ss] = sump_s(ss);
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
				candidatemax(0,1) = 1+begid-1;
				candidatemax(0,2) = i+begid-1;
			}
		}


		for (j = 1; j < (p - i + 1); j++)
		{
			for(ss = 0; ss < w_num; ss++)
			{
				testcal[ss] = 0;
			}

			for(ss = 0; ss < w_num; ss++)
			{
				w1[ss] = w1[ss] - Covw(ss*p + j - 1, j - 1) + Covw(ss*p + j + i - 1, j + i - 1);

				for (kk = 1; kk < i; kk++)
				{
					w2[ss] = w2[ss] - (Covw(ss*p + j - 1, j - 1)*Covw(ss*p + j - 1 + kk, j - 1 + kk) - Covw(ss*p + j - 1, j - 1 + kk)*Covw(ss*p + j - 1 + kk, j - 1));
					w2[ss] = w2[ss] + (Covw(ss*p + j + i - 1, j + i - 1)*Covw(ss*p + j - 1 + kk, j - 1 + kk) - Covw(ss*p + j + i - 1, j - 1 + kk)*Covw(ss*p + j - 1 + kk, j + i - 1));
				}

				c1[ss] = w1[ss];
				c2[ss] = pow(w1[ss],2) - 2*w2[ss];

				for(rr = 0; rr < Lmax; rr++)
				{
					eigenvals(rr,ss) = 0;
				}

			}

			for(ss = 0; ss < w_num; ss++)
			{
				sump_s(ss) = 1;

				// SKAT
				sum0_s(ss) = sum0_s(ss) - pow(x(j - 1), 2)*pow(weights(j - 1,ss), 2) + pow(x(j + i - 1), 2)*pow(weights(j + i - 1,ss), 2);
				tstar(ss) = (sum0_s(ss) - c1[ss])/sqrt(2*c2[ss]);

				if(tstar(ss) > filtervalue)
				{
					if(testcal[ss] < 1)
					{
						eigenvals(arma::span(0,i - 1),ss) = arma::eig_sym(Covw(arma::span(ss*p + j,ss*p + j + i - 1),arma::span(j, j + i - 1)));
						testcal[ss] = 1;
						for(rr = 0; rr < i; rr++)
						{
							if(abs(eigenvals(rr,ss)) < 1e-10)
							{
								eigenvals(rr,ss) = 0;
							}
						}

					}

					for(rr = 0; rr < i; rr++)
					{
						eigenvals_vec[rr] = eigenvals(rr,ss);
					}

					sump_s(ss) = Saddle(sum0_s(ss),eigenvals_vec(arma::span(0,i - 1)));
					if(sump_s(ss) == 2)
					{
						sump_s(ss) = Liumod(Covw(arma::span(ss*p + j,ss*p + j + i - 1),arma::span(j,j + i - 1)),tstar(ss));
					}
				}

			}


			Li = 0;

			for(ss = 0; ss < w_num; ss++)
			{
				if(sump_s(ss) < filter)
				{
					Li = 1;
				}
			}

			if(Li > 0)
			{
				for(ss = 0; ss < w_num; ss++)
				{
					// SKAT

					if(tstar(ss) <= filtervalue)
					{
						if(testcal[ss] < 1)
						{
							eigenvals(arma::span(0,i - 1),ss) = arma::eig_sym(Covw(arma::span(ss*p + j, ss*p + j + i - 1),arma::span(j, j + i - 1)));
							testcal[ss] = 1;
							for(rr = 0; rr < i; rr++)
							{
								if(eigenvals(rr,ss) < 1e-8)
								{
									eigenvals(rr,ss) = 0;
								}
							}

						}

						for(rr = 0; rr < i; rr++)
						{
							eigenvals_vec[rr] = eigenvals(rr,ss);
						}

						sump_s(ss) = Saddle(sum0_s(ss),eigenvals_vec(arma::span(0,i - 1)));
						if(sump_s(ss) == 2)
						{
							sump_s(ss) = Liumod(Covw(arma::span(ss*p + j,ss*p + j + i - 1),arma::span(j,j + i - 1)),tstar(ss));
						}
					}

				}

				for(ss = 0; ss < w_num; ss++)
				{
					CCT_p[ss] = sump_s(ss);
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
	for (kk = 0; kk < num; kk++)
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




