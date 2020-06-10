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

// [[Rcpp::export]]
List SCANG_O_Search_Relatedness(arma::sp_mat G, arma::mat P, arma::vec residuals, const double threshold_o, const double threshold_s, const double threshold_b, const int Lmax, const int Lmin, int steplength, arma::mat weights_B, arma::mat weights_S, const int begid, const double filter, const double f)
{
	int ij, i, j, k, ii, jj, kk, s, ss, rr;

	// Number of weights
	int w_num = weights_B.n_cols;

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
	//p_s
	double sump_os = 1.0;
	//p_b
	double sump_ob = 1.0;

	// Weights for ACAT
	// SCANG_O
	arma::vec CCT_weights;
	CCT_weights.ones(2*w_num);
	// SCANG_S
	arma::vec CCT_weights_s;
	CCT_weights_s.ones(w_num);
	//SCANG_B
	arma::vec CCT_weights_b;
	CCT_weights_b.ones(w_num);

	// p-value
	// SCANG_O
	arma::vec CCT_p;
	CCT_p.ones(2*w_num);
	// SCANG_S
	arma::vec CCT_ps;
	CCT_ps.ones(w_num);
	// SCANG_B
	arma::vec CCT_pb;
	CCT_pb.ones(w_num);

	// Searching Algorithm
	// SCANG_O
	int num = 0;
	arma::mat candidate(1, 4);
	candidate.zeros();

	double summax = -100000.0;
	arma::mat candidatemax(1,4);
	candidatemax.zeros();

	// SCANG_S
	int num_s = 0;
	arma::mat candidate_s(1, 4);
	candidate_s.zeros();

	double summax_s = -100000.0;
	arma::mat candidatemax_s(1,4);
	candidatemax_s.zeros();

	// SCANG_B
	int num_b = 0;
	arma::mat candidate_b(1, 4);
	candidate_b.zeros();

	double summax_b = -100000.0;
	arma::mat candidatemax_b(1,4);
	candidatemax_b.zeros();

	int p = G.n_cols;

	// U-scores
	arma::rowvec x = trans(residuals)*G;

	// Weights Matrix
	arma::mat W;
	W.zeros(p,p);

	// Cov
	arma::mat Cov;
	Cov.zeros(p,p);

	// Cov matrix with weights
	// Burden
	arma::mat Covw_B;
	Covw_B.zeros(w_num*p,p);
	// SKAT
	arma::mat Covw_S;
	Covw_S.zeros(w_num*p,p);

	Cov = trans(trans(P*G)*G);

	// weights_vec
	// Burden
	arma::vec weights_vec_B;
	weights_vec_B.zeros(p);
	// SKAT
	arma::vec weights_vec_S;
	weights_vec_S.zeros(p);



	for(k = 0; k < w_num; k++)
	{
		// Burden
		for(kk = 0; kk < p; kk++)
		{
			weights_vec_B[kk] = weights_B(kk,k);
		}

		W.each_col() = weights_vec_B;
		Covw_B(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1)) = W%Cov;
		W.each_row() = trans(weights_vec_B);
		Covw_B(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1)) = Covw_B(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1))%W;

		// SKAT
		for(kk = 0; kk < p; kk++)
		{
			weights_vec_S[kk] = weights_S(kk,k);
		}

		W.each_col() = weights_vec_S;
		Covw_S(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1)) = W%Cov;
		W.each_row() = trans(weights_vec_S);
		Covw_S(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1)) = Covw_S(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1))%W;

	}


	int lengthnum = (Lmax-Lmin)/steplength + 1;

	double filtervalue = 0.0;
	double filtervalue_b = R::qchisq(filter,1,false,false);


	for (ij = 0; ij < lengthnum; ij++)
	{
		i = Lmin + ij*steplength;

		filtervalue = (R::qchisq(filter,i,false,false) - i)/sqrt(2*i);

		// initial values
		for(ss = 0; ss < w_num; ss++)
		{
			sum0_s(ss) = 0;
			sum0_b(ss) = 0;
			sumx_b(ss) = 0;
			sump_b(ss) = 1;
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
			w[ss] = 0;
		}

		// Q_burden and Q_skat
		for(ss = 0; ss < w_num; ss++)
		{
			for (k = 0; k < i; k++)
			{
				sum0_s(ss) = sum0_s(ss) + pow(x(k),2)*pow(weights_S(k,ss),2);
				sum0_b(ss) = sum0_b(ss) + x(k)*weights_B(k,ss);
			}
		}

		// Moments
		for(k = 0; k < i; k++)
		{
			for(ss = 0; ss < w_num; ss++)
			{
				w1(ss) = w1(ss) + Covw_S(ss*p + k, k);
			}
		}



		for(ss = 0; ss < w_num; ss++)
		{
			for (ii = 0; ii < (i - 1); ii++)
			{
				for (jj = (ii + 1); jj < i; jj++)
				{
					w2(ss) = w2(ss) + Covw_S(ss*p + ii, ii)*Covw_S(ss*p + jj, jj) - Covw_S(ss*p + ii, jj)*Covw_S(ss*p + jj, ii);
				}
			}
		}


		for(ss = 0; ss < w_num; ss++)
		{
			c1[ss] = w1[ss];
			c2[ss] = pow(w1[ss],2) - 2*w2[ss];

			if(i > 1)
			{
				w[ss] = arma::accu(Covw_B(arma::span(ss*p, ss*p + i - 1), arma::span(0, i - 1)));
			}
			if(i == 1)
			{
				w[ss] = Covw_B(ss*p + i - 1, i - 1);
			}

		}

		for(ss = 0; ss < w_num; ss++)
		{
			// Burden
			sumx_b(ss) = pow(sum0_b(ss),2)/w[ss];

			// SKAT
			tstar(ss) = (sum0_s(ss) - c1(ss))/sqrt(2*c2(ss));
			if(tstar(ss) > filtervalue)
			{
				if(testcal(ss) < 1)
				{
					eigenvals(arma::span(0,i - 1),ss) = arma::eig_sym(Covw_S(arma::span(ss*p + 0,ss*p + i - 1),arma::span(0,i - 1)));
					testcal(ss) = 1;
					for(rr = 0; rr < i; rr++)
					{
						if(eigenvals(rr,ss) < 1e-8)
						{
							eigenvals(rr,ss) = 0.0;
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
					sump_s(ss) = Liumod(Covw_S(arma::span(ss*p + 0,ss*p + i - 1),arma::span(0,i-1)),tstar(ss));
				}
			}
		}

		Li = 0;

		for(ss = 0; ss < w_num; ss++)
		{
			if((sumx_b(ss) > filtervalue_b)||(sump_s(ss) < filter))
			{
				Li = 1;
			}
		}

		if(Li > 0)
		{
			for(ss = 0; ss < w_num; ss++)
			{
				sump_b(ss) = R::pchisq(sumx_b(ss),1,false,false);

				if(tstar(ss) <= filtervalue)
				{
					if(testcal(ss) < 1)
					{
						eigenvals(arma::span(0,i - 1),ss) = arma::eig_sym(Covw_S(arma::span(ss*p + 0,ss*p + i - 1),arma::span(0,i - 1)));
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
						sump_s(ss) = Liumod(Covw_S(arma::span(ss*p + 0,ss*p + i - 1),arma::span(0,i-1)),tstar(ss));
					}
				}
			}

			for(ss = 0; ss < w_num; ss++)
			{
				// SCANG_O
				CCT_p[ss*2] = sump_s(ss);
				CCT_p[ss*2 + 1] = sump_b(ss);
				// SCANG_S
				CCT_ps[ss] = sump_s(ss);
				// SCANG_B
				CCT_pb[ss] = sump_b(ss);
			}

			// SCANG_O
			sump_o = CCT_pval(CCT_p,CCT_weights);
			sump_o = -log(sump_o);
			// SCANG_S
			sump_os = CCT_pval(CCT_ps,CCT_weights_s);
			sump_os = -log(sump_os);
			// SCANG_B
			sump_ob = CCT_pval(CCT_pb,CCT_weights_b);
			sump_ob = -log(sump_ob);

			// Signal_regions
			// SCANG_O
			if(sump_o > threshold_o)
			{
				num = num + 1;
				candidate.resize(num, 4);
				candidate(num - 1, 0) = sump_o;
				candidate(num - 1, 1) = 1;
				candidate(num - 1, 2) = i;
			}
			// SCANG_S
			if(sump_os > threshold_s)
			{
				num_s = num_s + 1;
				candidate_s.resize(num, 4);
				candidate_s(num - 1, 0) = sump_os;
				candidate_s(num - 1, 1) = 1;
				candidate_s(num - 1, 2) = i;
			}
			// SCANG_B
			if(sump_ob > threshold_b)
			{
				num_b = num_b + 1;
				candidate_b.resize(num, 4);
				candidate_b(num - 1, 0) = sump_os;
				candidate_b(num - 1, 1) = 1;
				candidate_b(num - 1, 2) = i;
			}

			// TOP 1 Region
			// SCANG_O
			if(sump_o > summax)
			{
				summax = sump_o;
				candidatemax(0,0) = sump_o;
				candidatemax(0,1) = 1 + begid - 1;
				candidatemax(0,2) = i + begid - 1;
			}

			// SCANG_S
			if(sump_os > summax_s)
			{
				summax_s = sump_os;
				candidatemax_s(0,0) = sump_os;
				candidatemax_s(0,1) = 1 + begid - 1;
				candidatemax_s(0,2) = i + begid - 1;
			}

			// SCANG_B
			if(sump_ob > summax_b)
			{
				summax_b = sump_ob;
				candidatemax_b(0,0) = sump_ob;
				candidatemax_b(0,1) = 1 + begid - 1;
				candidatemax_b(0,2) = i + begid - 1;
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
				w1[ss] = w1[ss] - Covw_S(ss*p + j - 1, j - 1) + Covw_S(ss*p + j + i - 1, j + i - 1);

				for (kk = 1; kk < i; kk++)
				{
					w2[ss] = w2[ss] - (Covw_S(ss*p + j - 1, j - 1)*Covw_S(ss*p + j - 1 + kk, j - 1 + kk) - Covw_S(ss*p + j - 1, j - 1 + kk)*Covw_S(ss*p + j - 1 + kk, j - 1));
					w2[ss] = w2[ss] + (Covw_S(ss*p + j + i - 1, j + i - 1)*Covw_S(ss*p + j - 1 + kk, j - 1 + kk) - Covw_S(ss*p + j + i - 1, j - 1 + kk)*Covw_S(ss*p + j - 1 + kk, j + i - 1));
				}

				c1[ss] = w1[ss];
				c2[ss] = pow(w1[ss],2) - 2*w2[ss];

				for(rr = 0; rr < Lmax; rr++)
				{
					eigenvals(rr,ss) = 0;
				}

				if (i > 1)
				{
					w[ss] = w[ss] - arma::accu(Covw_B(arma::span(ss*p + j - 1, ss*p + j - 1), arma::span(j, j + i - 2))) - arma::accu(Covw_B(arma::span(ss*p +j, ss*p + j + i - 2), arma::span(j - 1, j - 1))) - Covw_B(ss*p +j - 1, j - 1);
					w[ss] = w[ss] + arma::accu(Covw_B(arma::span(ss*p + j + i - 1, ss*p + j + i - 1), arma::span(j, j + i - 2))) + arma::accu(Covw_B(arma::span(ss*p + j, ss*p + j + i - 2), arma::span(j + i - 1, j + i - 1))) + Covw_B(ss*p +j + i - 1, j + i - 1);
				}
				if (i == 1)
				{
					w[ss] = w[ss] - Covw_B(ss*p + j - 1, j - 1) + Covw_B(ss*p + j + i - 1, j + i - 1);
				}
			}

			for(ss = 0; ss < w_num; ss++)
			{
				sump_b(ss) = 1;
				sump_s(ss) = 1;

				// Burden
				sum0_b(ss) = sum0_b(ss) - x(j - 1)*weights_B(j - 1,ss) + x(j + i - 1)*weights_B(j + i - 1,ss);
				sumx_b(ss) = pow(sum0_b(ss),2)/w[ss];


				// SKAT
				sum0_s(ss) = sum0_s(ss) - pow(x(j - 1), 2)*pow(weights_S(j - 1,ss), 2) + pow(x(j + i - 1), 2)*pow(weights_S(j + i - 1,ss), 2);
				tstar(ss) = (sum0_s(ss) - c1[ss])/sqrt(2*c2[ss]);

				if(tstar(ss) > filtervalue)
				{
					if(testcal[ss] < 1)
					{
						eigenvals(arma::span(0,i - 1),ss) = arma::eig_sym(Covw_S(arma::span(ss*p + j,ss*p + j + i - 1),arma::span(j, j + i - 1)));
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
						sump_s(ss) = Liumod(Covw_S(arma::span(ss*p + j,ss*p + j + i - 1),arma::span(j,j + i - 1)),tstar(ss));
					}
				}

			}


			Li = 0;

			for(ss = 0; ss < w_num; ss++)
			{
				if((sumx_b(ss) > filtervalue_b)||(sump_s(ss) < filter))
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

					// SKAT

					if(tstar(ss) <= filtervalue)
					{
						if(testcal[ss] < 1)
						{
							eigenvals(arma::span(0,i - 1),ss) = arma::eig_sym(Covw_S(arma::span(ss*p + j, ss*p + j + i - 1),arma::span(j, j + i - 1)));
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
							sump_s(ss) = Liumod(Covw_S(arma::span(ss*p + j,ss*p + j + i - 1),arma::span(j,j + i - 1)),tstar(ss));
						}
					}

				}

				for(ss = 0; ss < w_num; ss++)
				{
					// SCANG_O
					CCT_p[ss*2] = sump_s(ss);
					CCT_p[ss*2 + 1] = sump_b(ss);
					// SCANG_S
					CCT_ps[ss] = sump_s(ss);
					// SCANG_B
					CCT_pb[ss] = sump_b(ss);
				}

				// SCANG_O
				sump_o = CCT_pval(CCT_p,CCT_weights);
				sump_o = -log(sump_o);
				// SCANG_S
				sump_os = CCT_pval(CCT_ps,CCT_weights_s);
				sump_os = -log(sump_os);
				// SCANG_B
				sump_ob = CCT_pval(CCT_pb,CCT_weights_b);
				sump_ob = -log(sump_ob);

				// Signal_regions
				// SCANG_O
				if(sump_o > threshold_o)
				{
					num = num + 1;
					candidate.resize(num, 4);
					candidate(num - 1, 0) = sump_o;
					candidate(num - 1, 1) = j + 1;
					candidate(num - 1, 2) = j + i;
				}
				// SCANG_S
				if(sump_os > threshold_s)
				{
					num_s = num_s + 1;
					candidate_s.resize(num_s, 4);
					candidate_s(num_s - 1, 0) = sump_os;
					candidate_s(num_s - 1, 1) = j + 1;
					candidate_s(num_s - 1, 2) = j + i;
				}
				// SCANG_B
				if(sump_ob > threshold_b)
				{
					num_b = num_b + 1;
					candidate_b.resize(num_b, 4);
					candidate_b(num_b - 1, 0) = sump_ob;
					candidate_b(num_b - 1, 1) = j + 1;
					candidate_b(num_b - 1, 2) = j + i;
				}

				// TOP 1 Region
				// SCANG_O
				if(sump_o > summax)
				{
					summax = sump_o;
					candidatemax(0,0) = sump_o;
					candidatemax(0,1) = j + 1 + begid - 1;
					candidatemax(0,2) = j + i + begid - 1;
				}

				// SCANG_S
				if(sump_os > summax_s)
				{
					summax_s = sump_os;
					candidatemax_s(0,0) = sump_os;
					candidatemax_s(0,1) = j + 1 + begid - 1;
					candidatemax_s(0,2) = j + i + begid - 1;
				}

				// SCANG_B
				if(sump_ob > summax_b)
				{
					summax_b = sump_ob;
					candidatemax_b(0,0) = sump_ob;
					candidatemax_b(0,1) = j + 1 + begid - 1;
					candidatemax_b(0,2) = j + i + begid - 1;
				}
			}

		}

	}

	// SCANG_O
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

	// SCANG_S
	arma::uvec indices_s = sort_index(-candidate_s.col(0));
	candidate_s = candidate_s.rows(indices_s);

	double loc_left_s = 0;
	double loc_right_s = 0;

	for (ii = 0; ii < (num_s-1); ii++)
	{
		if (candidate_s(ii,3) < 1)
		{
			for (jj = ii + 1; jj < num_s; jj++)
			{
				if(candidate_s(ii, 1) < candidate_s(jj, 1))
				{
					loc_left_s = candidate_s(jj, 1);
				}else
				{
					loc_left_s = candidate_s(ii, 1);
				}
				if(candidate_s(ii, 2) < candidate_s(jj, 2))
				{
					loc_right_s = candidate_s(ii, 2);
				}else
				{
					loc_right_s = candidate_s(jj, 2);
				}

				if (loc_right_s > loc_left_s - 1)
				{
					if((loc_right_s-loc_left_s + 1)/(candidate_s(jj,2) - candidate_s(jj,1) + 1) > f)
					{
						candidate_s(jj, 3) = 1;
					}


				}
			}
		}
	}

	int num1_s = 0;
	arma::mat res_s(1, 4);
	res_s.zeros();
	for (kk = 0; kk<num_s; kk++)
	{
		if (candidate_s(kk, 3)<1)
		{
			num1_s = num1_s + 1;
			res_s.resize(num1_s, 4);
			res_s(num1_s - 1, 0) = candidate_s(kk, 0);
			res_s(num1_s - 1, 1) = candidate_s(kk, 1) + begid - 1;
			res_s(num1_s - 1, 2) = candidate_s(kk, 2) + begid - 1;
		}
	}

	// SCANG_B
	arma::uvec indices_b = sort_index(-candidate_b.col(0));
	candidate_b = candidate_b.rows(indices_b);

	double loc_left_b = 0;
	double loc_right_b = 0;

	for (ii = 0; ii < (num_b-1); ii++)
	{
		if (candidate_b(ii,3) < 1)
		{
			for (jj = ii + 1; jj < num_b; jj++)
			{
				if(candidate_b(ii, 1) < candidate_b(jj, 1))
				{
					loc_left_b = candidate_b(jj, 1);
				}else
				{
					loc_left_b = candidate_b(ii, 1);
				}
				if(candidate_b(ii, 2) < candidate_b(jj, 2))
				{
					loc_right_b = candidate_b(ii, 2);
				}else
				{
					loc_right_b = candidate_b(jj, 2);
				}

				if (loc_right_b > loc_left_b - 1)
				{
					if((loc_right_b-loc_left_b + 1)/(candidate_b(jj,2) - candidate_b(jj,1) + 1) > f)
					{
						candidate_b(jj, 3) = 1;
					}


				}
			}
		}
	}

	int num1_b = 0;
	arma::mat res_b(1, 4);
	res_b.zeros();
	for (kk = 0; kk < num_b; kk++)
	{
		if (candidate_b(kk, 3)<1)
		{
			num1_b = num1_b + 1;
			res_b.resize(num1_b, 4);
			res_b(num1_b - 1, 0) = candidate_b(kk, 0);
			res_b(num1_b - 1, 1) = candidate_b(kk, 1) + begid - 1;
			res_b(num1_b - 1, 2) = candidate_b(kk, 2) + begid - 1;
		}
	}

	return List::create(Named("res_o") = res, Named("resmost_o") = candidatemax, Named("res_s") = res_s, Named("resmost_s") = candidatemax_s, Named("res_b") = res_b, Named("resmost_b") = candidatemax_b);


}





