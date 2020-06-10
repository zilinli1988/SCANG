// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_64BIT_WORD 1
#include <STAAR.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
using namespace STAAR;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat maxO(int p, int Lmax, int Lmin, int steplength, arma::mat x, arma::mat weights_B, arma::mat weights_S, arma::mat Cov, double filter, int times)
{
	int ij, i, j, k, kk, r, rr, s, ss;

	int w_num = weights_B.n_cols;

	int Li = 0;

	// Cov matrix with weights
	// Burden
	arma::mat Covw_B;
	Covw_B.zeros(w_num*p,p);

	// arma::mat Covw2_B;
	// Covw2_B.zeros(w_num*p,p);

	// SKAT
	arma::mat Covw_S;
	Covw_S.zeros(w_num*p,p);

	arma::mat Covw2_S;
	Covw2_S.zeros(w_num*p,p);

	// Weights Matrix
	arma::mat W;
	W.zeros(p,p);

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

		// Covw2_B(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1)) = Covw_B(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1))%Covw_B(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1));

		// SKAT
		for(kk = 0; kk < p; kk++)
		{
			weights_vec_S[kk] = weights_S(kk,k);
		}

		W.each_col() = weights_vec_S;
		Covw_S(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1)) = W%Cov;
		W.each_row() = trans(weights_vec_S);
		Covw_S(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1)) = Covw_S(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1))%W;

		Covw2_S(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1)) = Covw_S(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1))%Covw_S(arma::span(k*p, (k+1)*p - 1), arma::span(0, p - 1));
	}


	// Intermediate parameter for SKAT
	arma::vec c1;
	c1.zeros(w_num);

	arma::vec c2;
	c2.zeros(w_num);

	arma::vec c4;
	c4.zeros(w_num);

	arma::vec l;
	l.zeros(w_num);

	arma::vec testcov;
	testcov.zeros(w_num);

	arma::vec testcal;
	testcal.zeros(w_num);

	arma::mat Covsq;
	Covsq.zeros(Lmax,Lmax);

	arma::mat sum0_s;
	sum0_s.zeros(w_num,times);

	arma::mat sump_s;
	sump_s.ones(w_num,times);

	arma::mat tstar;
	tstar.zeros(w_num,times);

	arma::mat eigenvals;
	eigenvals.zeros(Lmax,w_num);

	arma::vec eigenvals_vec;
	eigenvals_vec.zeros(Lmax);

	// Intermediate parameter for Burden
	arma::vec w;
	w.zeros(w_num);

	arma::mat sum0_b;
	sum0_b.zeros(w_num,times);

	arma::mat sumx_b;
	sumx_b.zeros(w_num,times);

	arma::mat sump_b;
	sump_b.ones(w_num,times);

	// Intermediate parameter for Omnibus test (0:Omnibus; 1:SKAT; 2:Burden)
	arma::mat sump_o;
	sump_o.ones(3,times);

	arma::vec CCT_weights;
	CCT_weights.ones(2*w_num);

	arma::vec CCT_p;
	CCT_p.ones(2*w_num);

	// thres
	// row 0: SCANG-O
	// row 1: SCANG-S
	// row 2: SCANG-B

	arma::mat threshold;
	threshold.zeros(3,times);

	double filtervalue = 0.0;
	double filtervalue_b = R::qchisq(filter,1,false,false);

	bool lower1 = true;
	bool logp1 = false;

	int lengthnum = (Lmax-Lmin)/steplength + 1;

	for (ij = 0; ij < lengthnum; ij++)
	{
		i = Lmin + ij*steplength;

		filtervalue = (R::qchisq(1-filter,i,lower1,logp1) - i)/sqrt(2*i);

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
			testcov[ss] = 0;
			c1[ss] = 0;
			c2[ss] = 0;
			c4[ss] = 0;
			w[ss] = 0;
			l[ss] = 0;
		}

		for(r = 0; r < times; r++)
		{
			for(ss = 0; ss < w_num; ss++)
			{
				sum0_s(ss,r) = 0;
				sum0_b(ss,r) = 0;
			}
		}


		for(r = 0; r < times; r++)
		{
			for(ss = 0; ss < w_num; ss++)
			{
				for (k = 0; k < i; k++)
				{
					sum0_s(ss,r) = sum0_s(ss,r) + pow(x(k,r),2)*pow(weights_S(k,ss),2);
					sum0_b(ss,r) = sum0_b(ss,r) + x(k,r)*weights_B(k,ss);
				}
			}
		}

		for(k = 0; k < i; k++)
		{
			for(ss = 0; ss < w_num; ss++)
			{
				c1[ss] = c1[ss] + Covw_S(ss*p + k, k);
			}
		}


		for(ss = 0; ss < w_num; ss++)
		{
			if(i > 1)
			{
				w[ss] = arma::accu(Covw_B(arma::span(ss*p, ss*p + i - 1), arma::span(0, i - 1)));
				c2[ss] = arma::accu(Covw2_S(arma::span(ss*p, ss*p + i - 1), arma::span(0, i - 1)));
			}
			if(i == 1)
			{
				w[ss] = Covw_B(ss*p + i - 1, i - 1);
				c2[ss] = Covw2_S(ss*p + i - 1, i - 1);
			}

		}


		for(r = 0; r < times; r++)
		{
			for(ss = 0; ss < w_num; ss++)
			{
				// Burden
				sumx_b(ss,r) = pow(sum0_b(ss,r),2)/w[ss];

				// SKAT
				tstar(ss,r) = (sum0_s(ss,r) - c1[ss])/sqrt(2*c2[ss]);
				if(tstar(ss,r) > filtervalue)
				{
					if(testcov[ss] < 1)
					{
						Covsq(arma::span(0,i-1),arma::span(0,i-1)) = Covw_S(arma::span(ss*p, ss*p + i - 1), arma::span(0, i - 1))*Covw_S(arma::span(ss*p, ss*p + i - 1), arma::span(0, i - 1));
						Covsq(arma::span(0,i-1),arma::span(0,i-1)) = Covsq(arma::span(0,i-1),arma::span(0,i-1))%Covsq(arma::span(0,i-1),arma::span(0,i-1));
						c4[ss] = arma::accu(Covsq(arma::span(0,i-1),arma::span(0,i-1)));
						testcov[ss] = 1;
					}

					l[ss] = pow(c2[ss],2)/c4[ss];

					sump_s(ss,r) = R::pchisq(tstar(ss,r)*sqrt(2*l[ss])+l[ss],l[ss],false,false);

					if(sump_s(ss,r) < filter)
					{
						if(testcal[ss] < 1)
						{
							eigenvals(arma::span(0,i - 1),ss) = arma::eig_sym(Covw_S(arma::span(ss*p + 0,ss*p + i - 1),arma::span(0,i - 1)));
							testcal[ss] = 1;
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


						sump_s(ss,r) = Saddle(sum0_s(ss,r),eigenvals_vec(arma::span(0,i - 1)));
						if(sump_s(ss,r) == 2)
						{
							sump_s(ss,r) = R::pchisq(tstar(ss,r)*sqrt(2*l[ss])+l[ss],l[ss],false,false);
						}
					}

				}

			}

			Li = 0;

			for(ss = 0; ss < w_num; ss++)
			{
				if((sumx_b(ss,r) > filtervalue_b)||(sump_s(ss,r) < filter))
				{
					Li = 1;
				}
			}


			if(Li > 0)
			{
				for(ss = 0; ss < w_num; ss++)
				{
					sump_b(ss,r) = R::pchisq(sumx_b(ss,r),1,false,false);

					if(tstar(ss,r) <= filtervalue)
					{
						if(testcov[ss] < 1)
						{
							Covsq(arma::span(0,i-1),arma::span(0,i-1)) = Covw_S(arma::span(ss*p, ss*p + i - 1), arma::span(0, i - 1)) * Covw_S(arma::span(ss*p, ss*p + i - 1), arma::span(0, i - 1));
							Covsq(arma::span(0,i-1),arma::span(0,i-1)) = Covsq(arma::span(0,i-1),arma::span(0,i-1))%Covsq(arma::span(0,i-1),arma::span(0,i-1));
							c4[ss] = arma::accu(Covsq(arma::span(0,i-1),arma::span(0,i-1)));
							testcov[ss] = 1;
						}

						l[ss] = pow(c2[ss],2)/c4[ss];

						sump_s(ss,r) = R::pchisq(tstar(ss,r)*sqrt(2*l[ss])+l[ss],l[ss],false,false);

						if(sump_s(ss,r) < filter)
						{
							if(testcal[ss] < 1)
							{
								eigenvals(arma::span(0,i - 1),ss) = arma::eig_sym(Covw_S(arma::span(ss*p + 0,ss*p + i - 1),arma::span(0,i - 1)));
								testcal[ss] = 1;
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

							sump_s(ss,r) = Saddle(sum0_s(ss,r),eigenvals_vec(arma::span(0,i - 1)));
							if(sump_s(ss,r) == 2)
							{
								sump_s(ss,r) = R::pchisq(tstar(ss,r)*sqrt(2*l[ss])+l[ss],l[ss],false,false);
							}
						}

					}
				}

				for(ss = 0; ss < w_num; ss++)
				{
					CCT_p[ss] = sump_s(ss,r);
					CCT_p[w_num + ss] = sump_b(ss,r);
				}

				// SCANG-O
				sump_o(0,r) = CCT_pval(CCT_p,CCT_weights);
				sump_o(0,r) = -log(sump_o(0,r));

				if(sump_o(0,r) > threshold(0,r))
				{
					threshold(0,r) = sump_o(0,r);
				}

				// SCANG-S
				sump_o(1,r) = CCT_pval(CCT_p(arma::span(0,w_num - 1)),CCT_weights(arma::span(0,w_num - 1)));
				sump_o(1,r) = -log(sump_o(1,r));

				if(sump_o(1,r) > threshold(1,r))
				{
					threshold(1,r) = sump_o(1,r);
				}

				// SCANG-B
				sump_o(2,r) = CCT_pval(CCT_p(arma::span(w_num,2*w_num - 1)),CCT_weights(arma::span(w_num,2*w_num - 1)));
				sump_o(2,r) = -log(sump_o(2,r));

				if(sump_o(2,r) > threshold(2,r))
				{
					threshold(2,r) = sump_o(2,r);
				}
			}
		}


		for (j = 1; j < (p - i + 1); j++)
		{
			for(ss = 0; ss < w_num; ss++)
			{
				testcal[ss] = 0;
				testcov[ss] = 0;
			}

			for(ss = 0; ss < w_num; ss++)
			{
				c1[ss] = c1[ss] - Covw_S(ss*p + j - 1, j - 1) + Covw_S(ss*p + j + i - 1, j + i - 1);

				for(rr = 0; rr < Lmax; rr++)
				{
					eigenvals(rr,ss) = 0.0;
				}

				if (i > 1)
				{
					w[ss] = w[ss] - arma::accu(Covw_B(arma::span(ss*p + j - 1, ss*p + j - 1), arma::span(j, j + i - 2))) - arma::accu(Covw_B(arma::span(ss*p +j, ss*p + j + i - 2), arma::span(j - 1, j - 1))) - Covw_B(ss*p +j - 1, j - 1);
					w[ss] = w[ss] + arma::accu(Covw_B(arma::span(ss*p + j + i - 1, ss*p + j + i - 1), arma::span(j, j + i - 2))) + arma::accu(Covw_B(arma::span(ss*p + j, ss*p + j + i - 2), arma::span(j + i - 1, j + i - 1))) + Covw_B(ss*p +j + i - 1, j + i - 1);

					c2[ss] = c2[ss] - arma::accu(Covw2_S(arma::span(ss*p + j - 1, ss*p + j - 1), arma::span(j, j + i - 2))) - arma::accu(Covw2_S(arma::span(ss*p +j, ss*p + j + i - 2), arma::span(j - 1, j - 1))) - Covw2_S(ss*p +j - 1, j - 1);
					c2[ss] = c2[ss] + arma::accu(Covw2_S(arma::span(ss*p + j + i - 1, ss*p + j + i - 1), arma::span(j, j + i - 2))) + arma::accu(Covw2_S(arma::span(ss*p + j, ss*p + j + i - 2), arma::span(j + i - 1, j + i - 1))) + Covw2_S(ss*p +j + i - 1, j + i - 1);
				}
				if (i == 1)
				{
					w[ss] = w[ss] - Covw_B(ss*p + j - 1, j - 1) + Covw_B(ss*p + j + i - 1, j + i - 1);
					c2[ss] = c2[ss] - Covw2_S(ss*p + j - 1, j - 1) + Covw2_S(ss*p + j + i - 1, j + i - 1);
				}
			}


			for(r = 0; r < times; r++)
			{
				for(ss = 0; ss < w_num; ss++)
				{
					sump_b(ss,r) = 1;
					sump_s(ss,r) = 1;

					// Burden
					sum0_b(ss,r) = sum0_b(ss,r) - x(j - 1, r)*weights_B(j - 1,ss) + x(j + i - 1, r)*weights_B(j + i - 1,ss);
					sumx_b(ss,r) = pow(sum0_b(ss,r),2)/w[ss];


				    // SKAT
					sum0_s(ss,r) = sum0_s(ss,r) - pow(x(j - 1, r), 2)*pow(weights_S(j - 1,ss), 2) + pow(x(j + i - 1, r), 2)*pow(weights_S(j + i - 1,ss), 2);
					tstar(ss,r) = (sum0_s(ss,r) - c1[ss])/sqrt(2*c2[ss]);

					if(tstar(ss,r) > filtervalue)
					{
						if(testcov[ss] < 1)
						{
							Covsq(arma::span(0,i-1),arma::span(0,i-1)) = Covw_S(arma::span(ss*p + j, ss*p + j + i - 1), arma::span(j, j + i - 1)) * Covw_S(arma::span(ss*p + j, ss*p + j + i - 1), arma::span(j, j + i - 1));
							Covsq(arma::span(0,i-1),arma::span(0,i-1)) = Covsq(arma::span(0,i-1),arma::span(0,i-1))%Covsq(arma::span(0,i-1),arma::span(0,i-1));
							c4[ss] = arma::accu(Covsq(arma::span(0,i-1),arma::span(0,i-1)));
							testcov[ss] = 1;
						}

						l[ss] = pow(c2[ss],2)/c4[ss];

						sump_s(ss,r) = R::pchisq(tstar(ss,r)*sqrt(2*l[ss])+l[ss],l[ss],false,false);

						if(sump_s(ss,r) < filter)
						{
							if(testcal[ss] < 1)
							{
								eigenvals(arma::span(0,i - 1),ss) = arma::eig_sym(Covw_S(arma::span(ss*p + j,ss*p + j + i - 1),arma::span(j, j + i - 1)));
								testcal[ss] = 1;
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

							sump_s(ss,r) = Saddle(sum0_s(ss,r),eigenvals_vec(arma::span(0,i - 1)));
							if(sump_s(ss,r) == 2)
							{
								sump_s(ss,r) = R::pchisq(tstar(ss,r)*sqrt(2*l[ss])+l[ss],l[ss],false,false);
							}
						}
					}

				}

				Li = 0;

				for(ss = 0; ss < w_num; ss++)
				{
					if((sumx_b(ss,r) > filtervalue_b)||(sump_s(ss,r) < filter))
					{
						Li = 1;
					}
				}


				if(Li > 0)
				{
					for(ss = 0; ss < w_num; ss++)
					{
						// Burden
						sump_b(ss,r) = R::pchisq(sumx_b(ss,r),1,false,false);

						// SKAT

						if(tstar(ss,r) <= filtervalue)
						{
							if(testcov[ss] < 1)
							{
								Covsq(arma::span(0,i-1),arma::span(0,i-1)) = Covw_S(arma::span(ss*p + j, ss*p + j + i - 1), arma::span(j, j + i - 1)) * Covw_S(arma::span(ss*p + j, ss*p + j + i - 1), arma::span(j, j + i - 1));
								Covsq(arma::span(0,i-1),arma::span(0,i-1)) = Covsq(arma::span(0,i-1),arma::span(0,i-1))%Covsq(arma::span(0,i-1),arma::span(0,i-1));
								c4[ss] = arma::accu(Covsq(arma::span(0,i-1),arma::span(0,i-1)));
								testcov[ss] = 1;
							}

							l[ss] = pow(c2[ss],2)/c4[ss];

							sump_s(ss,r) = R::pchisq(tstar(ss,r)*sqrt(2*l[ss])+l[ss],l[ss],false,false);

							if(sump_s(ss,r) < filter)
							{
								if(testcal[ss] < 1)
								{
									eigenvals(arma::span(0,i - 1),ss) = arma::eig_sym(Covw_S(arma::span(ss*p + j, ss*p + j + i - 1),arma::span(j, j + i - 1)));
									testcal[ss] = 1;
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

								sump_s(ss,r) = Saddle(sum0_s(ss,r),eigenvals_vec(arma::span(0,i - 1)));
								if(sump_s(ss,r) == 2)
								{
									sump_s(ss,r) = R::pchisq(tstar(ss,r)*sqrt(2*l[ss])+l[ss],l[ss],false,false);
								}
							}
						}

					}

					for(ss = 0; ss < w_num; ss++)
					{
						CCT_p[ss] = sump_s(ss,r);
						CCT_p[w_num + ss] = sump_b(ss,r);
					}

					// SCANG-O
					sump_o(0,r) = CCT_pval(CCT_p,CCT_weights);
					sump_o(0,r) = -log(sump_o(0,r));

					if(sump_o(0,r) > threshold(0,r))
					{
						threshold(0,r) = sump_o(0,r);
					}

					// SCANG-S
					sump_o(1,r) = CCT_pval(CCT_p(arma::span(0,w_num - 1)),CCT_weights(arma::span(0,w_num - 1)));
					sump_o(1,r) = -log(sump_o(1,r));

					if(sump_o(1,r) > threshold(1,r))
					{
						threshold(1,r) = sump_o(1,r);
					}

					// SCANG-B
					sump_o(2,r) = CCT_pval(CCT_p(arma::span(w_num,2*w_num - 1)),CCT_weights(arma::span(w_num,2*w_num - 1)));
					sump_o(2,r) = -log(sump_o(2,r));

					if(sump_o(2,r) > threshold(2,r))
					{
						threshold(2,r) = sump_o(2,r);
					}
				}
			}
		}

	}

	return threshold;
}





