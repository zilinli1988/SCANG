#' SCANG procedure using omnibus test
#'
#' The \code{SCANG} function takes in genotype and the object from fitting the null
#' model and detect the association between a
#' quantitative/dichotomous phenotype and a variant-set in a sequence by using SCANG procedure, including SCANG-O, SCANG-B and SCANG-S.
#' For each region, the scan statistic of SCANG-O is the set-based p-value of STAAR-O, which is an omnibus test that aggregated p-values
#' across different types of multiple annotation-weighted variant-set tests SKAT(1,1), SKAT(1,25), Burden(1,1) and Burden(1,25) using ACAT method;
#' the scan statistic of SCANG-S is the set-based p-value of STAAR-S, which is an omnibus test that aggregated p-values
#' across multiple annotation-weighted variant-set tests SKAT(1,1) and SKAT(1,25) using ACAT method;
#' the scan statistic of SCANG-B is the set-based p-value of STAAR-B, which is an omnibus test that aggregated p-values
#' across multiple annotation-weighted variant-set tests Burden(1,1) and Burden(1,25) using ACAT method.
#' @param genotype an n*p genotype matrix (dosage matrix) of the target sequence, where n is the sample
#' size and p is the number of variants.
#' @param obj_nullmodel an object from fitting the null model, which is the
#' output from either \code{\link{fit_null_glm_SCANG}} function for unrelated samples or
#' \code{\link{fit_null_glmmkin_SCANG}} function for related samples. Note that \code{\link{fit_null_glmmkin_SCANG}}
#' is a wrapper of \code{\link{glmmkin}} function from the \code{\link{GMMAT}} package.
#' @param Lmax maximum number of variants in searching windows.
#' @param Lmin minimum number of variants in searching windows.
#' @param annotation_phred a data frame or matrix of functional annotation data
#' of dimension p*q (or a vector of a single annotation score with length p).
#' Continuous scores should be given in PHRED score scale, where the PHRED score
#' of j-th variant is defined to be -10*log10(rank(-score_j)/total) across the genome. (Binary)
#' categorical scores should be taking values 0 or 1, where 1 is functional and 0 is
#' non-functional. If not provided, SCANG will perform the original procedure without annotations (default = NULL).
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.05).
#' @param steplength difference of number of variants in searching windows, that is, the number of variants in
#' searching windows are Lmin, Lmin+steplength, Lmin+steplength,...,Lmax (default = 5).
#' @param alpha familty-wise/genome-wide significance level (default = 0.05).
#' @param filter a filtering threshold of screening method for SKAT. SKAT p-values are calculated for regions whose p-value
#' is possibly smaller than the filtering threshold (default = 1e-4).
#' @param f an overlap fraction, which controls for the overlapping proportion of of detected regions. For example,
#' when f=0, the detected regions are non-overlapped with each other,
#' and when f=1, we keep every susceptive region as detected regions (default = 0.5).
#' @param subseq_num number of variants run in each sub-sequence (default = 2000).
#' @return The function returns a list with the following members:
#' @return \code{SCANG_O_res}:   A matrix that summarized the significant region detected by SCANG-O.
#' The first column is the -log(p-value) of the detected region.
#' The next two columns are the location of the detected region (in sense of variants order).
#' The last column is the family-wise/genome-wide error rate of the detected region.
#' The result (0,0,0,1) means there is no significant region.
#' @return \code{SCANG_O_top1}:  A vector of length 4 which summarized the top 1 region detected by SCANG-O.
#' The first element is the -log(p-value) of the region.
#' The next two elements are the location of the detected region (in sense of variants order).
#' The last element is the family-wise/genome-wide p-value.
#' @return \code{SCANG_O_thres}: Empirical threshold of SCANG-O for controlling the family-wise type I error at alpha level.
#' @return \code{SCANG_O_thres_boot}: A vector of Monte Carlo simulation sample for generating the empirical threshold. The 1-alpha quantile of this vector is
#' the empirical threshold.
#' @return \code{SCANG_S_res, SCANG_S_thres, SCANG_S_top1, SCANG_S_thres_boot}: Analysis
#' results using SCANG-S. Details see SCANG-O.
#' @return \code{SCANG_B_res, SCANG_B_thres, SCANG_B_top1, SCANG_B_thres_boot}: Analysis
#' results using SCANG-B. Details see SCANG-O.
#' @references Li, Z., et al. (2019). Dynamic Scan Procedure for Detecting Rare-Variant
#' Association Regions in Whole-Genome Sequencing Studies.
#' \emph{The American Journal of Human Genetics 104}(5), 802-814.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/30982610/}{pub})
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in-silico functional annotations empowers rare variant association analysis of
#' large whole genome sequencing studies at scale.
#' \emph{Nature Genetics (in press)}.
#' @references Liu, Y., et al. (2019). Acat: A fast and powerful p value combination
#' method for rare-variant analysis in sequencing studies.
#' \emph{The American Journal of Human Genetics 104}(3), 410-421.
#' (\href{https://www.sciencedirect.com/science/article/pii/S0002929719300023}{pub})
#' @export

SCANG <- function(genotype,obj_nullmodel,Lmin,Lmax,annotation_phred=NULL,rare_maf_cutoff=0.05,steplength=5,alpha=0.05,filter=1e-4,f=0.5,subseq_num=2000)
{

	seed <- 666

	samplesize <- dim(genotype)[1]
	folds <- floor(dim(genotype)[2]/subseq_num)
	folds <- max(1,folds)
	times <- obj_nullmodel$times

	if(class(genotype)=="dgCMatrix")
	{
	  MAF <- colMeans(genotype)/2
	  genotype <- genotype[,(MAF<rare_maf_cutoff)&(MAF>0)]
	  RV_label <- (MAF<rare_maf_cutoff)&(MAF>0)
	  MAF <- MAF[(MAF<rare_maf_cutoff)&(MAF>0)]
	}else
	{
	  if(is.matrix(genotype) == FALSE){
	    stop(paste0("Number of rare variant in the set is less than 2!"))
	  }

	  genotype <- matrix_flip(genotype)
	  MAF <- genotype$MAF
	  genotype <- genotype$Geno[,(MAF<rare_maf_cutoff)&(MAF>0)]
	  RV_label <- (MAF<rare_maf_cutoff)&(MAF>0)
	  MAF <- MAF[(MAF<rare_maf_cutoff)&(MAF>0)]
	  genotype <- matsp(genotype)$Geno
	}


	## beta(1,25)
    w_1 <- dbeta(MAF,1,25)
    ## beta(1,1)
    w_2 <- dbeta(MAF,1,1)
    if(length(annotation_phred) == 0)
	{
      ## Burden, SKAT
      w_B <- as.matrix(cbind(w_1,w_2))
	  w_S <- as.matrix(cbind(w_1,w_2))
    }else
	{
	  annotation_phred <- annotation_phred[RV_label,,drop=FALSE]
	  annotation_rank <- 1 - 10^(-annotation_phred/10)

      ## Burden
      w_B_1 <- annotation_rank*w_1
      w_B_1 <- cbind(w_1,w_B_1)
      w_B_2 <- annotation_rank*w_2
      w_B_2 <- cbind(w_2,w_B_2)
      w_B <- cbind(w_B_1,w_B_2)
      w_B <- as.matrix(w_B)

      ## SKAT
      w_S_1 <- sqrt(annotation_rank)*w_1
      w_S_1 <- cbind(w_1,w_S_1)
      w_S_2 <- sqrt(annotation_rank)*w_2
      w_S_2 <- cbind(w_2,w_S_2)
      w_S <- cbind(w_S_1,w_S_2)
      w_S <- as.matrix(w_S)
    }


	## SCANG-O
	L20_O <- matrix(0,folds,times)
	res_O <- c()
	resmost_O <- c()

	## SCANG-S
	L20_S <- matrix(0,folds,times)
	res_S <- c()
	resmost_S <- c()


	## SCANG-B
	L20_B <- matrix(0,folds,times)
	res_B <- c()
	resmost_B <- c()



	subnum <- floor(dim(genotype)[2]/folds)

	for(i in 1:folds)
	{
		if(i<folds)
		{
			genotypesub <- genotype[,(subnum*(i-1)+1):(i*subnum+Lmax)]
			weightssub_B <- w_B[(subnum*(i-1)+1):(i*subnum+Lmax),]
			weightssub_S <- w_S[(subnum*(i-1)+1):(i*subnum+Lmax),]
		}
		if(i==folds)
		{
			genotypesub <- genotype[,(subnum*(i-1)+1):dim(genotype)[2]]
			weightssub_B <- w_B[(subnum*(i-1)+1):dim(genotype)[2],]
			weightssub_S <- w_S[(subnum*(i-1)+1):dim(genotype)[2],]
		}

		if(obj_nullmodel$relatedness)
		{
			if(!obj_nullmodel$sparse_kins)
			{
				residuals.phenotype <- obj_nullmodel$scaled.residuals

				set.seed(19880615+seed)
				threstemp <- SCANG_O_Thres_Relatedness(genotypesub,obj_nullmodel$P,obj_nullmodel$pseudo_residuals,times,Lmax,Lmin,steplength,weightssub_B,weightssub_S,filter)

				begid <- subnum*(i-1)+1

				##### SCANG-O
				L20_O[i,] <- threstemp[1,]

				emL20_O <- apply(L20_O,2,max)
				th0_O <- quantile(emL20_O,1-alpha)
				if(th0_O < -log(filter))
				{
					th0_O = -log(filter)
				}

				restemp <- SCANG_O_Search_Relatedness(genotypesub,obj_nullmodel$P,residuals.phenotype,th0_O,Lmax,Lmin,steplength,weightssub_B,weightssub_S,begid,filter,f)
				res_O <- rbind(res_O,restemp$res)
				resmost_O <- rbind(resmost_O,restemp$resmost)


				##### SCANG-S
				L20_S[i,] <- threstemp[2,]

				emL20_S <- apply(L20_S,2,max)
				th0_S <- quantile(emL20_S,1-alpha)
				if(th0_S < -log(filter))
				{
					th0_S = -log(filter)
				}

				restemp <- SCANG_S_Search_Relatedness(genotypesub,obj_nullmodel$P,residuals.phenotype,th0_S,Lmax,Lmin,steplength,weightssub_S,begid,filter,f)
				res_S <- rbind(res_S,restemp$res)
				resmost_S <- rbind(resmost_S,restemp$resmost)

				##### SCANG-B
				L20_B[i,] <- threstemp[3,]

				emL20_B <- apply(L20_B,2,max)
				th0_B <- quantile(emL20_B,1-alpha)
				if(th0_B < -log(filter))
				{
					th0_B = -log(filter)
				}

				restemp <- SCANG_B_Search_Relatedness(genotypesub,obj_nullmodel$P,residuals.phenotype,th0_B,Lmax,Lmin,steplength,weightssub_B,begid,filter,f)
				res_B <- rbind(res_B,restemp$res)
				resmost_B <- rbind(resmost_B,restemp$resmost)

			}else
			{
				Sigma_i <- obj_nullmodel$Sigma_i
				Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
				cov <- obj_nullmodel$cov

				residuals.phenotype <- obj_nullmodel$scaled.residuals

				set.seed(19880615+seed)
				threstemp <- SCANG_O_Thres_Relatedness_sp(genotypesub,Sigma_i,Sigma_iX,cov,obj_nullmodel$pseudo_residuals,times,Lmax,Lmin,steplength,weightssub_B,weightssub_S,filter)

				begid <- subnum*(i-1)+1

				##### SCANG-O-threshold
				L20_O[i,] <- threstemp[1,]

				emL20_O <- apply(L20_O,2,max)
				th0_O <- quantile(emL20_O,1-alpha)
				if(th0_O < -log(filter))
				{
					th0_O = -log(filter)
				}
				##### SCANG-S-threshold
				L20_S[i,] <- threstemp[2,]

				emL20_S <- apply(L20_S,2,max)
				th0_S <- quantile(emL20_S,1-alpha)
				if(th0_S < -log(filter))
				{
					th0_S = -log(filter)
				}
				##### SCANG-B-threshold
				L20_B[i,] <- threstemp[3,]

				emL20_B <- apply(L20_B,2,max)
				th0_B <- quantile(emL20_B,1-alpha)
				if(th0_B < -log(filter))
				{
					th0_B = -log(filter)
				}

				##### SCANG-O
				restemp <- SCANG_O_Search_Relatedness_sp(genotypesub,Sigma_i,Sigma_iX,cov,residuals.phenotype,th0_O,th0_S,th0_B,Lmax,Lmin,steplength,weightssub_B,weightssub_S,begid,filter,f)
				res_O <- rbind(res_O,restemp$res_o)
				resmost_O <- rbind(resmost_O,restemp$resmost_o)

				##### SCANG-S
				res_S <- rbind(res_S,restemp$res_s)
				resmost_S <- rbind(resmost_S,restemp$resmost_s)


				## SCANG-B
				res_B <- rbind(res_B,restemp$res_b)
				resmost_B <- rbind(resmost_B,restemp$resmost_b)
			}
		}else
		{
			### unrelated samples
			X <- model.matrix(obj_nullmodel)

			working <- obj_nullmodel$weights
			sigma <- sqrt(summary(obj_nullmodel)$dispersion)
			residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values

			if(obj_nullmodel$family[1] == "binomial")
			{
				fam <- 1
			}
			if(obj_nullmodel$family[1] == "gaussian")
			{
				fam <- 0
			}

			set.seed(19880615+seed)
			threstemp <- SCANG_O_Thres(genotypesub,X,working,sigma,fam,obj_nullmodel$times,Lmax,Lmin,steplength,weightssub_B,weightssub_S,filter)

			begid <- subnum*(i-1)+1

			##### SCANG-O
			L20_O[i,] <- threstemp[1,]

			emL20_O <- apply(L20_O,2,max)
			th0_O <- quantile(emL20_O,1-alpha)
			if(th0_O < -log(filter))
			{
				th0_O = -log(filter)
			}

			restemp <- SCANG_O_Search(genotypesub,X,working,sigma,fam,residuals.phenotype,th0_O,Lmax,Lmin,steplength,weightssub_B,weightssub_S,begid,filter,f)
			res_O <- rbind(res_O,restemp$res)
			resmost_O <- rbind(resmost_O,restemp$resmost)

			##### SCANG-S
			L20_S[i,] <- threstemp[2,]

			emL20_S <- apply(L20_S,2,max)
			th0_S <- quantile(emL20_S,1-alpha)
			if(th0_S < -log(filter))
			{
				th0_S = -log(filter)
			}

			restemp <- SCANG_S_Search(genotypesub,X,working,sigma,fam,residuals.phenotype,th0_S,Lmax,Lmin,steplength,weightssub_S,begid,filter,f)
			res_S <- rbind(res_S,restemp$res)
			resmost_S <- rbind(resmost_S,restemp$resmost)

			## SCANG-B
			L20_B[i,] <- threstemp[3,]

			emL20_B <- apply(L20_B,2,max)
			th0_B <- quantile(emL20_B,1-alpha)
			if(th0_B < -log(filter))
			{
				th0_B = -log(filter)
			}

			restemp <- SCANG_B_Search(genotypesub,X,working,sigma,fam,residuals.phenotype,th0_B,Lmax,Lmin,steplength,weightssub_B,begid,filter,f)
			res_B <- rbind(res_B,restemp$res)
			resmost_B <- rbind(resmost_B,restemp$resmost)
		}
	}

	rm(L20_O)
	gc()

	rm(L20_S)
	gc()

	rm(L20_B)
	gc()


	## SCANG-O
	## results
	res_O <- res_O[res_O[,1]>th0_O,]
	if(length(res_O)==0)
	{
		res_O <- c(0,0,0,0)
	}
	if(length(res_O)>4)
	{
		res_O <- regionfilter(res_O,f)
	}

	## calcualte empirical p-value
	if(length(res_O)==4)
	{
	  res_O[4] <- sum(emL20_O>res_O[1])/length(emL20_O)
	  res_O <- matrix(res_O,nrow=1)
	}
	if(length(res_O)!=4)
	{
	  res_O[,4] <- apply(res_O,1,function(z) sum(emL20_O>z[1])/length(emL20_O))
	}

	## Top 1 region
	mostnum <- which.max(resmost_O[,1])
	resmost_O <- resmost_O[mostnum,]
	resmost_O[4] <- sum(emL20_O>resmost_O[1])/length(emL20_O)


	## SCANG-S
	## results
	res_S <- res_S[res_S[,1]>th0_S,]
	if(length(res_S)==0)
	{
		res_S <- c(0,0,0,0)
	}
	if(length(res_S)>4)
	{
		res_S <- regionfilter(res_S,f)
	}

	## calcualte empirical p-value
	if(length(res_S)==4)
	{
	  res_S[4] <- sum(emL20_S>res_S[1])/length(emL20_S)
	  res_S <- matrix(res_S,nrow=1)
	}
	if(length(res_S)!=4)
	{
	  res_S[,4] <- apply(res_S,1,function(z) sum(emL20_S>z[1])/length(emL20_S))
	}

    ## Top 1 region
	mostnum <- which.max(resmost_S[,1])
	resmost_S <- resmost_S[mostnum,]
	resmost_S[4] <- sum(emL20_S>resmost_S[1])/length(emL20_S)


	## SCANG-B
	res_B <- res_B[res_B[,1]>th0_B,]
	if(length(res_B)==0)
	{
		res_B <- c(0,0,0,0)
	}
	if(length(res_B)>4)
	{
		res_B <- regionfilter(res_B,f)
	}

	## calcualte empirical p-value
	if(length(res_B)==4)
	{
	  res_B[4] <- sum(emL20_B>res_B[1])/length(emL20_B)
	  res_B <- matrix(res_B,nrow=1)
	}
	if(length(res_B)!=4)
	{
	  res_B[,4] <- apply(res_B,1,function(z) sum(emL20_B>z[1])/length(emL20_B))
	}

	## Top 1 region
	mostnum <- which.max(resmost_B[,1])
	resmost_B <- resmost_B[mostnum,]
	resmost_B[4] <- sum(emL20_B>resmost_B[1])/length(emL20_B)


	Lst <- list(SCANG_O_res=res_O,SCANG_O_top1=resmost_O,SCANG_O_thres=th0_O,SCANG_O_thres_boot=emL20_O,
	SCANG_S_res=res_S,SCANG_S_top1=resmost_S,SCANG_S_thres=th0_S,SCANG_S_thres_boot=emL20_S,
	SCANG_B_res=res_B,SCANG_B_top1=resmost_B,SCANG_B_thres=th0_B,SCANG_B_thres_boot=emL20_B,
	RV_label = RV_label)

	return(Lst)
}
