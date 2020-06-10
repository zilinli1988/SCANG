#' Fitting generalized linear mixed model with known relationship matrices
#' under the null hypothesis for related samples.
#'
#' The \code{fit_null_glmmkin_SCANG} function is a wrapper of the \code{\link{glmmkin}} function from
#' the \code{\link{GMMAT}} package that fits a regression model under the null hypothesis
#' for related samples, which provides the preliminary step for subsequent
#' variant-set tests in whole genome sequencing data analysis. More details
#' see \code{\link{glmmkin}}.
#' @param fixed an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the fixed effects model to be fitted.
#' @param data a data frame or list (or object coercible by \code{\link{as.data.frame}} to a data frame)
#' containing the variables in the model.
#' @param kins a known positive semi-definite relationship matrix
#' (e.g. kinship matrix in genetic association studies) or a list of known
#' positive semi-definite relationship matrices. The rownames and colnames of
#' these matrices must at least include all samples as specified in the \code{id} column
#' of the data frame \code{data}.
#' @param use_sparse a logical switch of whether the provided dense \code{kins} matrix should be
#' transformed to a sparse matrix (default = NULL).
#' @param kins_cutoff the cutoff of setting all entries with smaller values to 0 in \code{kins} matrix
#' (default = 0.022).
#' @param id a column in the data frame \code{data}, indicating the id of samples.
#' When there are duplicates in \code{id}, the data is assumed to be longitudinal with repeated measures.
#' @param random.slope an optional column indicating the random slope for time effect used
#' in a mixed effects model for longitudinal data. It must be included in the names of \code{data}.
#' There must be duplicates in \code{id} and \code{method.optim} must be "AI" (default = NULL).
#' @param groups an optional categorical variable indicating the groups used in a
#' heteroscedastic linear mixed model (allowing residual variances in different groups to be different).
#' This variable must be included in the names of \code{data}, and \code{family} must be "gaussian"
#' and \code{method.optim} must be "AI" (default = NULL).
#' @param family a description of the error distribution and link function to be used
#' in the model. This can be a character string naming a family function, a family
#' function or the result of a call to a family function. (See \code{\link{family}} for details of family functions).
#' @param method method of fitting the generalized linear mixed model. Either "REML" or "ML" (default = "REML").
#' @param method.optim optimization method of fitting the generalized linear mixed model.
#' Either "AI", "Brent" or "Nelder-Mead" (default = "AI").
#' @param maxiter a positive integer specifying the maximum number of iterations when
#' fitting the generalized linear mixed model (default = 500).
#' @param tol a positive number specifying tolerance, the difference threshold for parameter
#' estimates below which iterations should be stopped (default = 1e-5).
#' @param taumin the lower bound of search space for the variance component parameter \eqn{\tau} (default = 1e-5),
#' used when \code{method.optim} = "Brent". See Details.
#' @param taumax the upper bound of search space for the variance component parameter \eqn{\tau} (default = 1e5),
#' used when \code{method.optim} = "Brent". See Details.
#' @param tauregion the number of search intervals for the REML or ML estimate of the variance component
#' parameter \eqn{\tau} (default = 10), used when \code{method.optim} = "Brent". See Details.
#' @param times a number of pesudo-residuals (default = 2000).
#' @param verbose a logical switch for printing detailed information (parameter estimates in each iteration)
#' for testing and debugging purpose (default = FALSE).
#' @param ... additional arguments that could be passed to \code{\link{glm}}.
#' @return The function returns an object of the model fit from \code{\link{glmmkin}} (\code{obj_nullmodel}),
#' with additional elements indicating the samples are related (\code{obj_nullmodel$relatedness = TRUE}),
#' whether the \code{kins} matrix is sparse when fitting the null model, and the matrix of pseudo residuals. See \code{\link{glmmkin}} for more details.
#' @references Chen, H. et al. (2016). Control for population structure and relatedness for binary traits
#' in genetic association studies via logistic mixed models. \emph{American Journal of Human Genetics 98(4), 653-666}.
#' (\href{https://www.sciencedirect.com/science/article/pii/S000292971600063X}{pub})
#' @references Chen, H. et al. (2019). Efficient variant set mixed model association tests for continuous and
#' binary traits in large-scale whole genome sequencing studies. \emph{American Journal of Human Genetics 104(2), 260-274}.
#' (\href{https://www.sciencedirect.com/science/article/pii/S0002929718304658}{pub})
#' @references Chen, H. & Conomos, M.P. (2019). GMMAT-package: Generalized Linear Mixed Model Association Tests.
#' (\href{https://rdrr.io/github/hanchenphd/GMMAT/man/GMMAT-package.html}{web})
#' @export

fit_null_glmmkin_SCANG <- function(fixed, data = parent.frame(), kins, use_sparse = NULL,
                             kins_cutoff = 0.022, id, random.slope = NULL, groups = NULL,
                             family = binomial(link = "logit"), method = "REML",
                             method.optim = "AI", maxiter = 500, tol = 1e-5,
                             taumin = 1e-5, taumax = 1e5, tauregion = 10, times = 2000,
                             verbose = FALSE, ...){
  if(class(kins) != "matrix" && !(!is.null(attr(class(kins), "package")) && attr(class(kins), "package") == "Matrix")){
    stop("kins is not a matrix!")
  }
  else if(!is.null(attr(class(kins), "package")) && attr(class(kins), "package") == "Matrix"){
    print("kins is a sparse matrix.")
    obj_nullmodel <- glmmkin(fixed = fixed, data = data, kins = kins, id = id,
                             random.slope = random.slope, groups = groups,
                             family = family, method = method,
                             method.optim = method.optim, maxiter = maxiter,
                             tol = tol, taumin = taumin, taumax = taumax,
                             tauregion = tauregion, verbose = verbose, ...)
	obj_nullmodel$sparse_kins <- TRUE

	Sigma_i <- obj_nullmodel$Sigma_i
	Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
	cov <- obj_nullmodel$cov

	R <- chol(Sigma_i)
	S <- chol(cov)
	R_inverse <- solve(R)

	K <- crossprod(t(S),crossprod(Sigma_iX,R_inverse))
	KS <- svd_c(as.matrix(t(K)))
	rm(K)

	eigen_diff <- rep(1,dim(Sigma_i)[1]) - c((KS$s)^2,rep(0,dim(Sigma_i)[1]-length(KS$s)))
	eigen_diff[abs(eigen_diff)<1e-10] <- 0
	eigen_diff <- sqrt(eigen_diff)

	### residuals
	set.seed(19880615+666)
	y <- matrix(rnorm(times*dim(KS$U)[2]),times,dim(KS$U)[2])
	y <- y%*%diag(eigen_diff)
	ytU <- tcrossprod(y,KS$U)

	rm(y)

	obj_nullmodel$pseudo_residuals <- ytU%*%R
	obj_nullmodel$pseudo_residuals <- as.matrix(obj_nullmodel$pseudo_residuals)

  }else if(!is.null(use_sparse) && use_sparse){
    print(paste0("kins is a dense matrix, transforming it into a sparse matrix using cutoff ", kins_cutoff,"."))
    kins_sp <- makeSparseMatrix(kins, thresh = kins_cutoff)
    if(class(kins_sp) == "dsyMatrix" || kins_cutoff <= min(kins)){
      stop(paste0("kins is still a dense matrix using cutoff ", kins_cutoff,". Please try a larger kins_cutoff or use_sparse = FALSE"))
    }
    rm(kins)
    obj_nullmodel <- glmmkin(fixed = fixed, data = data, kins = kins_sp, id = id,
                             random.slope = random.slope, groups = groups,
                             family = family, method = method,
                             method.optim = method.optim, maxiter = maxiter,
                             tol = tol, taumin = taumin, taumax = taumax,
                             tauregion = tauregion, verbose = verbose, ...)
    obj_nullmodel$sparse_kins <- TRUE

	Sigma_i <- obj_nullmodel$Sigma_i
	Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
	cov <- obj_nullmodel$cov

	R <- chol(Sigma_i)
	S <- chol(cov)
	R_inverse <- solve(R)

	K <- crossprod(t(S),crossprod(Sigma_iX,R_inverse))
	KS <- svd_c(as.matrix(t(K)))
	rm(K)

	eigen_diff <- rep(1,dim(Sigma_i)[1]) - c((KS$s)^2,rep(0,dim(Sigma_i)[1]-length(KS$s)))
	eigen_diff[abs(eigen_diff)<1e-10] <- 0
	eigen_diff <- sqrt(eigen_diff)

	### residuals
	set.seed(19880615+666)
	y <- matrix(rnorm(times*dim(KS$U)[2]),times,dim(KS$U)[2])
	y <- y%*%diag(eigen_diff)
	ytU <- tcrossprod(y,KS$U)

	rm(y)

	obj_nullmodel$pseudo_residuals <- ytU%*%R
	obj_nullmodel$pseudo_residuals <- as.matrix(obj_nullmodel$pseudo_residuals)

  }else{
    print("kins is a dense matrix.")
    obj_nullmodel <- glmmkin(fixed = fixed, data = data, kins = kins, id = id,
                             random.slope = random.slope, groups = groups,
                             family = family, method = method,
                             method.optim = method.optim, maxiter = maxiter,
                             tol = tol, taumin = taumin, taumax = taumax,
                             tauregion = tauregion, verbose = verbose, ...)
	  obj_nullmodel$sparse_kins <- FALSE
	  P_eigen <- eigen(obj_nullmodel$P)
	  P_eigen$values[abs(P_eigen$values) < 1e-12] <- 0
	  Omega <- sqrt(P_eigen$values)*t(P_eigen$vectors)

	  set.seed(19880615+666)
	  y <- matrix(rnorm(times*dim(Omega)[1]),times,dim(Omega)[1])
	  obj_nullmodel$pseudo_residuals <- y%*%Omega
	  obj_nullmodel$pseudo_residuals <- as.matrix(obj_nullmodel$pseudo_residuals)

  }
  obj_nullmodel$times <- times
  obj_nullmodel$relatedness <- TRUE
  return(obj_nullmodel)
}


