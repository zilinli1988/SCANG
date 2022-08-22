#' Fit generalized linear model under the null hypothesis for unrelated samples.
#'
#' The \code{fit_null_glm_SCANG} function is a wrapper of the \code{\link{glm}} function from the
#' \code{\link{stats}} package that fits a regression model under the null hypothesis for
#' unrelated samples, which provides the preliminary step for subsequent
#' variant-set tests in whole
#' genome sequencing data analysis.
#' @param fixed an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the fixed effects model to be fitted.
#' @param data a data frame or list (or object coercible by \code{\link{as.data.frame}} to a data frame)
#' containing the variables in the model.
#' @param family a description of the error distribution and link function to be used
#' in the model. This can be a character string naming a family function, a family
#' function or the result of a call to a family function. (See \code{\link{family}} for details of family functions).
#' Can be either "gaussian" for continuous phenotype or "binomial" for binary phenotype.
#' @param times a number of pseudo-residuals (default = 2000).
#' @param ... additional arguments that could be passed to \code{\link{glm}}.
#' @return The function returns an object of the model fit from \code{\link{glm}} (\code{obj_nullmodel}),
#' with an additional element indicating the samples are unrelated
#' (\code{obj_nullmodel$relatedness = FALSE}). See \code{\link{glm}} for more details.
#' @export

fit_null_glm_SCANG <- function(fixed, data, family = binomial(link = "logit"), times = 2000, ...){

  obj_nullmodel <- glm(formula = fixed, data = data, family = family, ...)

  obj_nullmodel$relatedness <- FALSE
  obj_nullmodel$times <- times
  return(obj_nullmodel)
}


