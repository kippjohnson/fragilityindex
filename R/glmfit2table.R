#' @title glmfit2table
#' @description Return beta and p value for selected coefficent from glm.fit()
#'
#' @param a A glm.fit() object
#' @param covariate A covariate name
#'
#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response
#' @importFrom stats pt

#' @return Returns a given coefficent's beta and p value
#'
#' @export glmfit2table


glmfit2table <- function(a, covariate){
  # Compute P values and extract results from the GLM
  df.r <- a$df.residual
  coef.p <- a$coefficients
  covmat <- chol2inv(a$qr$qr)
  dimnames(covmat) <- list(names(coef.p), names(coef.p))
  var.cf <- diag(covmat)
  s.err <- sqrt(var.cf)
  tvalue <- coef.p/s.err
  dn <- c("Estimate", "Std. Error")
  pvalue <- 2 * pt(-abs(tvalue), df.r)
  coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
  dimnames(coef.table) <- list(names(coef.p), c(dn,"t value", "Pr(>|t|)"))
  #print(coef.table)
  model.beta <- coef.table[covariate,1]
  model.pval <- coef.table[covariate,4]
  return(list(beta=model.beta, pval=model.pval))
}

