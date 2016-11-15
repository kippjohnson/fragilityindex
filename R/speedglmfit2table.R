speedglmfit2table <- function(a, covariate){
  # Compute P values and extract results from the GLM
  df.r <- a$df
  coef.p <- a$coefficients
  #covmat <- chol2inv(a$XTX)
  covmat <- chol2inv(qr(m0$XTX)$qr)
  #print(head(covmat))
  dimnames(covmat) <- list(names(coef.p), names(coef.p))
  var.cf <- diag(covmat)
  #print(head(var.cf))
  s.err <- sqrt(var.cf)
  #print(s.err)
  tvalue <- coef.p/s.err
  #print(tvalue)
  dn <- c("Estimate", "Std. Error")
  #print(abs(tvalue))
  #print(df.r)
  pvalue <- 2 * pt(-abs(tvalue), df.r)
  coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
  dimnames(coef.table) <- list(names(coef.p), c(dn,"t value", "Pr(>|t|)"))
  #print(coef.table)
  model.beta <- coef.table[covariate,1]
  model.pval <- coef.table[covariate,4]
  return(list(beta=model.beta, pval=model.pval))
}


out <- speedglmfit2table(m0)

sort(round(speedglmfit2table(m0)$pval,2))
