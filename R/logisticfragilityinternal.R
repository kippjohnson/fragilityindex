#' @title Logistic fragility internal function
#' @description Compute the fragility of a coefficient in a logistic regression for dichotomous outcomes. This function is called by the wrapper function logistic.fragility(), which replicates this function many times to obtain more stable estimates.
#'
#' @param formula Model formula which will be evaluated by glm()
#' @param data Dataframe with values for model forma, passed to glm()
#' @param covariate Covarite name (string) whose fragility you would like to test
#' @param conf.level Significance level, set by default to 95\%
#'
#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response
#' @importFrom stats binomial
#' @importFrom stats binomial
#'
#' @examples
#' mydata <- read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")
#' mydata$rank <- factor(mydata$rank)
#' logisticfragilityinternal(admit ~ gre + gpa + rank, data = mydata, covariate="gre")
#'
#' @return Returns the fragility index for a single run
#' @export logisticfragilityinternal

logisticfragilityinternal <- function(formula, data, covariate, conf.level=0.95){

  alpha <- (1 - conf.level)

  # Pass in data to initial glm.fit
  formula=formula
  indata=data

  y0 <- model.response(model.frame(formula=formula, data = indata))
  x0 <- model.matrix(formula, data=indata)

  if(dim(table(y0))!=2){stop('Need a binary outcome')}

  fragility.index <- 0
  m0 <- glm.fit(x=x0, y=y0, family=binomial())
  model.beta <- glmfit2table(m0, covariate)$beta
  model.pval <- glmfit2table(m0, covariate)$pval

  if(model.pval>alpha){ # If covariate is not significant to start...
    return(list(index=fragility.index))
  }else{

    while(model.pval<alpha){
      # Case 1: Beta > 0 --> switch 1s to 0s
      if(model.beta>0){
        y0[sample(which(y0==1), size=1)] <- 0
        m0 <- glm.fit(x=x0, y=y0, family=binomial())
        model.pval <- glmfit2table(m0, covariate)$pval
        fragility.index <- fragility.index + 1
      }
      # Case 2: Beta < 0 --> switch 0s to 1s
      if(model.beta<0){
        y0[sample(which(y0==0), size=1)] <- 1
        m0 <- glm.fit(x=x0, y=y0, family=binomial())
        model.pval <- glmfit2table(m0, covariate)$pval
        fragility.index <- fragility.index + 1
      }
    }
    }
  return(list(index=fragility.index))
}


#logisticfragilityinternal(admit ~ gpa + gre + rank, data = mydata, covariate="all")
#logisticfragilityinternal(admit ~ gpa + gre + rank, data = mydata, covariate="gre")
#logisticfragilityinternal(admit ~ gpa + gre + rank, data = mydata, covariate="gpa")

