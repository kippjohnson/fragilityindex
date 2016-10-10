#' @title Logistic fragility function
#' @description Compute the fragility of a coefficient in a logistic regression for dichotomous outcomes. Wrapper function for logisticfragilityinternal
#'
#' @param formula Model formula which will be evaluated by glm()
#' @param data Dataframe with values for model forma, passed to glm()
#' @param covariate Covariate name (string) whose fragility you would like to test
#' @param conf.level Significance level, set by default to 95%
#'
#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response
#' @examples
#' mydata <- read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")
#' mydata$rank <- factor(mydata$rank)
#' # Suggest much higher value for niter (1000+)
#' logisticfragility(admit ~ gre + gpa + rank, data = mydata, covariate="gre", niter=25)
#'
#' @return Returns the fragility index for a single run
#' @export logisticfragility


logisticfragility <- function(formula, data, covariate, niter){
  res <- replicate(niter, logisticfragilityinternal(formula, data=data, covariate=covariate)$index)
  #hist(res)
  mean.fragility <- mean(res)
  return(index=mean.fragility)
}
