#' @title Logistic fragility function
#' @description Compute the fragility of a coefficient in a logistic regression for dichotomous outcomes. Wrapper function for logisticfragilityinternal
#'
#' @param formula Model formula which will be evaluated by glm()
#' @param data Dataframe with values for model forma, passed to glm()
#' @param covariate Covariate name (string) whose fragility you would like to test
#' @param niter Number of iterations of algorithm to run
#' @param conf.level Significance level, set by default to 95\%
#' @param progress.bar Print a progress bar?

#' @importFrom stats glm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response
#' @importFrom graphics hist
#' @importFrom pbapply pbreplicate

#'
#' @examples
#' # Import and format example data
#' mydata <- read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")
#' mydata$rank <- factor(mydata$rank)
#'
#' # Suggest much higher value for niter in practice (1000+)
#' logisticfragility(admit ~ gre + gpa + rank, data = mydata, covariate="gre", niter=5)
#'
#' logisticfragility(admit ~ gre + gpa + rank, data = mydata, covariate="all", niter=5,
#'                  progress.bar = TRUE)
#'
#' @return Returns the fragility index for a single run
#' @export logisticfragility

logisticfragility <- function(formula, data, covariate, niter, conf.level=0.95, progress.bar=FALSE){

  if(covariate=="all"){
    covariate.names <- colnames(model.matrix(formula, data=data))
    result.store <- as.data.frame(matrix(NA, nrow=length(covariate.names), ncol=2))
    for(i in 1:length(covariate.names)){
        xi <- covariate.names[i]
        if(progress.bar==TRUE){
        print(paste0("Doing ", xi, "..."))
        res <- pbreplicate(niter, logisticfragilityinternal(formula=formula, data=data, covariate=xi, conf.level=conf.level)$index)
        }else{
        res <- replicate(niter, logisticfragilityinternal(formula=formula, data=data, covariate=xi, conf.level=conf.level)$index)
        }
        result.store[i, 1] <- xi
        result.store[i, 2] <- mean(res)
    }

  colnames(result.store) <- c("coefficient","fragility.index")
  return(result.store)

  }else{
  if(progress.bar==TRUE){
  res <- pbreplicate(niter, logisticfragilityinternal(formula=formula, data=data, covariate=covariate, conf.level=conf.level)$index)
  }else{
  res <- replicate(niter, logisticfragilityinternal(formula=formula, data=data, covariate=covariate, conf.level=conf.level)$index)
    }
    mean.fragility <- mean(res)
    return(fragility.index=mean.fragility)
  }
}


