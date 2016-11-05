#' @title Survival fragility function
#' @description Compute the fragility of a coefficient in a survival test. Wrapper function for logisticfragilityinternal. Uses the survdiff() function from the survival package.
#'
#' @param formula Model formula which will be evaluated by glm()
#' @param data Dataframe with values for model forma, passed to glm()
#' @param niter Number of iterations of algorithm to run
#' @param conf.level Significance level, set by default to 95\%
#' @param progress.bar Print a progress bar?
#' @importFrom stringr str_extract
#' @importFrom survival survdiff
#' @importFrom survival Surv
#' @importFrom graphics hist
#' @importFrom pbapply pbreplicate
#'
#' @examples
#' library(survival)
#' lung$status = lung$status - 1 # recode status as a 0/1 variable (mandatory)
#'
#' survivalfragility(Surv(time, status) ~ pat.karno + strata(inst),
#'                   data=lung, niter=10, progress.bar = FALSE)
#'
#' @return Returns the mean fragility index for a niter calculations
#' @export survivalfragility

survivalfragility <- function(formula, data, niter, conf.level=0.95, progress.bar=FALSE){

     if(progress.bar==TRUE){

        res <- pbreplicate(niter, survivalfragilityinternal(formula=formula, data=data, conf.level=conf.level))
        mean.fragility <- mean(unlist(res))
        return(fragility.index=mean.fragility)

     }else{

       res <- replicate(niter, survivalfragilityinternal(formula=formula, data=data, conf.level=conf.level))
       mean.fragility <- mean(unlist(res))
       return(index=mean.fragility)

     }
}
