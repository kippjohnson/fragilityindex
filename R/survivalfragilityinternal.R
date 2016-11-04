#' @title Survival analyysis fragility internal function
#' @description Compute the fragility of a coefficient in a log-rank survival test. This function is called by the wrapper function survivalfragility(), which replicates this function many times to obtain stable estimates.
#'
#' @param formula Model formula which will be evaluated by survdiff()
#' @param data Dataframe with values for model forma, passed to survdiff()
#' @param conf.level Significance level, set by default to 95\%
#'
#' @importFrom stringr str_extract
#' @importFrom survival survdiff
#' @importFrom survival Surv
#' @importFrom stats pchisq
#'
#' @examples
#' library(survival)
#' lung$status <- lung$status - 1 # recode status as a 0/1 variable (mandatory)
#'
#' survivalfragilityinternal(Surv(time, status) ~ pat.karno + strata(inst), data=lung)
#'
#' @return Returns the fragility index for a single run of a survival difference
#' @export survivalfragilityinternal

survivalfragilityinternal <- function(formula, data, conf.level=0.95){

  alpha <- (1 - conf.level)

  # Pass in data to initial survdiff() test
  formula=formula
  indata=data
  fragility.index <- 0

  m0 <- survdiff(formula, indata)
  summary(m0)
  model.pval <- ( 1 - pchisq(q=m0$chisq, df=nrow(m0$obs)-1) )

  # Extract the Surv object (i.e " Surv(time, event) ")
  so.str <- str_extract(as.character(formula)[2], pattern="^(.+?)\\)")
  so <- with(indata, eval(parse(text=so.str)))
  #print(paste("so.str:", so.str))
  #print(so)

  response.name <- names(which(sapply(indata, identical, so[1:nrow(indata),2])))
  #print(paste("response.name:", response.name))
  #print(paste("nchar(response.name)", nzchar(response.name)))

  a <- indata[,response.name]
  #print(a)
  #print(model.pval)

  if(model.pval>alpha){ # If covariate is not significant to start with...
    return(list(index=fragility.index))
  }else{

    while(model.pval<alpha){
      if(sum(levels(factor(a))==c("0","1"))==2){

        x <- sample(indata[,response.name], size=1)
        # print(c("x", x))

        if(x==0){
          # swap 0 for 1 in event column of indata
          index <- sample(which(indata[,response.name]==0),size=1)
          indata[,response.name][index] <- 1
        }

        if(x==1){
          # swap 1 for 0 in event column of indata
          index <- sample(which(indata[,response.name]==1),size=1)
          indata[,response.name][index] <- 0
        }

        # print(sum( indata[,response.name]  ))

        m1 <- survdiff(formula, indata)
        model.pval.new <- ( 1 - pchisq(q=m1$chisq, df=nrow(m1$obs)-1) )

        if(model.pval.new>model.pval){
          fragility.index <- fragility.index + 1
          model.pval <- model.pval.new
        }

        # print(c("0/1 Outcome:", model.pval, fragility.index))

      }else{
        stop("Something wrong with response variable")
      }
    }
  }
  return(list(index=fragility.index))
}

