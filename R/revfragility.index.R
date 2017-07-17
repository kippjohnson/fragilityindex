#' @title Reverse Fragility Index Calculation
#' @description Compute the reverse fragility index for a dichotomous outcome, i.e. the number of flipped cases it would take to make a non-significant result significant.
#'
#' @param intervention_event Number of events in intervention group
#' @param control_event Number of events in control group
#' @param intervention_n Total number of patients in intervention group
#' @param control_n Total number of patients in the control group
#' @param conf.level Significance level
#' @param verbose Logical indicating if function will return verbose results or only fragility index
#' @param print.mat Logical indicating if 2x2 matrices should be printed for each iteration of algorithm
#'
#' @examples
#' revfragility.index(6,5,50,50, verbose=TRUE, print.mat=FALSE)
#'
#' @return If verbose is FALSE, returns a list with fragility index. If
#' verbose is TRUE, returns a list with p-values for each fragility index
#' at each iteration of the algorithm.
#'
#' @importFrom stats fisher.test
#' @importFrom stats chisq.test
#' @export revfragility.index

revfragility.index <- function(intervention_event, control_event, intervention_n, control_n, conf.level=0.95, verbose=FALSE, print.mat=FALSE){

  if(control_event>intervention_event){

    tmp_event <- intervention_event
    tmp_n <- intervention_n

    intervention_event <- control_event
    intervention_n <- control_n

    control_event <- tmp_event
    control_n <- tmp_n

  }

  alpha <- (1 - conf.level)

  mat <- matrix(c(intervention_event, control_event, intervention_n-intervention_event, control_n-control_event),nrow=2)
  if(print.mat==TRUE){ print(mat) }
  fragility.index <- 0

  test <- fisher.test(mat)
  test2 <- chisq.test(mat)

  # Cells with 0 give Chi-square test a problem; if this is the case,
  # then just use the fisher test p-value
  if(is.na(test2$p.value)){test2$p.value <- test$p.value}

  if(verbose==FALSE){
    if(test$p.value<alpha | test2$p.value<alpha){
      return(list(index=fragility.index))

    }else{

      while(test$p.value > alpha){
        fragility.index <- fragility.index + 1

        if(control_event>0){ # cannot have negative events
          control_event = control_event - 1
        }

        intervention_event = intervention_event + 1
        if(print.mat==TRUE){ print(mat) }
        mat <- matrix(c(intervention_event, control_event, intervention_n-intervention_event, control_n-control_event),nrow=2)
        if(print.mat==TRUE){ print(mat) }
        test <- fisher.test(mat)

      }

      return(list(index=fragility.index))
    }
  }

  if(verbose==TRUE){
    outdf <- as.data.frame(matrix(NA, nrow=0, ncol=2))

    if(test$p.value<alpha | test2$p.value<alpha){
      res <- c(fragility.index, test$p.value)
      outdf <- rbind(outdf, res)
      return(outdf)

    }else{
      res <- c(fragility.index, test$p.value)
      outdf <- rbind(outdf, res)

      while(test$p.value > alpha){
        fragility.index <- fragility.index + 1

        if(control_event>1){ # cannot have negative events
          control_event = control_event - 1
        }

        intervention_event = intervention_event + 1

        mat <- matrix(c(intervention_event, control_event, intervention_n-intervention_event, control_n-control_event),nrow=2)
        if(print.mat==TRUE){ print(mat) }
        test <- fisher.test(mat)
        res <- c(fragility.index, test$p.value)
        outdf <- rbind(outdf, res)
      }
    }
  }
  names(outdf) <- c("reverse.fragility.index","p.value")
  outdf$p.value <- round(outdf$p.value, digits=3)
  return(outdf)
}

