#' @title Survival Fragility Function
#' @description Compute the fragility of a coefficient in a Cox P-H regression for survival analysis, i.e. the number of removed observations it would take to make a significant-result non-significant. Uses the coxph() function from the survival package.
#'
#' @param formula Model formula which will be evaluated by coxph()
#' @param data Dataframe with values for model forma, passed to coxph()
#' @param covariate Vector of covariates to find fragility index for. Default is all covariates in formula
#' @param conf.level Significance level
#' @param verbose Logical indicating if function will return verbose results or only fragility index
#'
#' @importFrom survival coxph
#' @importFrom survival Surv
#' @importFrom stats glm
#' @importFrom stats terms.formula
#' @importFrom stats terms
#' @importFrom stats update
#' @importFrom stats as.formula
#' @importFrom stats residuals
#' @importFrom stats anova
#'
#' @examples
#' library(survival); data <- lung
#' data$status = lung$status - 1 # recode status as a 0/1 variable
#'
#' survivalfragility(Surv(time, status) ~ pat.karno + strata(inst),
#'                   data, covariate = "pat.karno")
#'
#' survivalfragility(Surv(time, status) ~ pat.karno + ph.karno + strata(inst),
#'                   data, verbose = TRUE)
#'                   #algorithm does not converge for strata(inst)
#'
#' survivalfragility(Surv(time, status) ~ pat.karno + ph.karno + strata(inst),
#'                   data, covariate = c("pat.karno","ph.karno"))
#'
#'
#'
#'
#' @return If verbose is FALSE, returns a list with fragility indices for selected covariates. If
#' verbose is TRUE, returns a list with p-values for each fragility index
#' at each iteration of the algorithm.
#'
#' @export survivalfragility


#include in examples later. makes checking take a long time (approx 3 min)
#survivalfragility(Surv(futime,death)~sex+age,data=flchain)

survivalfragility <- function(formula, data, covariate = "all.factors.default", conf.level = 0.95, verbose = FALSE) {

  if ("all.factors.default" %in% covariate) {
    object <- terms.formula(formula)
    covariate.names <- attr(object, "term.labels")
  } else {
    covariate.names <- covariate
  }

  result.store <- vector("list", length(covariate.names))
  names(result.store) = covariate.names

  na.model <- coxph(formula,data)
  na.data <- na.model$na.action
  if (!is.null(na.data)){
    data <- data[-na.data,]
  }

  for (i in 1:length(covariate.names)) {
    if (verbose == TRUE) {
      result <- survivalfragilityinternal(formula, data, covariate.names[i], conf.level)
    } else{
      result <- survivalfragilityinternal(formula, data, covariate.names[i], conf.level)
      result <- result[1]
    }
    result.store[[paste(covariate.names[i])]] <- result
  }
  return(result.store)
  }

survivalfragilityinternal <- function(formula, data, covariate, conf.level) {

  alpha <- (1 - conf.level)
  model <- coxph(formula, data)
  nullmodel <- update(model, as.formula(paste(".~.-", covariate)))

  delta.resid <- residuals(model, type = "deviance") - residuals(nullmodel, type = "deviance") #which residuals to use?
  index <- c(1:length(delta.resid))

  response <- formula[[2]]
  if (length(response) == 4){
    event = response[[4]]
  } else{
    event = response[[3]]
  }
  y <- event

  ordering <- cbind(index, delta.resid, data[ ,paste(y)])

  ordering <- cbind(ordering, (ordering[ ,2] - ordering[ ,3] * 2 * ordering[ ,2])) #m1

  #ordering[,4] <- ordering[,2] #m2

  #ordering[,4] <- ordering[,2]*-1 #m3


  #ordering[,4] <- abs(ordering[,2]) #m4

  #m5
  #resid.ratio <- residuals(nullmodel,type = "martingale") / residuals(model,type="martingale")

  #ordering[,4] <- resid.ratio

  #m6
  #resid.ratio <- residuals(nullmodel,type="deviance") / residuals(model,type = "deviance") * -1
  #ordering[,4] <- resid.ratio


  #ordering <- ordering[order(-ordering[ ,4]), ]
  #ordering <- ordering[sample(nrow(ordering)),]

  #m7 log parallel
  #prediction <- predict(model,type="expected")
  #nullprediction <- predict(nullmodel,type="expected")
  #residual <- data[,paste(event)] - exp(-prediction)
  #nullresid <- data[,paste(event)] - exp(-nullprediction)

  #delta.resid <- residual - nullresid
  #response <- formula[[2]]

  #if (length(response) == 4){
  #  event = response[[4]]
  #} else{
  #  event = response[[3]]
  #}
  #y <- event
  #ordering <- cbind(index, delta.resid, data[ ,paste(y)])
  #ordering <- cbind(ordering, (ordering[ ,2] - ordering[ ,3]*2*ordering[ ,2]))
  #ordering <- cbind(ordering,abs(ordering[,2]))
  ordering <- ordering[order(-ordering[ ,4]), ]

  #ordering <- ordering[order(-ordering[ ,4]), ]

  # status = all.vars(formula(model))[2]



  pval <- anova(model, nullmodel)$'P(>|Chi|)'[2]
  if (is.na(pval)) {
    return(list(fragility.index = 0, point.diagnostics =
                "No points removed. Covariate already not significant at confidence level" ))
  }
  index <- 1
  indices <- c()
  fragility.index <- 0
  iter <- 0

  pvalues <- pval
  not.significant <- FALSE
  if (pval > alpha) {
    not.significant <- TRUE
  }

  while (pval <= alpha & (iter <= nrow(data) - 3 * length(terms(formula)))) {

    indices.new <- c(indices,index)
    point <- ordering[indices.new,1]
    modified.data <- data[-point, ]
    newmodel <- coxph(formula, modified.data)
    newnullmodel <- update(newmodel, as.formula(paste(".~.-", covariate)))
    pval.new <- anova(newmodel, newnullmodel)$'P(>|Chi|)'[2]

    if (is.na(pval.new)) {
      return(list(fragility.index = NA, point.diagnostics = "algorithm did not converge"))
    }

    if (pval.new >= pval) {
      pval <- pval.new
      indices <- indices.new
      fragility.index <- fragility.index + 1
      pvalues <- append(pvalues, pval)
    }

    index <- index + 1
    iter <- iter + 1
  }
  if (iter >= nrow(data) - 3 * length(terms(formula))) {
    return(list(fragility.index = NA, point.diagnostics = "algorithm did not converge"))
  }

  if (not.significant) {
    resulting.pval <- anova(model, nullmodel)$'P(>|Chi|)'[2]
    point.diagnostics <- paste("No points removed. Covariate already not significant at confidence level", conf.level)
  } else{
    resulting.pval <- pvalues[-1]
    point.diagnostics <- data[indices, ]
  }
  point.diagnostics <- cbind(point.diagnostics, resulting.pval)

  return(list(fragility.index = fragility.index, point.diagnostics = point.diagnostics))
}


