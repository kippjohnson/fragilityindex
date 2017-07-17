#' @title Survival Fragility Function
#' @description Compute the fragility of a coefficient in a survival test, i.e. the number of removed observations it would take to make a significant-result non-significant. Uses the coxph() function from the survival package.
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
#' @importFrom stats terms
#' @importFrom stats terms.formula
#' @importFrom stats anova
#' @importFrom stats update
#' @importFrom stats as.formula
#' @importFrom stats residuals
#' @importFrom stats anova

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
#' @return If verbose is FALSE, returns a list with fragility indices for selected covariates. If
#' verbose is TRUE, returns a list with p-values for each fragility index
#' at each iteration of the algorithm.
#'
#' @export survivalfragility

survivalfragility <- function(formula, data, covariate = "all.factors.default", conf.level = 0.95, verbose = FALSE){

  if ("all.factors.default" %in% covariate) {
    object <- terms.formula(formula)
    covariate.names <- attr(object, "term.labels")
  } else {
    covariate.names <- covariate
  }

  result.store <- vector("list", length(covariate.names))
  names(result.store) = covariate.names

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
  na.data <- model$na.action
  indata <- data[-na.data, ]
  model <- coxph(formula, indata)

  formula <- model$formula

  nullmodel <- update(model, as.formula(paste(".~.-", covariate)))



  delta.resid <- residuals(model, type="deviance") - residuals(nullmodel, type="deviance") #which residuals to use?
  index <- c(1:length(delta.resid))
  surv.status = formula[2][[1]][[3]]
  ordering <- cbind(index, delta.resid, indata[ ,paste(surv.status)])
  ordering <- cbind(ordering, (ordering[ ,2] - ordering[ ,3] * 2 * ordering[ ,2]))
  ordering <- ordering[order(-ordering[ ,4]), ]

  pval <- anova(model, nullmodel)$'P(>|Chi|)'[2]
  index <- 1
  iter <- 0

  pvalues <- pval
  not.significant <- FALSE
  if (pval > alpha) {
    not.significant <- TRUE
  }
  while (pval <= alpha & (iter <= nrow(indata) - 3 * length(terms(formula)))) {
    points <- ordering[1:index,1]

    modified.data <- indata[-points, ]

    newmodel <- coxph(formula, modified.data)
    newnullmodel <- update(newmodel, as.formula(paste(".~.-", covariate)))

    pval <- anova(newmodel, newnullmodel)$'P(>|Chi|)'[2]
    if (is.na(pval)) {
      return(list(fragility.index = NA, point.diagnostics = "algorithm did not converge"))
    }

    index <- index + 1
    iter <- iter + 1
    pvalues <- append(pvalues, pval)
  }
  if (iter >= nrow(indata) - 3 * length(terms(formula))) {
    return(list(fragility.index = NA, point.diagnostics = "algorithm did not converge"))
  }
  resulting.pval <- pvalues[-1]
  if (not.significant) {
    resulting.pval <- anova(model, nullmodel)$'P(>|Chi|)'[2]
  }
  fragility.index <- index - 1
  point.diagnostics <- indata[1:(index-1), ]
  if (not.significant) {
    point.diagnostics <- paste("No points removed. Covariate already not significant at confidence level", conf.level)
  }
  point.diagnostics <- cbind(point.diagnostics, resulting.pval)

  return(list(fragility.index = fragility.index, point.diagnostics = point.diagnostics))
}


