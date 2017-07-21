#' @title Linear Fragility Function
#' @description Compute the fragility of a coefficient in a linear regression, i.e. the number of removed observations it would take to make a significant-result non-significant. Uses the lm() function from the stats package.
#'
#' @param formula Model formula which will be evaluated by lm()
#' @param data Dataframe with values for model forma, passed to lm()
#' @param covariate Vector of covariates to find fragility index for. Default is all covariates in formula
#' @param conf.level Significance level
#' @param verbose Logical indicating if function will return verbose results or only fragility index
#'
#' @importFrom stats lm
#' @importFrom stats terms.formula
#' @importFrom stats complete.cases
#' @importFrom stats anova
#'
#' @examples
#'
#'
#' @return If verbose is FALSE, returns a list with fragility indices for selected covariates. If
#' verbose is TRUE, returns a list with p-values for each fragility index at each iteration
#' of the algorithm.
#'
#' @export linearfragility

linearfragility <- function(formula, data, covariate = "all.factors.default", conf.level = 0.95, verbose = FALSE) {

  if ("all.factors.default" %in% covariate) {
    object <- terms.formula(formula)
    terms <- attr(object, "term.labels")
    factors <- attr(attr(object, "factors"), "dimnames")[[1]]
    covariate.names <- intersect(factors, terms)
  } else {
    covariate.names <- covariate
    terms <- covariate
  }

  if (!identical(sort(terms), sort(covariate.names))) {
    stop("Error: Formula has predictors which are not covariates!")
  }

  result.store <- vector("list", length(covariate.names))
  names(result.store) <- covariate.names

  data <- data[complete.cases(data[ ,covariate.names]), ]

  for (i in 1:length(covariate.names)) {
    result <- linearfragilityinternal(formula, data, covariate.names[i], conf.level)
    if (verbose == FALSE) {
      result <- result[1]
    }
    result.store[[paste(covariate.names[i])]] <- result
  }
  return(result.store)
}

