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
#' @importFrom stats as.formula
#' 
#' @examples
#' # Import example data
#' ad <- "https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv"
#' mydata <- read.csv(file = ad, header= TRUE, sep=";")
#'   
#' 
#' formula = quality ~ fixed.acidity + citric.acid + residual.sugar + free.sulfur.dioxide + 
#'   total.sulfur.dioxide + pH + sulphates + alcohol
#' linearfragility(formula, data = mydata, covariate = c("citric.acid", 
#'                 "total.sulfur.dioxide", "free.sulfur.dioxide"))
#' \donttest{
#' # citric acid nonsignificant at 197 points removed and 
#' # residual.sugar is not significant at 0 points removed
#' linearfragility(quality ~ citric.acid + residual.sugar, data = mydata, verbose = TRUE)
#' }
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

linearfragilityinternal <- function(formula, data, covariate, conf.level) {
  alpha <- 1 - conf.level
  model <- lm(formula, data)
  nullmodel <- update(model, as.formula(paste(".~.-", covariate)))

  delta.resid <- model$residuals - nullmodel$residuals
  index <- c(1:length(delta.resid))
  y <- formula[[2]]
  ordering <- cbind(index, delta.resid, data[ ,paste(y)])
  ordering <- cbind(ordering, ordering[,2])
  ordering <- ordering[order(-ordering[,4]), ]

  pval <- anova(model, nullmodel, test = "LRT")$`Pr(>Chi)`[2]
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

  while (pval <= alpha & (iter < nrow(data))) {

    indices.new <- c(indices, index)
    point <- ordering[indices.new, 1]
    modified.data <- data[-point, ]
    newmodel <- lm(formula, modified.data)
    newnullmodel <- update(newmodel, as.formula(paste(".~.-", covariate)))
    pval.new <- anova(newmodel, newnullmodel, test = "LRT")$`Pr(>Chi)`[2]

    if (is.na(pval.new)) {
      return(list(fragility.index = NA, point.diagnostics = "algorithm did not converge"))
    }
    if (pval.new > pval) {
      pval <- pval.new
      indices <- indices.new
      fragility.index <- fragility.index + 1
      pvalues <- append(pvalues, pval)
    }

    index <- index + 1
    iter <- iter + 1
  }

  if (iter >= nrow(data)) {
    return(list(fragility.index = NA, point.diagnostics = "algorithm did not converge"))
  }

  if (not.significant) {
    resulting.pval <- anova(model, nullmodel, test = "LRT")$`Pr(>Chi)`[2]
    point.diagnostics <- paste("No points removed. Covariate already not significant at confidence level", conf.level)
  } else{
    resulting.pval <- pvalues[-1]
    removed <-  ordering[indices,1]
    point.diagnostics <- data[removed, ]
  }
  point.diagnostics <- cbind(point.diagnostics, resulting.pval)

  return(list(fragility.index = fragility.index, point.diagnostics = point.diagnostics))

}

