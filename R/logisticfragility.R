#' @title Logistic Fragility Function
#' @description Compute the fragility of a coefficient in a logistic regression for dichotomous outcomes, i.e. the number of removed observations it would take to make a significant-result non-significant. Uses the glm() function from the stats package.
#'
#' @param formula Model formula which will be evaluated by glm()
#' @param data Dataframe with values for model forma, passed to glm()
#' @param covariate Vector of covariates to find fragility index for. Default is all covariates in formula
#' @param conf.level Significance level
#' @param verbose Logical indicating if function will return verbose results or only fragility index
#'
#' @importFrom stats glm
#' @importFrom stats terms
#' @importFrom stats terms.formula
#' @importFrom stats anova
#' @importFrom stats update
#' @importFrom stats as.formula
#' @importFrom stats residuals
#' @importFrom stats anova
#'
#' @examples
#' # Import and format example data
#' mydata <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
#' mydata$rank <- factor(mydata$rank)
#'
#' logisticfragility(admit ~ gre + gpa + rank, data = mydata, covariate="gpa", verbose = TRUE)
#'
#' logisticfragility(admit ~ gre + gpa + rank, data = mydata)
#'
#' @return If verbose is FALSE, returns a list with fragility indices for selected covariates. If
#' verbose is TRUE, returns a list with p-values for each fragility index
#' at each iteration of the algorithm.
#'
#' @export logisticfragility



logisticfragility <- function(formula, data, covariate = "all.factors.default", conf.level = 0.95, verbose = FALSE){

  if ("all.factors.default" %in% covariate) {
    object <- terms.formula(formula)
    covariate.names <- attr(object, "term.labels")
  } else {
    covariate.names <- covariate
  }

  result.store <- vector("list", length(covariate.names))
  names(result.store) = covariate.names

  for (i in 1:length(covariate.names)) {
    if (verbose==TRUE) {
      result <- logisticfragilityinternal(formula, data, covariate.names[i], conf.level)
    } else{
      result <- logisticfragilityinternal(formula, data, covariate.names[i], conf.level)
      result <- result[1]
    }
    result.store[[paste(covariate.names[i])]] <- result
    }
  return(result.store)
  }

logisticfragilityinternal <- function(formula, data, covariate, conf.level) {

  alpha <- (1 - conf.level)

  model <- glm(formula, data, family = "binomial")
  model <- update(model, .~.-1)
  formula <- model$formula

  nullmodel <- update(model, as.formula(paste(".~.-", covariate)))

  delta.resid <- model$residuals - nullmodel$residuals
  index <- c(1:length(delta.resid))
  y = formula[[2]]
  ordering <- cbind(index, delta.resid, data[,paste(y)])
  ordering <- cbind(ordering, (ordering[ ,2] - ordering[ ,3]*2*ordering[ ,2]))
  ordering <- ordering[order(-ordering[ ,4]), ]

  pval <- anova(model, nullmodel, test = "LRT")$`Pr(>Chi)`[2]
  index <- 1
  iter <- 0

  pvalues <- pval
  not.significant <- FALSE
  if (pval > alpha) {
    not.significant <- TRUE
  }
  while (pval <= alpha & (iter < nrow(data))) {
    points <- ordering[1:index,1]

    modified.data <- data[-points, ]

    newmodel <- glm(formula, modified.data, family = "binomial")
    newnullmodel <- update(newmodel, as.formula(paste(".~.-", covariate)))

    pval <- anova(newmodel, newnullmodel, test = "LRT")$`Pr(>Chi)`[2]
    if (is.na(pval)) {
      stop("algorithm did not converge")
    }

    index <- index + 1
    iter <- iter + 1
    pvalues <- append(pvalues, pval)
  }
  if (iter >= nrow(data)) {
    stop(list("algorithm did not converge"))
  }
  resulting.pval <- pvalues[-1]
  if (not.significant) {
    resulting.pval <- anova(model, nullmodel, test = "LRT")$`Pr(>Chi)`[2]
  }
  fragility.index <- index - 1
  point.diagnostics <- data[1:(index-1), ]
  if (not.significant) {
    point.diagnostics <- paste("No points removed. Covariate already not significant at confidence level", conf.level)
  }
  point.diagnostics <- cbind(point.diagnostics, resulting.pval)

  return(list(fragility.index = fragility.index, point.diagnostics = point.diagnostics))
}
