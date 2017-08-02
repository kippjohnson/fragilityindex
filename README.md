# Repository for the R Package Fragility Index

Kipp Johnson and Eli Rapoport

Implements and extends the fragility index calculation as described in Walsh M, Srinathan SK, McAuley DF, et al. _The statistical significance of randomized controlled trial results is frequently fragile: a case for a Fragility Index_. Journal of clinical epidemiology. 67(6):622-8. 2014.

## Introduction

As originally defined, the fragility index is the number of patients with a different outcome it would require to change a result from significant to non-significant. Consider for example the following example situation: In a clinical trial, there are two groups of patients. In group 1, 15/40 patients have an adverse event. In group 2, 5/40 patients have the adverse event. We can test this for statistical significance in the following way:

```
> mat1 <- matrix(c(15,6,25,34), nrow=2)
> mat2
     [,1] [,2]
[1,]   15   25
[2,]    6   34
> fisher.test(mat2)$p.value
[1] 0.04060921
```

However, what if a single additional patient in the second group had an event?

```
> mat2 <- matrix(c(15,7,25,33), nrow=2)
> mat2
     [,1] [,2]
[1,]   15   25
[2,]    7   33
> fisher.test(mat2)$p.value
[1] 0.07836101
```

This result is no longer statistically significant at the alpha=0.05 level! This is a "fragile" statistical result, despite the moderately low initial p value. Because it took only a single additional patient to make the result non-significant, we can say this clinical trial has a fragility index of 1. If it had taken two patients, this test would would have a fragility index of 2, and so on. 

This package contains functions to automatically calculate fragility indices in several situations, as explained below.

## Installation Instructions

### Installation from github (recommended)

We recommend installing from Github to ensure you have the latest version of the R package. The most up-to-date version can be installed from this github repository with the following commands:

```
install.packages("devtools") # If you do not have the devtools package
library(devtools)

install_github('kippjohnson/fragilityindex')
library(fragilityindex)
```

### Installation from CRAN

CRAN accepts only periodic submissions, and thus R packages cannot be frequently updated. Older versions of the package can also be installed from CRAN as follows:

```
install.packages("fragilityindex")
library(fragilityindex)
```

## Functions

### Fragility Index

~~~
fragility.index(15, 5, 40, 40)
~~~

For a dichotomous outcome, fragility index is the additional number of patients with an event it would take to make a significant result at a given P-value non-significant. A smaller index means the observed result is more fragile to small variations.

### Reverse Fragility Index

~~~
revfragility.index(6,5,50,50, verbose=TRUE, print.mat=FALSE)
~~~

This package also contains a function to compute the "reverse fragility index," or the number of patients it would require who if they did not experience an event would take a conclusion from non-significant to significant. This may be applied to analyze the sensitivity of non-inferiority trials to small differences in event counts. A smaller index means the observed result is more fragile to small variations.

### Logistic Beta Coefficient Regression Fragility

We present a new method to calculate logistic regression coefficient fragility, or how many events it would take to change a significant logistic regression coefficient to non-significant at the given confidence level. To do this, we remove responses (which should be binary events, i.e. 0 or 1) until the coefficient is nonsignificant. Responses which support the significance of the coefficient are removed in order of importance. 

We then count the number of events that must be removed to obtain a fragility index.

Examining fragility of a single covariate:

~~~~
mydata <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")

logisticfragility(admit ~ gre + gpa + rank, data = mydata, covariate="gre")
~~~~

Or looking at all covariates in one step:

~~~
logisticfragility(admit ~ gre + gpa + rank, data = mydata)
~~~

Example output:
~~~
$gre
$gre$fragility.index
[1] 2


$gpa
$gpa$fragility.index
[1] 4


$rank
$rank$fragility.index
[1] 23
~~~

### Fragility index for survival analysis and linear regression

We also present a new method to calculate fragility index for survival analysis and for linear regressions in a similar way to the calculation performed for the logistic regression fragility index calculations. We use the coxph() function from the survival package in R to compute survival tests from the G-rho family of tests and the lm() function from the stats package in R to compute linear regressions.



Example:

~~~
library(survival)
head(lung) # example survival data
lung$status <- lung$status - 1 # we require 0/1 outcomes; this variable originally is coded as 1/2

survivalfragility(Surv(time, status) ~ pat.karno + strata(inst), data=lung, covariate="pat.karno")
~~~

