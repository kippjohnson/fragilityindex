# Repository for the R Package Fragility Index

Kipp Johnson

Implements and extends the fragility index calculation as described in Walsh M, Srinathan SK, McAuley DF, et al. _The statistical significance of randomized controlled trial results is frequently fragile: a case for a Fragility Index_. Journal of clinical epidemiology. 67(6):622-8. 2014.

## Introduction

As originally defined, the fragility index is the number of patients with a different outcome it would require to change a result from significant to non-significant. Consider for example the following example situation: In a clinical trial, there are two groups of patients. In group 1, 15/40 patients have an adverse event. In group 2, 5/40 patients have the adverse event. We can test this for statistical significance in the following way:

```
> mat1 <- matrix(c(15,5,25,35), nrow=2)
> mat1
     [,1] [,2]
[1,]   15   25
[2,]    5   35
> fisher.test(mat1)$p.value
[1] 0.01877238
```

However, what if we move a single patient from Group 1 to Group 2?

```
> mat2 <- matrix(c(14,6,26,34), nrow=2)
> mat2
     [,1] [,2]
[1,]   14   26
[2,]    6   34
> fisher.test(mat2)$p.value
[1] 0.06916506
```

This result is no longer statistically significant at the alpha=0.05 level! This is a "fragile" statistical result, despite the moderately low initial p value. Because it took only a single patient moving from Group 1 to Group 2 to make the result non-significant, we can say this clinical trial has a fragility index of 1. If it had taken two patients with swapped outcomes, it would woul have a fragility index of 2, and so on. 

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

For a dichotomous outcome, fragility index is the number of patients who if they were switched from one result to another would make a significant result at a given P-value non-significant. A smaller index means the observed result is more fragile to small variations.

### Reverse Fragility Index

~~~
revfragility.index(6,5,50,50, verbose=TRUE, print.mat=FALSE)
~~~

This package also contains a function to compute the "reverse fragility index," or the number of patients it would require who if their results were swapped would take a conclusion from non-significant to significant. This may be applied to analyze the sensitivity of non-inferiority trials to small differences in event counts. A smaller index means the observed result is more fragile to small variations.

### Logistic Beta Coefficient Regression Fragility

We present a new method to calculate logistic regression coefficient fragility, or how many events it would take to change a signficiant logistic regression coefficient to non-significant at the given confidence level. To do this, we replace responses (which should be binary events, i.e. 0 or 1) with the opposite event until the coefficient is nonsignificant. If the initial regression coefficient (beta) is positive, we change a 1 event to a 0. If the regression coefficient is negative, we change a 0 event to a 1. 

We then count the number of times this replacement must be and then obtain a fragility index. To account for variability in the response replacement, we then repeat this process a number of times and take the mean of all of the computed fragility indices to obtain a single fragility index.

Examining fragility of a single covariate:

~~~~
mydata <- read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")
mydata$rank <- factor(mydata$rank)
logisticfragility(admit ~ gre + gpa + rank, data = mydata, covariate="gre", niter=100)
~~~~

Or looking at all covariates in one step:

~~~
logisticfragility(admit ~ gre + gpa + rank, data = mydata, covariate="all", niter=100, progress.bar=TRUE)
~~~

Example output:
~~~
[1] "Doing (Intercept)..."
   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 15s
[1] "Doing gre..."
   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 03s
[1] "Doing gpa..."
   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 08s
[1] "Doing rank2..."
   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 05s
[1] "Doing rank3..."
   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 28s
[1] "Doing rank4..."
   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 28s
  coefficient fragility.index
1 (Intercept)           65.05
2         gre           11.25
3         gpa           29.65
4       rank2           19.02
5       rank3          120.62
6       rank4          123.35
~~~


