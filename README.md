# fragilityindex
Kipp Johnson

## Repository for the R Package Fragility Index

Implements the fragility index calculation as described in Walsh M, Srinathan SK, McAuley DF, et al. The statistical significance of randomized controlled trial results is frequently fragile: a case for a Fragility Index. Journal of clinical epidemiology. 67(6):622-8. 2014.

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

### Logistic Regression Fragility

~~~~
mydata <- read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")
mydata$rank <- factor(mydata$rank)
logisticfragility(admit ~ gre + gpa + rank, data = mydata, covariate="gre", niter=100)
~~~~

We present a new method to calculate logistic regression coefficient fragility, or how many events will it take to change a signficiant logistic regression coefficient to non-significant at the given confidence level. To do this, we replace responses (which should be binary events, i.e. 0 or 1) with the opposite event until the event is nonsignificant. If the regression coefficient (beta) is positive, we change a 1 event to a 0. If the regression coefficient is negative, we change a 0 event to a 1. 

We then count the number of times this replacement must be done randomly and then obtain a fragility index. To account for variability, we then repeat this process a great number of times and take the mean of all of the computed fragility indices to obtain a single fragility index.

