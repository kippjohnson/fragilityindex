# fragilityindex
Kipp Johnson

## Repository for the R Package Fragility Index

Implements the fragility index calculation as described in Walsh M, Srinathan SK, McAuley DF, et al. The statistical significance of randomized controlled trial results is frequently fragile: a case for a Fragility Index. Journal of clinical epidemiology. 67(6):622-8. 2014.

### Fragility Index

For a dichotomous outcome, fragility index is the number of patients who if they were switched from one result to another would make a significant result at a given P-value non-significant. A smaller index means the observed result is more fragile to small variations.

### Reverse Fragility Index

This package also contains a function to compute the "reverse fragility index," or the number of patients it would require who if their results were swapped would take a conclusion from non-significant to significant. This may be applied to analyze the sensitivity of non-inferiority trials to small differences in event counts. A smaller index means the observed result is more fragile to small variations.
