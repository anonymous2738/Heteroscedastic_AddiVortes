AddiVortes
===========

An implementation of the (Bayesian) Additive Voronoi Tessellation (AddiVortes) algorithm in R.

Copyright (C) 2024

Adam Stone & John Paul Gosling

Department of Mathematical Sciences, Durham University
 
Setup Instructions
------------------

To install the AddiVortes functions in R, you can run the following code: 

```r

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
library(devtools)

source_url("https://raw.githubusercontent.com/anonymous2738/Heteroscedastic_AddiVortes/main/Algorithm.R")

```
The following function can then be used in Rstudio.

```r
 AddiVortes_Algorithm(y, x, m = 200, max_iter = 1200, burn_in = 200,
                       nu = 6, q =0.85, k = 3, sd = 0.8, Omega = 3,
                       lambda_rate = 25, YTest, XTest, IntialSigma = "Linear", thinning = 1)
```
### Arguments

`y`- Dependent variable for training (in sample) data. A numerical vector with index of output value corresponding to the row of observation in x.

`x`- Explanatory variables for training (in sample) data. A matrix with (as usual) rows corresponding to observations and columns to variables.

`m`- Number of tessellations.

`max_iter`- Number of iterations of the MCMC backfitting algorithm.

`Burn_in`- Number of iterations discarded before sampling posterior. 0 < burn_in < max_iter.

`nu`- Degrees of freedom for error variance prior.

`q`- The quantile of the prior that the rough estimate of σ is placed at. The closer the quantile is to 1, the more aggressive the fit will be as you are putting more prior weight on error standard deviations (σ) less than the rough estimate.

`k`- Number of prior standard deviations E(Y∣x)=f(x) is away from +/-.5. The response (y) is internally scaled to range from -.5 to .5. 

`sd` - Standard deviation of the normal distribution for new coordinates in tessellations.

`Omega`- Probability of covariate being included as a dimension in Tessellation prior. 

`lambda_rate`- The rate of the number of centres in a tessellation for Poisson distribution in tessellation prior.

`YTest`- Dependent variable for test (out of sample) data. These should have the same structure as y.

`XTest`- Explanatory variables for test (out of sample) data. These should have the same structure as x.

`IntialSigma`- Either "Linear" or "Naive". When “Naive”, a rough estimate of σ corresponds to the sample standard deviation of the transformed training response values.
If “Linear”, the rough estimate of σ is based on the residual standard deviation from a least-squares linear regression of Y on the original X variables.

`thinning`- Retaining every posterior post burn in sample equal to thinning. Default retains every posterior sample post burn in.

`plot_qq` - If TRUE will plot a predictive qqplot and give the estatistic.

### Values

`In_sample_RMSE`- The RMSE of the in-sample estimations.

`Out_of_sample_RMSE`- The RMSE of the out-of-sample estimations.

`e_stat` - If plot_qq is TRUE, then will give the e-statistic.

Real Datasets
-----------------------------

To import the real-world benchmark datasets used in the paper in Rstudio one can run the following code:

```r
source_url("https://raw.githubusercontent.com/Adam-Stone2/AddiVortes/main/Datasets.R")

```

For each dataset, this imports the full datasets and the feature matrix (X_dataset) and output variable (Y_dataset) so that the Addivortes algorithm can be implemented. For example, to run the AddiVortes algorithm for the Boston dataset, use the following code:

```r
Boston #Full Boston dataset

n <- length(Y_Boston)
TrainSet <- sort(sample.int(n,5*n/6))
TestSet <- 1:n
TestSet <- TestSet[! TestSet %in% TrainSet]

 AddiVortes_Algorithm(Y_Boston[TrainSet],X_Boston[TrainSet,],
                      200,2000,200,6,0.85,3,0.8,3,25,
                      Y_Boston[TestSet],X_Boston[TestSet,],
                      IntialSigma = "Linear")

```
Reproducing Figures in the paper 
---------------------------

To reproduce the figures in the paper, source the following Github page by running:

```r

source_url("https://raw.githubusercontent.com/Adam-Stone2/AddiVortes/main/CodeForFigures.R")

```

**Warning:** This figure uses parallel processing using 10 cores at a time, producing this figure is only recommended if you have 12+ cores.

