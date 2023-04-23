# covmat
Covmat is a collection of techniques for estimating convariance matrices. Covariance matrices can be built using missing data. Stambaugh Estimation and FMMC methods can be used to construct such matrices. Covariance matrices can be built by denoising or shrinking the eigenvalues of a sample covariance matrix. Such techniques work by exploiting the tools in Random Matrix Theory to analyse the distribution of eigenvalues. Covariance matrices can also be built assuming that data has many underlying regimes. Each regime is allowed to follow a Dynamic Conditional Correlation model. Robust covariance matrices can be constructed by multivariate cleaning and smoothing of noisy data.

Installation
------------

To get started, you can install the package from github using `devtools`.

``` r
library(devtools)
install_github("https://github.com/arorar/PortfolioConstruction.git")
```
