
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sgmm

<!-- badges: start -->
<!--
[![R-CMD-check](https://github.com/SGDinference-Lab/SGDinference/workflows/R-CMD-check/badge.svg)](https://github.com/SGDinference-Lab/SGDinference/actions)
[![codecov](https://codecov.io/gh/SGDinference-Lab/SGDinference/branch/master/graph/badge.svg?token=YTBY15IXEP)](https://app.codecov.io/gh/SGDinference-Lab/SGDinference)
-->
<!-- badges: end -->
<!-- Gentle introduction TO BE ADDED 
__SGDinference__ is an R package that provides estimation and inference methods for large-scale mean and quantile regression models via stochastic (sub-)gradient descent (S-subGD) algorithms. The inference procedure handles cross-sectional data sequentially: 
&#10;  (i) updating the parameter estimate with each incoming "new observation", 
  (ii) aggregating it as a Polyak-Ruppert average, and 
  (iii) computing an asymptotically pivotal statistic for inference through random scaling.
&#10;The methodology used in the SGDinference package is described in detail in the following papers:
&#10;- Lee, S., Liao, Y., Seo, M.H. and Shin, Y., 2022. Fast and robust online inference with stochastic gradient descent via random scaling. In _Proceedings of the AAAI Conference on Artificial Intelligence_ (Vol. 36, No. 7, pp. 7381-7389).
<https://doi.org/10.1609/aaai.v36i7.20701>.
&#10;- Lee, S., Liao, Y., Seo, M.H. and Shin, Y., 2023. Fast Inference for Quantile Regression with Tens of Millions of Observations.    arXiv:2209.14502 [econ.EM] <https://doi.org/10.48550/arXiv.2209.14502>.
-->

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools") # if you have not installed "devtools" package
devtools::install_github("SGDinference-Lab/sgmm")
```

We begin by calling the **sgmm** package.

``` r
library(sgmm)
```
