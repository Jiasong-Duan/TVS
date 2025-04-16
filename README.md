
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TVS

<!-- badges: start -->
<!-- badges: end -->

**TVS** is an R package for *Testing-driven Variable Selection in Modal
Regression*. It implements efficient Bayesian variable selection methods
based on a permutation-based hypothesis testing.

## Installation

You can install the development version of TVS from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("Jiasong-Duan/TVS")
```

Or,

``` r
# install.packages("devtools")
devtools::install_github("Jiasong-Duan/TVS")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(TVS)
# Load example data containing:
# - Response: matrix with 100 subjects (rows) and 1 column
# - Predictors: matrix with 100 subjects (rows) and 8 predictors (columns)
# Note: Among the 8 predictors, only the first and third are truly important predictors
data(data_tvs)
head(data_tvs$Y)
#>           [,1]
#> [1,]  6.040409
#> [2,] -1.285562
#> [3,]  1.061651
#> [4,]  2.965520
#> [5,]  1.003568
#> [6,]  6.905395
head(data_tvs$X)
#>            [,1]        [,2]       [,3]        [,4]       [,5]        [,6]
#> [1,]  0.3563288 -0.74948617  0.8232656 -1.17007162  1.0684467 -0.02063172
#> [2,] -1.9805344  0.02754626 -1.2848831 -0.34758835  0.9163480 -0.90366751
#> [3,] -0.8664830 -0.35761430  0.4064376 -0.56418173  0.7135874  2.14675010
#> [4,]  0.6639846  0.78216497 -1.0981723 -0.75721919  1.2618043  0.52784216
#> [5,] -2.1843897 -1.00216420  0.2883779  0.70201472  0.3897322 -1.16153602
#> [6,]  2.1740915 -0.13433870  0.5846641 -0.04805797 -1.3447321  0.32667616
#>             [,7]       [,8]
#> [1,]  0.08314367 -0.5458808
#> [2,]  2.78494703  0.5365853
#> [3,]  0.59528232  0.4196231
#> [4,]  0.07766635 -0.5836272
#> [5,] -0.55869661  0.8474600
#> [6,] -0.80676768  0.2660220
#Parameter estimation
par_est <- TVS_EM(wrapper_beta, wrapper_beta0, wrapper_nu, wrapper_gamma, data_tvs)
par_est[1:7]
#> $beta
#>               [,1]
#> [1,]  1.571404e+00
#> [2,]  5.007751e-02
#> [3,]  9.326624e-01
#> [4,] -4.979870e-04
#> [5,] -1.344290e-06
#> [6,] -1.361048e-11
#> [7,] -1.298470e-11
#> [8,]  1.215912e-11
#> 
#> $beta0
#> [1] 2.482194
#> 
#> $nu
#> [1] 2.305821
#> 
#> $gamma
#> [1] 1.707573
#> 
#> $theta
#> [1] 0.1393816
#> 
#> $iter
#> [1] 29
#> 
#> $converged
#> [1] TRUE
#Variable selection with pre-screening. Requires about 30 seconds.
VS_withscreening <- TVS_multi_stage(wrapper_beta, wrapper_beta0, wrapper_nu, wrapper_gamma, data_tvs)
#> [Info] Group screening done.
#> [Info] Individual screening done.
#> [Info] Final step done.
# The output contains selected variables and the p-values for each predictor
# while a p-value of "-1" denotes a variable that was not selected.
print(VS_withscreening)
#> $selected_indices
#>      [,1]
#> [1,]    1
#> [2,]    3
#> 
#> $p_values
#>             [,1]
#> [1,]  0.00000000
#> [2,]  0.08666667
#> [3,]  0.00000000
#> [4,] -1.00000000
#> [5,] -1.00000000
#> [6,] -1.00000000
#> [7,] -1.00000000
#> [8,] -1.00000000
#Variable selection without pre-screening. Requires about 80 seconds.
VS_noscreening <- TVS(wrapper_beta, wrapper_beta0, wrapper_nu, wrapper_gamma, data_tvs)
#> Predictor 1
#> Predictor 2
#> Predictor 3
#> Predictor 4
#> Predictor 5
#> Predictor 6
#> Predictor 7
#> Predictor 8
print(VS_noscreening)
#>           [,1]
#> [1,] 0.0000000
#> [2,] 0.1500000
#> [3,] 0.0000000
#> [4,] 0.2500000
#> [5,] 0.3766667
#> [6,] 0.7566667
#> [7,] 0.7733333
#> [8,] 0.6433333
```
