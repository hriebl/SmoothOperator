
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SmoothOperator

<!-- badges: start -->

<!-- badges: end -->

> Easy Feature Engineering with mgcv Penalized Spline Smooths

A convenient API for working with penalized spline smooths, based on the
mgcv package. The smooths can easily be integrated with user-defined
models and algorithms.

## Installation

You can install the development version of SmoothOperator from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hriebl/SmoothOperator")
```

## Example

Define and construct a simple cubic regression spline:

``` r
library(SmoothOperator)
set.seed(1337)

df <- data.frame(x = runif(100))
smt <- Smooth$new("x", df, bs = "cr", k = 10)
smt
#> <Smooth>
#>   Public:
#>     absorb_constraints: TRUE
#>     add_centering_constraint: function () 
#>     add_knots: function () 
#>     add_point_constraints: function (data) 
#>     bs: cr
#>     clone: function (deep = FALSE) 
#>     constraints: NULL
#>     construct: function (data = NULL) 
#>     data: data.frame
#>     identity_penalty: FALSE
#>     initialize: function (terms, data, bs = "tp", k = 10) 
#>     k: 10
#>     knots: NULL
#>     m: NA
#>     remove_all_constraints: function () 
#>     scale_penalty: TRUE
#>     terms: x
#>     xt: NULL
#>   Private:
#>     add_constraints: function (new) 
#>     initialize_constraints: function () 
#>     s: function ()

mat <- smt$construct()
str(mat)
#> List of 3
#>  $ design_matrix   : num [1:100, 1:9] -0.1101 -0.1111 0.4374 -0.0771 -0.0406 ...
#>  $ penalty_matrices:List of 1
#>   ..$ : num [1:9, 1:9] 0.54038 -0.45371 0.26777 0.00597 0.07403 ...
#>  $ ranks           : num 8
```

Add and inspect knots and constraints:

``` r
smt$add_knots()
smt$knots
#>             x
#> 1  0.02522857
#> 2  0.14127644
#> 3  0.24540405
#> 4  0.33913968
#> 5  0.47802501
#> 6  0.59722110
#> 7  0.73213332
#> 8  0.84627031
#> 9  0.94763002
#> 10 0.99971561

smt$add_centering_constraint()
smt$constraints
#>            [,1]      [,2]      [,3]       [,4]      [,5]      [,6]      [,7]
#> [1,] 0.04813571 0.1139502 0.1262563 0.09708947 0.1047601 0.1184304 0.1267345
#>            [,8]     [,9]      [,10]
#> [1,] 0.08647797 0.124193 0.05397238

smt$add_point_constraints(data.frame(x = 0.5))
smt$constraints
#>             [,1]         [,2]       [,3]        [,4]      [,5]      [,6]
#> [1,] 0.048135705  0.113950156 0.12625629  0.09708947 0.1047601 0.1184304
#> [2,] 0.001184691 -0.008088538 0.03766261 -0.10002457 0.9166500 0.1864093
#>             [,7]       [,8]         [,9]       [,10]
#> [1,]  0.12673449 0.08647797  0.124193020 0.053972380
#> [2,] -0.04534949 0.01541751 -0.005727267 0.001865725
```
