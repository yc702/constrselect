
# constrselect

<!-- badges: start -->
[![R-CMD-check](https://github.com/yc702/constrselect/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yc702/constrselect/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/yc702/constrselect/branch/master/graph/badge.svg)](https://app.codecov.io/gh/yc702/constrselect?branch=master)
<!-- badges: end -->

The goal of constrselect is implement methods in paper **Randomized Phase II Selection Design with Order Constrained Strata**

## Installation

You can install the development version of constrselect from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("yc702/constrselect")
```

## Example

This is a basic example which shows you how to solve a common problem:

### Binary outcome

Suppose a clinical trial has two treatment arms: treatment A versus treatment B. The primary outcome of the study is response rate and the stratification is based on lymph node only metastasis versus metastasis of other sites. Previous studies reported that there are around 30\% of patients with lymph node only metastasis having better response rates. 

Suppose the lymph node only group will have higher response rate than the other group. We assume two strata of the inferior treatment arm have response rates $\bm \pi_b = (0.4, 0.5)$, $\bm \theta^*= (0.2, 0.2)$. We see that the sample size $N$ per arm derived from our method is around 25 to achieve $\lambda = \rho \times P_{amb}+P_{corr}= 0.8$ with $\rho = 0.5$ and $\bm \theta = (0.05,0.05)$. 
``` r
library(constrselect)
## basic example code
pickwin_bin_exact(n = 25, p1 = 0.4, strata_diff = 0.1, D=c(0.2,0.2),d=c(0.05,0.05),prop.strat=0.3,study="Constrained")


```

