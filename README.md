
# constrselect

<!-- badges: start -->
[![R-CMD-check](https://github.com/yc702/constrselect/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yc702/constrselect/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of *constrselect* is to implement methods in paper **Randomized Phase II Selection Design with Order Constrained Strata**

## Installation

You can install the development version of `constrselect` from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("yc702/constrselect")
```

## Introduction

In phase II trials, it is quite common to include heterogeneous patient subgroups with different prognoses in the same trial. Incorporating such patient heterogeneity or stratification into statistical calculation for sample size can improve efficiency and reduce sample sizes. This package focuses on the calculation of probability of correct selection ($\lambda$) for stratified randomized phase II selection designs comparing two treatment arms. We consider both binary and time-to-event outcomes in our development. Compared with methods that do not use order constraints, our method is shown to improve the probabilities of correct selection or reduce sample size . 


## Example

The followings are basic examples which show you how to use the functions to produce examples in the paper:

### Binary outcome

#### Two strata for randomization

Suppose a clinical trial has two treatment arms comparing treatment A versus treatment B. The primary outcome of the study is response rate and the stratification is based on lymph node only metastasis versus metastasis of other sites. Previous studies indicated that around 30\% of patients with lymph node only metastasis have better response rates. 

Suppose the lymph node only group will have a higher response rate than the other group. We assume two strata of the inferior treatment arm have response rates (0.4, 0.5) while the better treatment arm has (0.6, 0.7). Using exact binomial method,`pickwin_bin_exact()` function could calculate the $P_{corr}$ and $P_{amb}$. We see that the sample size $N$ per arm derived from our method is around 20 to achieve $\lambda = \rho \times P_{amb}+P_{corr}= 0.8$ with $\rho = 0.5$ and ambiguous regions (0.05,0.05). 
``` r
library(constrselect)
## basic example code
result = pickwin_bin_exact(n = 20, p_inf = c(0.4,0.5),
D=c(0.2,0.2),d=c(0.05,0.05),prop.strat=0.7,study="Constrained")

## lambda calculation with rho = 0.5
result[1]+0.5*result[2]

```

#### More than two strata for randomization

Suppose a clinical trial has three treatment arms comparing treatment A versus treatment B. The primary outcome of the study is response rate and the stratification is based on cancer stage 1,2,3. Previous studies indicated that the larger the cancer stage the worse the prognosis, with the sample proportion of 4:3:3. 

We assume three strata of the inferior treatment arm have response rates (0.3,0.4,0.5) while the better treatment arm have (0.45,0.55,0.65) for cancer stage 3,2,1. Using 5000 Monte Carlo simulations, `pickwin_bin_multiple()` function could calculate $P_{corr}$ and $P_{amb}$. For this ordering, we need to specify `order_list` to be `list(1,2,3)` which indicate the order constraints strata1 < strata2 < strata3. If we want to specify partial ordering, eg. strata1 < strata2 and strata1 < strata3, we could specify `order_list=list(1,c(2,3))`. We see that the sample size $N$ per arm derived from our method is around 58 to achieve $\lambda = \rho \times P_{amb}+P_{corr}= 0.8$ with $\rho = 0.5$ and ambiguous regions (0.02,0.02,0.02). 
``` r
library(constrselect)
## basic example code
result <- pickwin_bin_multiple(n = 58, p_inf = c(0.3,0.4,0.5), D=c(0.15,0.15,0.15),d=c(0.02,0.02,0.02),
prop.strat=c(0.3,0.3,0.4),study="Constrained",S = 5000,cluster=6,order_list=list(1,2,3))

## lambda calculation with rho = 0.5
(sum(result$Corr)+0.5*(5000-sum(result$Corr)-sum(result$Wrong)))/5000

```



### Survival outcome

Suppose two treatment arms are examined with the primary outcome of the study to be EFS and the stratification of this study to be based on nodal status. Previous studies showed that the node positive group had lower 2-year EFS than node negative group. The prevalence of node positive is around 30\%. 

Suppose the two strata of the inferior treatment arm have 2-year EFS (0.6, 0.7), and the sample size is determined based on an improvement of 0.15 for the better treatment arm. Suppose patients enroll according to a Poisson process with an accrual rate of 8 patients per year for each of the treatment arm stratum. We will continue follow-up for an additional 2 years after the last patient is enrolled for each stratum. Suppose the survival time follows exponential distribution and we are constraining and comparing survival probabilities at 2 years. Based on 8000 Monte Carlo simulations, we only need a sample size of 30 to achieve $\lambda =\rho \times P_{amb}+P_{corr}= 0.8$ with $\rho = 0.5$ and ambiguous region (0.02, 0.02). 

``` r
## basic example code

result <- pickwin_surv_fun(maxn=25,prop=c(0.3,0.7),
surv_inf=c(0.6,0.7),
surv_sup=c(0.75,0.85),
d=c(0.02,0.02), arrival_rate=8,FUP=2,x=2,S=8000,
study = "Constrained",cluster=2,order_list=list(1,2),with_seed = 111)
                           
## Pamb
pamb=8000-sum(result$Corr)-sum(result$Wrong)

## lambda calculation with rho = 0.5
(sum(result$Corr)+(pamb)/2)/8000



```
