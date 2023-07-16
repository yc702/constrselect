
# constrselect

<!-- badges: start -->
[![R-CMD-check](https://github.com/yc702/constrselect/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yc702/constrselect/actions/workflows/R-CMD-check.yaml)
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



### Survival outcome

Suppose two treatment arms are examined for patients with triple negative breast cancer and residual disease. The primary outcome of the study is EFS and the stratification of this study is based on nodal status. Previous studies reported that the node positive group had lower 2-year EFS than node negative group. The prevalence of node positive is around 30\%. 

Suppose the two strata of the inferior treatment arm have 2-year EFS $\bm S_b(2) = (0.6, 0.7)$, and sample size is determined based on an improvement of $\bm \theta^*= (0.15, 0.15)$ for the better treatment arm. Suppose patients enroll according to a Poisson process with accrual rate of 8 patients per year for each of the treatment arm stratum. We will continue follow up for additional 2 years after the last patient is enrolled for each stratum. Suppose the survival time follows exponential distribution and we are constraining and comparing survival probabilities at 2 years. Based on 5000 Monte Carlo simulations, we only need sample size 30 to achieve $\lambda = 0.8$ with $\rho = 0.5$ and $\bm \theta = (0.02,0.02)$. 

``` r
## basic example code

result <- pickwin_surv_fun(maxn=30,prop=c(0.3,0.7),event_rate_A=c(-log(0.85)/2,-log(0.75)/2),
                           trt_diff=c(log(0.75)/2-log(0.6)/2,log(0.85)/2-log(0.7)/2),
                           d=c(0.02,0.02), arrival_rate=8,FUP=2,x=2,S=5000,
                           study = "Constrained",cluster=2,
                           order_list=list(1,2),with_seed = 111)
                           
(sum(result$Corr)+(5000-sum(result$Corr)-sum(result$Error))/2)/5000



```
