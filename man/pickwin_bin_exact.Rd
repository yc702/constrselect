% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{pickwin_bin_exact}
\alias{pickwin_bin_exact}
\title{Binary outcome estimation for two treatments with two strata.}
\usage{
pickwin_bin_exact(
  n,
  p_inf,
  D = c(0.15, 0.15),
  d = c(0.05, 0.05),
  prop.strat = 0.3,
  study = NULL,
  order_list
)
}
\arguments{
\item{n}{Total sample size for all strata.}

\item{p_inf}{A vector of true probability for inferior treatment stratum arm.}

\item{D}{A vector of two treatment arms differences for each stratum for the superior treatment, Default: c(0.15, 0.15).}

\item{d}{A vector of ambiguous region for each stratum, Default: c(0.05, 0.05).}

\item{prop.strat}{The sample size proportion for the inferior strata, Default: 0.3.}

\item{study}{Could be either "Constrained" or "Origin" for the two types of study design, Default: NULL.}

\item{order_list}{A list of strata order allowing for partial ordering grouped in a vector within a list, default list(1,2).}
}
\value{
Return a vector of pcorr and pamb representing $P_{corr}$ and $P_{amb}$.
}
\description{
Estimation of P_corr and P_amb for selection design between
two treatments with two strata using exact binomial method.
}
\details{
DETAILS
}
\examples{
library(constrselect)
pickwin_bin_exact(n = 50, p_inf = c(0.25,0.3), D=c(0.15,0.15),d=c(0.05,0.05),prop.strat=0.4,study="Constrained",order_list = list(1,2))
}
