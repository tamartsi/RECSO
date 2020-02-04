\name{secondary.outcome.effect.log}
\alias{secondary.outcome.effect.log}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The effect of exposure on secondary outcome, modeled via the log link
}
\description{
This function is called by the main function secondary.outcome.effect.
}
\usage{
secondary.outcome.effect.log(Y, design.mat.main.model, b.0, D, 
design.mat.bias.model, weights, p.Dx.pop, p.Dx.cc, max.iter = 500, eps = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
Vector of the secondary outcome of interest.
}
  \item{design.mat.main.model}{
Design matrix for the effect of exposure (and covariates) on the secondary outcome
}
  \item{b.0}{
Starting value for main effect estimates
}
  \item{D}{
Vector of disease status
}
  \item{design.mat.bias.model}{
Design matrix for the effect of covariates on the function modeling the bias resulting from case-control sampling
}
  \item{weights}{
Inverse probability weights (probability of being in the sample given disease status)
}
  \item{p.Dx.pop}{
Estimated probability of disease given covariates in the population (calculated by secondary.outcome.effect)
}
  \item{p.Dx.cc}{
Estimated probability of disease given covariates in the case-control study (calculated by secondary.outcome.effect)
}
  \item{max.iter}{
Maximum number of iterations of the Newton-Raphson algorithm. Default is 500.
}
  \item{eps}{
Criterion for convergence of the estimated vector. (Maximum absolute value of the difference between every two entries of the vector). 
}
}
