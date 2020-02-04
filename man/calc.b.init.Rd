\name{calc.b.init}
\alias{calc.b.init}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Starting valued of the main effect parameter
}
\description{
Calculates a starting point for parameter estimates. This function is called by the main function "secondary.outcome.effect". 
}
\usage{
calc.b.init(Y, design.mat.main.model, weights, b.init.type = "ipw", link = "identity")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
Secondary outcome of interest. 
}
  \item{design.mat.main.model}{
Design matrix for the effect of exposure (and covariates) on the secondary outcome
}
  \item{weights}{
Inverse probability weights (probability of being in the sample given disease status)
}
  \item{b.init.type}{
type of initial regression parameter estimator. Could be "ipw" (default), "pooled"  (based on all observations, ignoring case-control sampling) or "controls" (estimated from controls only data). 
}
  \item{link}{
Type of link functions. Either "identity" or "log". 
}
}
\value{
Initial estimator of the effect of covariates on the secondary outcome. 
}
