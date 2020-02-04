\name{RECSO-package}
\alias{RECSO-package}
\alias{RECSO}
\docType{package}
\title{
Robust and Efficient analysis using Control function, of Secondary Outcomes. 
}
\description{
Performs robust and efficient analysis of the effect of exposure on a secondary outcome in a case-control study. 
}
\details{
\tabular{ll}{
Package: \tab RECSO\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-08-21\cr
License: \tab What license is it under?\cr
}
The main function in this package is secondary.outcome.effect. It calculates the effect in the population of exposure covariates
on a secondary outcome in a case-control study. It uses the disease status, prevalence of disease in the population. It either models
the probability of disease given covariates using logistic regression, or accepts estimated probabilities for each subjects. It also models the ``selection bias" 
function resulting from the case-control sampling. 
}
\author{
Tamar Sofer

Maintainer: Tamar Sofer <tamar.sofer@post.harvard.edu>
}
\references{
Sofer, T. and Tchetgen Tchetgen EJ. ``Control function assisted IPW estimation with asecondary outcome in case-control studies" (In preparation). 
}
\examples{

### Identity link example:



n.case <- 500 
n.cont <- 500
n <- n.case + n.cont
D <- c(rep(1, n.case), rep(0, n.cont))
X1 <- c(rnorm(n.case, 0,4), rnorm(n.cont, 2,4))
X2 <- rnorm(n, 2,2)
p.D.pop <- 0.12
p.D.cc.samp <- n.case/(n.case + n.cont)
pi.X <- fitted(glm(D ~ X1, family = "binomial"))
p.DX <- expit( logit(pi.X) + 
log(p.D.pop*(1-p.D.cc.samp)/(p.D.cc.samp*(1-p.D.pop))) )			

mean.Y.DX <- 50 + 3*X1 - (D - p.DX)*(3 + 2*X1 - 2*X1*X2)
Y <- mean.Y.DX + rnorm(n, 0,4)
X <- data.frame(X1, X2)
X$X1X2 <- X$X1*X$X2

## misspecified model for selection bias function:
result <- secondary.outcome.effect(Y, X, D,  main.effect.vars = c("X1", "X2"), 
disease.vars = c("X1", "X2"), selection.bias.function.vars = c("X2"), 
p.D.pop = 0.12, b.init.type = "ipw", link = "identity")

## correctly specified model for selection bias function:
result <- secondary.outcome.effect(Y, X, D,  main.effect.vars = c("X1", "X2"), 
disease.vars = c("X1", "X2"), selection.bias.function.vars = c("X1", "X2", "X1X2"), 
p.D.pop = 0.12, b.init.type = "ipw", link = "identity")
  

### log link example:
beta <- c(3, 0.7, 0.4, 0.4)
alpha <- c(0.5, 0.3, 0.3,0)
p.D.pop <- 0.12
n.case <- 500
n.cont <- 500
n <- n.case + n.cont
D <- c(rep(1, n.case), rep(0, n.cont))
X1 <- c(rnorm(n.case, 1, 0.2), rnorm(n.cont, 1.5, 0.2))
X2 <- rnorm(n,1, 0.2)
p.D.cc.samp <- n.case/(n.case + n.cont)
pi.X <- fitted(glm(D ~ X1, family = "binomial"))
p.DX <- expit( logit(pi.X) + 
log(p.D.pop*(1-p.D.cc.samp)/(p.D.cc.samp*(1-p.D.pop))) )

X <- cbind(1,X1, X2, X1*X2) 
colnames(X) <- c("intercept", "X1", "X2", "X1X2")
mean.Y.DX <- exp(X \%*\% beta +  D*(X \%*\% alpha) - log(exp(X \%*\% alpha)*p.DX + 1-p.DX) ) 
Y <- rpois(n, mean.Y.DX)  

## misspecified model for selection bias function:
result <- secondary.outcome.effect(Y, X, D,  main.effect.vars = c("X1", "X2", "X1X2"), 
disease.vars = c("X1", "X2"), selection.bias.function.vars = c("X2"), 
p.D.pop = p.D.pop, b.init.type = "ipw", link = "log")

## correctly specified model for selection bias function:
result <- secondary.outcome.effect(Y, X, D,  main.effect.vars = c("X1", "X2", "X1X2"), 
disease.vars = c("X1", "X2"), selection.bias.function.vars = c("X1", "X2", "X1X2"), 
p.D.pop = p.D.pop, b.init.type = "ipw", link = "log")



}
