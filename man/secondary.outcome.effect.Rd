\name{secondary.outcome.effect}
\alias{secondary.outcome.effect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The effect of exposure on secondary outcome
}
\description{
The main function of the package. 
}
\usage{
secondary.outcome.effect(Y, X, D, main.effect.vars, disease.vars, 
selection.bias.function.vars, p.D.pop, b.init.type = "ipw", 
link = "identity", p.Dx.cc = NULL, max.iter = 500, eps = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
Vector of the secondary outcome of interest. 
}
  \item{X}{
Matrix of all covariates and exposure variables needed (could have more variables than actually used). 
}
  \item{D}{
Vector of disease status
}
  \item{main.effect.vars}{
Exposure/covariates to estimate their effect on the secondary outcome, in the general populations. Must be column names of X. 
}
  \item{disease.vars}{
Exposure/covariates to estimate their effect on the disease, in the general/case-control study populations. Must be column names of X. 
}
  \item{selection.bias.function.vars}{
Exposure/covariates to estimate their effect on the selection bias, in the general populations. Must be column names of X. 
}
  \item{p.D.pop}{
Disease prevalence in the general population
}
  \item{b.init.type}{
Type of initial estimator for the main effect? Could be "ipw" (default), "pooled"  (based on all observations, ignoring case-control sampling) or "controls" (estimated from controls only data).
}
  \item{link}{
Type of link functions. Either "identity" or "log". 
}
  \item{p.Dx.cc}{
Probability of disease given covariates, in the case control study, for each subject. If not given, this is modeled in a logistic regression. 
}
  \item{max.iter}{
Maximum number of iterations of the Newton-Raphson algorithm. Default is 500.
}
  \item{eps}{
Criterion for convergence of the estimated vector. (Maximum absolute value of the difference between every two entries of the vector). 
}
}
\value{
\item{beta.hat}{Vector giving the estimated coefficients of the effect of main.effect.vars (and intercept) on Y.}
\item{cov.beta.hat}{Estimated covariance matrix of the beta.hat variables.}
\item{alpha.hat}{Vector giving the estimated coefficients of the effect of selection.bias.function.vars (and intercept) on the selection bias function.}
\item{cov.alpha.hat}{Estimated covariance matrix of the alpha.hat variables.}
}
\references{
Sofer, T., Tchetgen Tchetgen, EJ. ``Control function assisted IPW estimation with a
secondary outcome in case-control studies" (In preparation). 
}
\author{
Tamar Sofer
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
