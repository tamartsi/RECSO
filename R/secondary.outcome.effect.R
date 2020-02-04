secondary.outcome.effect <-
function(Y, X, D, main.effect.vars, disease.vars, selection.bias.function.vars, p.D.pop, b.init.type = "ipw", link = "identity", p.Dx.cc = NULL, max.iter = 500, eps = 1e-6){
	# Y is the secondary outcome of interest
	# X is a matrix/data.frame with all covariates to use in the various models
	# main.effect.vars are column names of X, corresponding to covariates in the population mean model E(Y|X)
	# disease.vars are column names of X, corresponding to covariates in the disease model p(D=1|X)
	# selection.bias.function.vars are column names of X, corresponding to covariates in the model of the "selection bias function" gamma(X,D)  (identity link) or nu(X,D) (log link)
	# p.D.pop is the disease prevalence in the general population
	# link is the link function to use (identity or log)
	
	stopifnot(length(Y) == nrow(X), length(Y) == length(D), all(main.effect.vars %in% colnames(X)), all(disease.vars %in% colnames(X)), all(selection.bias.function.vars %in% colnames(X)), p.D.pop > 0, p.D.pop <1, is.element(link, c("identity", "log")), all(is.element(D, c(0,1))), is.element(b.init.type, c("ipw", "pooled", "controls")))
	

	
	p.D.cc <- sum(D)/length(D) ## probabiliy of having the disease in the CC sample
	
	if (!is.element("intercept", colnames(X))) {
		X <- cbind(1, X)
		colnames(X)[1] <- "intercept"	
	}
	
	## setting up the design matrices
	if (!is.element("intercept", disease.vars)) disease.vars <- c("intercept", disease.vars)
	if (!is.element("intercept", main.effect.vars)) main.effect.vars <- c("intercept", main.effect.vars)
	if (!is.element("intercept", selection.bias.function.vars)) selection.bias.function.vars <- c("intercept", selection.bias.function.vars)
	design.mat.D.model <- as.matrix(X[, disease.vars])
	design.mat.main.model <- as.matrix(X[, main.effect.vars])
	design.mat.bias.model <- as.matrix(X[, selection.bias.function.vars])
	
	
	##### First step: estimate the disease model using logistic regression unless given by the user
	if (is.null(p.Dx.cc))   p.Dx.cc <- fitted(glm(D ~ design.mat.D.model, family = "binomial")) ## probability of D given X in CC sample
	p.Dx.pop <- expit( logit(p.Dx.cc) + log(p.D.pop*(1-p.D.cc)/(p.D.cc*(1-p.D.pop))) ) ## probability of D given X in the population
	
	# ##### Second step: calculate the IPW for cases and controls
	# w.d1 <- p.D.pop/p.D.cc
	# w.d0 <- (1 - p.D.pop)/(1-p.D.cc)
	# weights <- w.d1*D + w.d0*(1 - D)
	
	p.d1 <- 1/(p.D.pop*1000)*p.D.cc
	p.d0 <- 1/((1-p.D.pop)*1000)*(1-p.D.cc)
	pi.D <-  p.d1*D + p.d0*(1 - D)
	weights <- 1/pi.D
	
	b.0 <- calc.b.init(Y,design.mat.main.model, weights = weights, b.init.type = b.init.type, link = link)
	
	
	if (link == "identity") res <- secondary.outcome.effect.identity(Y, design.mat.main.model, b.0, D, design.mat.bias.model, weights = weights, p.Dx.pop = p.Dx.pop, p.Dx.cc = p.Dx.cc, max.iter = max.iter, eps = eps) else ## log link
		res <- secondary.outcome.effect.log(Y = Y, design.mat.main.model = design.mat.main.model, D = D, b.0 = b.0, design.mat.bias.model = design.mat.bias.model, weights = weights, p.Dx.pop = p.Dx.pop, p.Dx.cc = p.Dx.cc, max.iter = max.iter, eps = eps)
	
	return(res)
	
}
