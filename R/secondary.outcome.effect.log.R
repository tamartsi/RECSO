secondary.outcome.effect.log <-
function(Y, design.mat.main.model, b.0, D, design.mat.bias.model, weights, p.Dx.pop, p.Dx.cc, max.iter = 500, eps = 1e-6){
	
	
	 m.case  <- summary(glm(Y[which(D == 1)] ~ design.mat.bias.model[which(D == 1),]  , family = "poisson"))$coef
	 m.cont  <- summary(glm(Y[which(D == 0)] ~ design.mat.bias.model[which(D == 0),] , family = "poisson"))$coef	

	 exp.nu.x.d.0 <- rep(1, length(Y))
	 exp.nu.x.d.1 <- exp(as.matrix(design.mat.bias.model) %*% m.case[,1])/exp(as.matrix(design.mat.bias.model) %*% m.cont[,1])
	 exp.nu.x.d    <- exp.nu.x.d.0*(1-D) + exp.nu.x.d.1*D
	 exp.bar.nu.x <- (1- p.Dx.pop) + exp.nu.x.d.1*p.Dx.pop
	 g <- exp.nu.x.d/exp.bar.nu.x

	
	
	
	### inds.b are indices of beta (the parameter of interest). 
	### inds.a are the indices of the parameter of nu(X,D).
	
	
	converge <- F

	b.old <- b.0	

	means.y.x.hopt <- exp(design.mat.main.model %*% b.old)
	means.y.xd.hopt <- exp(design.mat.main.model %*% b.old)*(D*exp.nu.x.d.1 + (1-D)*exp.nu.x.d.0)/exp.bar.nu.x
	var.yx.pi.hopt  <- means.y.x.hopt/exp.bar.nu.x*(weights*exp.nu.x.d.1 + weights/( weights - 1)*exp.nu.x.d.0)
	
	tot.model.mat <- cbind(design.mat.main.model, diag(as.numeric(D - p.Dx.pop/exp.bar.nu.x*exp.nu.x.d.1)) %*% design.mat.bias.model)
	
	
	iter <- 1	
	
	temp.m <- summary(lm(log( exp.nu.x.d[which( D == 1)]) ~ design.mat.bias.model[which( D == 1),] ))$coef

	inds.b <- 1:length(b.old)
	b.old <- c(b.old, temp.m[,1])
	inds.a <- setdiff(1:length(b.old), inds.b)

	b.init <- b.old
		
		while (!converge & iter < max.iter){

			
				exp.nu.xd <- exp(as.matrix(design.mat.bias.model) %*% b.old[inds.a]* D)
				exp.bar.nu.x <- exp(as.matrix(design.mat.bias.model) %*% b.old[inds.a])*p.Dx.pop + (1 - p.Dx.pop)
				
				means.y.xd <- exp(design.mat.main.model %*% b.old[inds.b])*exp.nu.xd/exp.bar.nu.x
				
				term.1 <- ( Y - means.y.xd)*weights
				term.2 <- means.y.xd.hopt/var.yx.pi.hopt ## this is h(x), up to the derivative by beta and "delta" (parameter of nu(X,D)).
				term.3 <- term.1*term.2
				
				newt.1 <- t(1/length(Y)*t(term.3) %*% tot.model.mat)
				newt.2 <- 1/length(Y)*t(tot.model.mat) %*% diag(as.numeric(- term.2*weights*means.y.xd)) %*% tot.model.mat  
				
				b.new <- b.old - solve(newt.2) %*% newt.1
			
				if (max(abs(b.new - b.old)) <= eps) converge <- T else{
					b.old <- b.new
					iter <- iter + 1	
				}
				
			}

	
	## now estimate the variance: (accounting for the estimation of p(X)) 
		exp.nu.xd <- exp(as.matrix(design.mat.bias.model) %*% b.new[inds.a]* D)
		exp.bar.nu.x <- exp(as.matrix(design.mat.bias.model) %*% b.new[inds.a])*p.Dx.pop + (1 - p.Dx.pop)
		means.y.xd <- exp(design.mat.main.model %*% b.new[inds.b])*exp.nu.xd/exp.bar.nu.x 
				
		term.1 <- ( Y - means.y.xd)*weights
		term.2 <- means.y.xd.hopt/var.yx.pi.hopt ## this is h(x), up to the derivative of mu.
		term.3 <- term.1*term.2

		mat.1 <-  1/length(Y)*t(tot.model.mat) %*% diag(as.numeric(- term.2*weights*means.y.xd)) %*% tot.model.mat
		mat.2 <- diag(as.numeric(term.3)) %*% tot.model.mat/length(Y)
	
		#### terms accounting for the estimationg of p(X) = p(D=1|X, S=1). Suppose the parameter is called delta. (alpha in the paper!!!)


		### with est.p.D.s:

		## This is the vector 
		E.du.ddelta <- term.2*weights/exp.bar.nu.x*means.y.xd*p.Dx.pop*(1-p.Dx.pop)*(exp(as.matrix(design.mat.bias.model) %*% b.new[inds.a]) -1)
		E.du.ddelta <- t(tot.model.mat) %*%  diag(as.numeric(E.du.ddelta)) %*% design.mat.main.model/length(Y)
		
		E.dv.ddelta <- - t(design.mat.main.model) %*%  diag(as.numeric(p.Dx.cc*(1-p.Dx.cc))) %*% design.mat.main.model/length(Y)
		E.dv.ddelta.inv <- solve(E.dv.ddelta)	
		v.delta <- diag(as.numeric( D - p.Dx.cc)) %*% design.mat.main.model/length(Y)	
		term.account.delta.2 <-  v.delta %*% E.dv.ddelta.inv %*%  t(E.du.ddelta)


	
		mat.2 <- mat.2 - term.account.delta.2
	

	
		cov.mat <- solve(mat.1) %*% t(mat.2) %*% mat.2 %*% t(solve(mat.1))

	if (!is.null(colnames(design.mat.bias.model)))
		names(b.new) <- colnames(cov.mat) <- rownames(cov.mat) <-  c(colnames(design.mat.main.model) , colnames(design.mat.bias.model)) else 
			names(b.new) <- colnames(cov.mat) <- rownames(cov.mat) <-  c(colnames(design.mat.main.model) , "")
	
	return(list(beta.hat = b.new[inds.b], cov.beta.hat = cov.mat[inds.b, inds.b], alpha.hat = b.new[inds.a], cov.alpha.hat = cov.mat[inds.a, inds.a]))

	
}
