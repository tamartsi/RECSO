secondary.outcome.effect.identity <-
function(Y, design.mat.main.model, b.0, D, design.mat.bias.model, weights, p.Dx.pop, p.Dx.cc, max.iter = 500, eps = 1e-6){
	
	iter <- 1
	converge <- F

	m.case  <- summary(lm(Y[which(D == 1)] ~ design.mat.bias.model[which(D == 1),] ))$coef
	m.cont  <- summary(lm(Y[which(D == 0)] ~ design.mat.bias.model[which(D == 0),]))$coef	
	gam.x <- as.matrix(design.mat.bias.model) %*% (m.case[,1] - m.cont[,1])

	
	means.y.xd.hopt <-  design.mat.main.model %*% b.0 - (D - p.Dx.pop)*gam.x
	var.yx.pi.hopt  <- mean(((Y - means.y.xd.hopt)*weights)^2)
	
	tot.model.mat <- cbind(design.mat.main.model, diag(as.numeric(p.Dx.pop - D)) %*% design.mat.bias.model)
	
	
	temp.m <- summary(lm(gam.x[which(D == 1),] ~ design.mat.bias.model[which(D == 1),] ))$coef

	b.old <- b.0
	inds.b <- 1:length(b.old)
	b.old <- c(b.old, temp.m[,1])
	inds.a <- setdiff(1:length(b.old), inds.b)

		
	while (!converge & iter < max.iter){
			
				means.y.xd <- design.mat.main.model %*% b.old[inds.b] - (D - p.Dx.pop)*(as.matrix(design.mat.bias.model) %*% b.old[inds.a])
				term.1 <- Y - means.y.xd
				term.2 <- term.1*weights
				term.3 <- term.2/var.yx.pi.hopt
				

				
				newt.1 <- t(1/length(Y)*t(term.3) %*% tot.model.mat)  
				newt.2 <- -1/length(Y)*t(tot.model.mat) %*% diag(as.numeric(weights/var.yx.pi.hopt)) %*% tot.model.mat
				
				b.new <- b.old - solve(newt.2) %*% newt.1
			
				if (max(abs(b.new - b.old)) <= eps) converge <- T else{
					b.old <- b.new
					iter <- iter + 1	
				}
				
	}
	## now estimate the variance: (accounting for the estimation of p(X)) 
		means.y.xd <- design.mat.main.model %*% b.new[inds.b] - (D - p.Dx.pop)*(as.matrix(design.mat.bias.model) %*% b.new[inds.a])
		term.1 <- Y - means.y.xd
		term.2 <- term.1*weights
		term.3 <- term.2/var.yx.pi.hopt
		

		mat.1 <-  - 1/length(Y)*t(tot.model.mat) %*% diag(as.numeric(weights/var.yx.pi.hopt)) %*% tot.model.mat
		mat.2 <- diag(as.numeric(term.3)) %*% tot.model.mat/length(Y)
	
		#### terms accounting for the estimationg of p(X) = p(D=1|X). Suppose the parameter is called delta. (alpha in the paper!!!)
	


		### with p.Dx.cc:

		## This is the vector 
		E.du.ddelta <- p.Dx.pop*(1-p.Dx.pop)*(as.matrix(design.mat.bias.model) %*% b.new[inds.a])*weights/var.yx.pi.hopt
		E.du.ddelta <- t(tot.model.mat) %*%  diag(as.numeric(E.du.ddelta)) %*% design.mat.main.model/length(Y)
		
		E.dv.ddelta <- - t(design.mat.main.model) %*%  diag(as.numeric(p.Dx.cc*(1-p.Dx.cc))) %*% design.mat.main.model/length(Y)
		E.dv.ddelta.inv <- solve(E.dv.ddelta)	
		v.delta <- diag(as.numeric(D - p.Dx.cc)) %*% design.mat.main.model/length(Y)	
		term.account.delta.2 <-  v.delta %*% E.dv.ddelta.inv %*%  t(E.du.ddelta)


	
		mat.2 <- mat.2 - term.account.delta.2
	

	
		cov.mat <- solve(mat.1) %*% t(mat.2) %*% mat.2 %*% t(solve(mat.1))

		if (!is.null(colnames(design.mat.bias.model)))
		names(b.new) <- colnames(cov.mat) <- rownames(cov.mat) <-  c(colnames(design.mat.main.model) , colnames(design.mat.bias.model)) else 
			names(b.new) <- colnames(cov.mat) <- rownames(cov.mat) <-  c(colnames(design.mat.main.model) , "")
	
	
	return(list(beta.hat = b.new[inds.b], cov.beta.hat = cov.mat[inds.b, inds.b], alpha.hat = b.new[inds.a], cov.alpha.hat = cov.mat[inds.a, inds.a]))	
	
	
}
