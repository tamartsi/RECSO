calc.b.init <-
function(Y, design.mat.main.model, weights, b.init.type = "ipw", link = "identity"){
	
	stopifnot(length(Y) == nrow(design.mat.main.model), length(Y) == length(weights), is.element(link, c("identity", "log")), is.element(b.init.type, c("ipw", "pooled", "controls")))
	
	if (link == "identity"){
		if (b.init.type == "ipw") b.init <- summary(lm(Y ~ design.mat.main.model, weights = weights))$coef[,1] else
			if (b.init.type == "pooled") b.init <- summary(lm(Y ~ design.mat.main.model))$coef[,1] else  ## controls only estimator
				b.init <- summary(lm(Y[which(D == 0)] ~ design.mat.main.model[which(D == 0),]))$coef[,1]
		
	} else{ ## link == log
		if (b.init.type == "ipw") b.init <- summary(glm(Y ~ design.mat.main.model, weights = weights, family = "poisson"))$coef[,1] else
			if (b.init.type == "pooled") b.init <- summary(glm(Y ~ design.mat.main.model, family = "poisson"))$coef[,1] else  ## controls only estimator
				b.init <- summary(glm(Y[which(D == 0)] ~ design.mat.main.model[which(D == 0),], family = "poisson"))$coef[,1]

	}
	
	return(b.init)
}
