# copyright by the authors
compute.ptve.corr.coef.with.Omega <- function (test.set.indexes = NULL, true.model, data, simu.data.ptve, tot.var, res.var.ptve, lambda) {
  
  # This function computes correction coefficients that can be used
  # to match the relative weights of different sources of variance
  # (covariance, structured noise, residual noise) to the specifications
  # (simu.data.ptve, res.var.ptve). Parameter lambda is the proportion
  # of independent structured noise of total structured noise.

	var.H.LambdaT <- sum(apply(true.model$fa$context$Eta[test.set.indexes,] %*% t(true.model$fa$context$Lambda),2,var))
	

  var.Sigma <- sum(true.model$fa$context$variances)
  
  var.X.Psi.Gamma <- sum(apply(data$genotypes[test.set.indexes,] %*% true.model$brr$context$Psi %*% true.model$brr$context$Gamma,2,var))


	var.Omega.Gamma <- sum(apply(true.model$brr$context$Omega[test.set.indexes,] %*% true.model$brr$context$Gamma, 2, var))	


	if (is.na(tot.var)) {
	  # if the total variance is not specified, keep the current 
	  # weights. Otherwise the total variance will be scaled to the
	  # given value later.
	 
		tot.var <- var.H.LambdaT + var.Sigma + var.X.Psi.Gamma + var.Omega.Gamma	
	}
  
	
	
	var.portions <- c(var.X.Psi.Gamma, var.Omega.Gamma, var.H.LambdaT, var.Sigma)
	names(var.portions) <- c('var.X.Psi.Gamma', 'var.Omega.Gamma', 'var.H.LambdaT', 'var.Sigma')
	
	
	correction.coefs <- list()
	correction.coefs$Gamma <- (simu.data.ptve * tot.var) / var.X.Psi.Gamma
	correction.coefs$Gamma <- sqrt(correction.coefs$Gamma)
	
	correction.coefs$Sigma <- (res.var.ptve * tot.var)/var.Sigma

	
	if (var.Omega.Gamma != 0) {
		correction.coefs$Omega <- tot.var * (1 - res.var.ptve - simu.data.ptve) * (1-lambda) / (var.Omega.Gamma*correction.coefs$Gamma^2)
		correction.coefs$Omega <- sqrt(correction.coefs$Omega)	
	} else {
		correction.coefs$Omega <- 0
	}
	
	if (var.H.LambdaT != 0) {
  	correction.coefs$Lambda <- tot.var * (1 - res.var.ptve - simu.data.ptve) * lambda / var.H.LambdaT
  	correction.coefs$Lambda <- sqrt(correction.coefs$Lambda)
	} else {
		correction.coefs$Lambda <- 0
	}
	
	
	correction.coefs$orig.var.comps <- var.portions
	
	return(correction.coefs)
}

