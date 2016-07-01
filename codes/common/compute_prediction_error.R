# copyright by the authors
compute.prediction.error <- function(test.data, mcmc.output, burnin = 0.5) {
	# this function computes the prediction error for the low-rank BRRR model


	# remove burn-in
	n.iter <- length(mcmc.output$traces$Gamma)	
	tmp.output <- remove.burnin(mcmc.output, burnin = round(burnin * n.iter))

	
	
	# temporary matrix initialization
	preds <- test.data$phenotypes
	preds[, ] <- 0
	
	
	# compute predictions
	for (i in 1:length(tmp.output$traces$Psi)) {
		
		# prediction
		preds	<- preds + 1/n.iter * test.data$genotypes %*% tmp.output$traces$Psi[[i]] %*% tmp.output$traces$Gamma[[i]]
		
	}
	
	
		
	# subtract prediction mean, compute variance
	tmp <- (test.data$phenotypes - preds)
	test.data.ptve <- 1 - sum(apply(tmp,2,var)) / sum(apply(test.data$phenotypes, 2, var))
		
		
		
	# traditional MSE
	MSE <- mean((tmp)^2) 
	
	return(list(MSE=MSE, test.data.ptve = test.data.ptve))
	
}
