# copyright by the authors
compute.ptve <- function(data, true.model, test.set.indexes = NULL) {

  # this function computes the proportion of variance explained by
  # the covariates. Note that when the test set indexes are given,
  # this is done in the test set to allow for exact oracle performance.
	if (is.null(test.set.indexes)) {
		
		total.variation.explained <- compute.amount.total.variance.explained(genotypes=data$genotypes, Psi=true.model$brr$context$Psi, Gamma=true.model$brr$context$Gamma)
		
		total.variation.in.data <- sum(apply(data$phenotypes, 2, var))
		
		ptve <- total.variation.explained / total.variation.in.data
		
	} else {
		total.variation.explained <- compute.amount.total.variance.explained(genotypes=data$genotypes[test.set.indexes,], Psi=true.model$brr$context$Psi, Gamma=true.model$brr$context$Gamma)
		
		total.variation.in.data <- sum(apply(data$phenotypes[test.set.indexes,], 2, var))
		
		ptve <- total.variation.explained / total.variation.in.data
		
	}
	
	return(ptve)
	
}

