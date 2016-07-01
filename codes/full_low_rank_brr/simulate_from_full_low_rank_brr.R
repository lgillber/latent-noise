# copyright by the authors
simulate.from.full.low.rank.brr <- function(true.model, mafs=NULL, genotypes = NULL) {

  # This function simulates data (genotypes, phenotypes) from a
  # low-rank BRRR model
  
  # parameters:
  #   true.model: the low-rank BRRR model to be used for generating data
  #   genotypes (optional): a matrix to be used as the covariates
  #                         if not provided, will be generated according
  #                         to mafs to resemble SNP data
  #   mafs: when covariates (parameter genotypes) not provided,
  #         this parameter must contain the minor-allele frequencies 
  #         for the simulated SNP data
	

	# Variables related to the low-rank covariance
	variances <- true.model$fa$context$variances
	Eta <- true.model$fa$context$Eta
	Lambda <- true.model$fa$context$Lambda
	
	# Variables related to the reduced-rank regression
	Gamma <- true.model$brr$context$Gamma
	Psi <- true.model$brr$context$Psi
	
	n.patients <- nrow(Eta)
	n.pheno <- length(variances)
	
	if (is.null(genotypes)) {
		
		n.snps <- nrow(Psi)	
		# Simulate genotypes
		genotypes <- simulate.genotypes(n.patients=n.patients, n.snps=n.snps, mafs = mafs)	
	}
	
	crossprod.genotypes <- crossprod(genotypes)	
	

	# Simulate the phenotypes
	noise.covariance <- diag(variances)
	noise <- mvrnorm(n=n.patients, mu=rep(0,n.pheno), Sigma=noise.covariance)
	
	
	# Regression coefficient matrix for the confounders
	if (!is.null(true.model$A)) {
		A <- true.model$A	
		n.confounders <- nrow(A)
		confounders <- matrix(rnorm(n=n.patients*n.confounders, mean=0, sd=1), nrow=n.patients, ncol=n.confounders)
		
		# simulate phenotypes with confounders
		phenotypes <- genotypes %*% Psi %*% Gamma + confounders %*% A + Eta %*% t(Lambda) + noise
	} else {
		
		confounders <- NA
		
		# simulate phenotypes WITHOUT confounders
		phenotypes <- genotypes %*% Psi %*% Gamma + Eta %*% t(Lambda) + noise
	}
	
	 if (!is.null(true.model$brr$context$Omega)) {
	 	phenotypes <- phenotypes + true.model$brr$context$Omega %*% Gamma
	 } 

	# Return the simulated genotypes, confounders and phenotypes, 
	# and pre-computed cross product for genotypes.
	return(list(genotypes=genotypes, phenotypes=phenotypes, confounders=confounders, crossprod.genotypes=crossprod.genotypes))
}


