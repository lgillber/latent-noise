# copyright by the authors
initialize.from.prior <- function(n.pheno=10, n.snps=15, n.patients=200, fa.rank=10, brr.rank=3, alpha0=-1, alpha1=-0.0005, fa.relevance.cutoff=0.01, brr.relevance.cutoff=0.01, local.shrinkage.nu=3, n.confounders=2, a3.shape=3, a3.rate=0.28, a3.lower.bound=2.1, a4.shape=4.1, a4.rate=0.31, a4.lower.bound=3.1, independent.noise = TRUE, a.sigma = 2.2, b.sigma=0.3, prior.var.Eta=1) {
	
  ## Initialization of the model
  #
  # parameters:
  #   n.pheno: the number of target variables (paper notation: K)
  #   n.snps=15: the number of covariates (paper notation: P)
  #   n.patients: the number of observations (paper notation: N)
  #   fa.rank: rank of the independent-noise part
  #   brr.rank: rank of the regression/latent-noise part
  #
  #   alpha0, alpha1: parameters that define the update schedule 
  #                   of the rank
  #   fa.relevance.cutoff, brr.relevance.cutoff: parameters that define
  #                                              the threshold used in
  #                                              shutting down model
  #                                              components
  #   local.shrinkage.nu: element-wise shrinkage parameter for Gamma and
  #                       Lambda parameters
  #   n.confounders: the number observed confounders in the model
  #   a3.shape, a3.rate, a3.lower.bound: parameters of the regression
  #                                      weight variance distribution 
  #                    NOTE: a3 is a1 in paper notation and vice versa!
  #   a4.shape, a4.rate=0.31, a4.lower.bound: parameters of the regression
  #                                           weight variance distribution   #                    NOTE: a4 is a2 in paper notation and vice versa!
  # 
  # independent.noise: include the terms for the independent structured 
  #              noise (H \Lambda)
  # a.sigma, b.sigma: parameters of the residual noise parameter 
  #                   distribution
  # prior.var.Eta: the variance of the latent factors for the 
  #                independent-noise model (H)
  # 
	# n.confounders: number of known confounders. NA = no confounders
  
	model <- list()

	
	if (independent.noise) {
		# Simulate the low-rank covariance:
			
		fa <- initialize.fa.from.prior(a1.shape=3, a1.rate=1, a1.lower.bound=2.1, a2.shape=3.1, a2.rate=1, a2.lower.bound=3.1, a.sigma=a.sigma, b.sigma=b.sigma, local.shrinkage.nu=local.shrinkage.nu, factor.relevance.cutoff=fa.relevance.cutoff, alpha0=alpha0, alpha1=alpha1, rank=fa.rank, n.patients=n.patients, n.pheno=n.pheno, SHRINKAGE = TRUE, prior.var.Eta=prior.var.Eta)	
		
		
	} else {
		
		# only simulate the variances
		fa <- initialize.fa.from.prior(a1.shape=18, a1.rate=2, a1.lower.bound=2, a2.shape=18, a2.rate=2, a2.lower.bound=3, a.sigma=2.2, b.sigma=0.3, local.shrinkage.nu=local.shrinkage.nu, factor.relevance.cutoff=fa.relevance.cutoff, alpha0=alpha0, alpha1=alpha1, rank=fa.rank, n.patients=n.patients, n.pheno=n.pheno, only.variances = TRUE)
		
		
	}
	

	
	model$fa$context <-fa$context
	model$fa$prior <- fa$prior	
	
	# Simulate the reduced-rank regression coefficients:
	inf.brr <- initialize.infinite.brr(local.shrinkage.nu=local.shrinkage.nu, a3.shape=a3.shape, a3.rate=a3.rate, a3.lower.bound=a3.lower.bound, a4.shape=a4.shape, a4.rate=a4.rate, a4.lower.bound=a4.lower.bound, brr.factor.relevance.cutoff=brr.relevance.cutoff, alpha0=alpha0, alpha1=alpha1, a.sigma=NA, b.sigma=NA, brr.rank=brr.rank, n.snps=n.snps, n.pheno=n.pheno)
	

	model$brr$context <- inf.brr$context
	model$brr$prior <- inf.brr$prior
	# NOTE: a.sigma is NA, therefore, variances are not simulated
	# (they are already in model$fa$context)


	# Simulate the regression coefficients for the confounders:
	# Just some random initialization.
	# Formally, A has improper prior.

	if (!is.na(n.confounders)) {
		model$A <- matrix(rnorm(n=n.confounders*n.pheno, mean=0, sd=1), nrow=n.confounders, ncol=n.pheno)  	
	}
	
	
	return(model)
}

