# copyright by the authors
generate.simulation.study.model <- function(n.pheno=10, n.snps=15, n.patients=200, fa.rank=10, brr.rank=3, alpha0=-1, alpha1=-0.0005, fa.relevance.cutoff=0.01, brr.relevance.cutoff=0.01, local.shrinkage.nu=3, n.confounders=2, a3.shape=30, a3.rate=0.28, a3.lower.bound=2, a4.shape=4.1, a4.rate=0.31, a4.lower.bound=3) {
	
	# this function generates the model in the simulation study 
  # 
	# see function 'initialize.from.prior' for the explanations 
  # of the parameters
	
	model <- list()

	
	# Simulate the low-rank covariance:
	#output.clusters <- sort(unique(output.clustering))
	#n.clusters <- length(output.clusters)
		
			
			
	fa <- initialize.fa.from.prior(a.sigma=2.2, b.sigma=0.5, a1.shape=3, a1.rate=1, a1.lower.bound=3, a2.shape=3, a2.rate=1, a2.lower.bound=3, local.shrinkage.nu=3, factor.relevance.cutoff=fa.relevance.cutoff, alpha0=-1, alpha1=-0.005, rank=fa.rank, n.patients=n.patients, n.pheno=n.pheno, SHRINKAGE = FALSE)
	# note: no shrinkage here, as \Lambda will be orthogonalized using
	# the gramm-schmidt process to generate orthogonal Lambda and Gamma
	# parameters
			

	model$fa$context <-fa$context
	model$fa$prior <- fa$prior	
	
	
	
	# Simulate the reduced-rank regression coefficients:
	# NOTE: a.sigma is NA, therefore, variances are not simulated
	# (they are already in model$fa$context)
	inf.brr <- initialize.infinite.brr(local.shrinkage.nu=local.shrinkage.nu, a3.shape=a3.shape, a3.rate=a3.rate, a3.lower.bound=a3.lower.bound, a4.shape=a4.shape, a4.rate=a4.rate, a4.lower.bound=a4.lower.bound, brr.factor.relevance.cutoff=brr.relevance.cutoff, alpha0=alpha0, alpha1=alpha1, a.sigma=NA, b.sigma=NA, brr.rank=brr.rank, n.snps=n.snps, n.pheno=n.pheno)
	
	

	inf.brr$context$Omega <- matrix(rnorm(n.patients*brr.rank), nrow=n.patients, ncol=brr.rank)	
	
	
	model$brr$context <- inf.brr$context
	model$brr$prior <- inf.brr$prior
	
	# the parameters of the model generated here
	# will be scaled later to get the desired proportions of 
	# variance
	model$brr$prior$latent.noise.var <- 1
	


	
	
	
	# in addition, it is rotated according to rotation
	rotation.matrix <- matrix(rnorm(n.pheno^2), n.pheno, n.pheno)
	rotation.matrix <- qr.Q(qr(rotation.matrix))
	
	
	sds.orig <- apply(model$fa$context$Lambda,2,sd)
	model$fa$context$Lambda <- t(rotation.matrix[1:fa.rank,])
	model$brr$context$Gamma <- rotation.matrix[(fa.rank+1:brr.rank),]
	
	model$fa$context$Lambda <- t(t(model$fa$context$Lambda)/(apply(model$fa$context$Lambda,2,sd) / sds.orig))
	
	model$brr$context$Gamma <- model$brr$context$Gamma/(apply(model$brr$context$Gamma,1,sd) / sds.orig)


	
	
	
	# Simulate the regression coefficients for the confounders:
	if (!is.na(n.confounders)) {
	  
	  # Some random initialization.
	  # Formally, A has improper prior.
	  
		model$A <- matrix(rnorm(n=n.confounders*n.pheno, mean=0, sd=1), nrow=n.confounders, ncol=n.pheno)  	
	}
	
	#  minor allele frequencies for the simulated SNP data
	model$mafs <- runif(n=n.snps, min=0.05, 0.5)
	
	
	
	return(model)
}
