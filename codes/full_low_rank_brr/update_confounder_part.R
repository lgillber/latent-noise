# copyright by the authors
update.confounder.part <- function(data, residuals, model) {
  
  # if there are known confounders, update the part of the model
  # that accounts for them

	t1 <- proc.time()
	n.confounders <- ncol(data$confounders)
	n.pheno <- ncol(data$phenotypes)

	variances <- model$fa$context$variances

	A <- rep(NA, n.confounders, n.pheno)

	pars <- fit.bayes.lm.ref.prior(X=data$confounders, y=residuals, noise.var=variances, crossprod.X=data$crossprod.X)

	for (i in 1:n.pheno) {
		# Regression for each phenotype can be
		# fitted separately, as independent error
		# terms are assumed (after removing the
		# effects of the latent factors).

		
		model$A[,i] <- mvr.norm.own(mu=pars$posterior.mean[[i]], Sigma=pars$posterior.cov[[i]])
	}
	t2 <- proc.time()
	cpu.times <- t2[3]-t1[3]
	names(cpu.times) <- 'A'
	return(list(model=model, cpu.times=cpu.times))
}

