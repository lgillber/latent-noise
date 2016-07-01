# copyright by the authors
update.brr.part <- function(data, residuals, model, fixed.brr.rank, brr.vars.to.fix, current.iter, brr.vars.to.update = NULL, latent.noise = TRUE) {


	
	# BRR updates require the variance parameters
	# that are contained in FA context. These are given as prior
	# information for the BRR updates.
	prior <- model$brr$prior
	prior$variances <- model$fa$context$variances
	
	
	if (is.null(brr.vars.to.update)) {
	  
	  # the default parameters to update
	  vars.to.update <- c('Psi', 'Gamma','Gamma.local.shrinkage', 'star.deltas','a3a4')
	  
	  # Update rank if fixed rank not given
		if (is.null(fixed.brr.rank)) vars.to.update <- c('brr.rank', vars.to.update)
		
	  # latent noise parameters
		if (latent.noise) vars.to.update <- c(vars.to.update, 'Omega')
		
	}
		

  # exclude variables to be fixed
	vars.to.update <- setdiff(vars.to.update, brr.vars.to.fix)
	
	
	# Update once all variables
	brr.res <- infinite.brr.gibbs(n.iter=1, vars.to.update=vars.to.update, context=model$brr$context, prior=prior, genotypes=data$genotypes, phenotypes=residuals, crossprod.genotypes=data$crossprod.genotypes)
	
	cpu.times <- brr.res$cpu.times
	model$brr$context <- brr.res$updated.context

	

	return(list(model=model, cpu.times=cpu.times))
}