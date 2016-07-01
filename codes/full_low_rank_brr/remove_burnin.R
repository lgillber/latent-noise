# copyright by the authors
remove.burnin <- function(mcmc.output, burnin, n.to.use = NULL) {
  
  # this function removes the burnin from a trace list
  # parameters: 
  #   mcmc.output: the original (full) mcmc.output from 
  #                gibbs.full.low.rank.brr
  #   burnin: the number of burnin samples
  
	if (is.null(n.to.use)) {
		for (name in names(mcmc.output$traces)) {
			mcmc.output$traces[[name]] <- mcmc.output$traces[[name]][-seq(1:burnin)]
		}
	} 


	return(mcmc.output)
}

