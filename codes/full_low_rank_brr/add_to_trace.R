# copyright by the authors
add.to.trace <- function(model, data, trace=NULL, iter=NA, brr.vars.to.record=c('Psi','Gamma','Gamma.local.shrinkage','star.deltas','a3a4','brr.rank'), fa.vars.to.record=c('variances','local.shrinkage','rank','Eta','a1a2','deltas','Lambda','A'), print.interesting = FALSE) {
  
  # this function adds the current samples from the posteriors of model
  # paramters to the trace list
  
  if (is.null(trace)) {
  	# Initialize trace
  	trace <- list()
  
		var.names.to.add <- brr.vars.to.record
		if (!is.null(model$fa)) {
			var.names.to.add <- c(var.names.to.add, fa.vars.to.record)
		}
		
		for (name in var.names.to.add) {
			trace[[name]] <- list()
		}
  
  	trace$A <- list()
  	trace$iter <- NULL
  	trace$tpve <- list() # Always included
  
  	
  } else {
    # Add to existing trace
  	
  	curr.length <- length(trace$iter)
  
  	for (name in brr.vars.to.record) {
  		
  			trace[[name]][[curr.length+1]] <- model[['brr']]$context[[name]]
  	}
  
  	if (!is.null(model$fa)) {
  		for (name in fa.vars.to.record) {
  			trace[[name]][[curr.length+1]] <- model[['fa']]$context[[name]]
  		}
  	}
  	
  	trace$A[[curr.length+1]] <- model$A
  	trace$iter[curr.length+1] <- iter
  
  	total.variation.explained <- compute.amount.total.variance.explained(genotypes=data$genotypes, Psi=model$brr$context$Psi, Gamma=model$brr$context$Gamma)
  
  	total.variation.in.data <- sum(apply(data$phenotypes, 2, var))
  
  	trace$tpve[[curr.length+1]] <- total.variation.explained / total.variation.in.data
  	
  
  }
  
  return(trace)
}
