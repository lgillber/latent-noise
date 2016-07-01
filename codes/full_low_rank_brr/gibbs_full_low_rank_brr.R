# copyright by the authors
gibbs.full.low.rank.brr <- function(model, data, n.iter, thin, fixed.brr.rank=NULL,brr.vars.to.record=c('Gamma.local.shrinkage', 'star.deltas', 'Psi', 'Gamma', 'brr.rank', 'a3a4', 'latent.noise.var', 'Omega'), fa.vars.to.record=c('a1a2','rank', 'Lambda', 'variances'), brr.vars.to.fix=c('Psi.local.shrinkage','a3a4'), ind.struct.noise = FALSE, fa.vars.to.update = c('Lambda','variances','Eta','local.shrinkage','deltas','a1a2','rank'), latent.noise = TRUE) {
	#
	# Function for simulating from the posterior distribution of 
	# a the infinite Bayesian low-rank reduced-rank regression model, 
	#
	# Inputs:
	#	n.iter: The number of iterations to run
	#
	#	model: This is model from which the MCMC is started.
	#          This contains fields:
	#          fa$context: Context of the covariance part.
	#          fa$prior: Priors for the covariance part.
	#          brr$context: Context for the reduced-rank coef matrix.
	#          brr$prior: Priors for the reduced-rank coef matrix.
	#          A: Confounder to phenotype coefficient matrix.
	#
	#	data: this list has elements genotypes, phenotypes,
	#          confounders, crossprod.genotypes
	#
	#	thin: thinning parameter for the MCMC
	#
	# fixed.brr.rank: If this is a numeric value, this value 
	#          will be used as the fixed brr.rank
	#
	#
	#
	#	brr.vars.to.fix: these brr parameters are not updated
	#	brr/fa.vars.to.record: these parameters are included in the trace
  #	fa.vars.to.update: these fa parameters are updated
  #
	# ind.struct.noise: include the independent structured noise 
  #                   model (H \Lambda)
  #
  # latent.noise: include the latent-noise model (\Omega)
	#
	# Outputs, a list containing:
	#	model: updated model
	#
	#	traces: MCMC traces of the variables.
	#
	


	
	
	if (!ind.struct.noise) {
	  
		# if independent-noise structured noise not used, only variances
	  # are computed in sparse fa code
		fa.vars.to.record <- c('variances')
		fa.vars.to.update <- c('variances')
		model$fa$context$Lambda[,] <- 0
	}
	
	

	if (!is.null(fixed.brr.rank)) {
		# If brr.rank has been fixed, check
		# that the correct value is in the
		# brr context.
		if (model$brr$context$brr.rank != fixed.brr.rank) {
			stop('Initial brr rank not equal to the fixed value')
		}
	}
	

	# Initialize variables for storing MCMC outputs:
	traces <- add.to.trace(model=model, brr.vars.to.record = brr.vars.to.record, fa.vars.to.record=fa.vars.to.record)
	
	# computation times
	cpu.times <- NULL
	t.res <- 0

	

	
	
	# use global variables to avoid passing things back and forth between functions
	GLOB.crosspr.genotypes <<- crossprod(data$genotypes)
	GLOB.eigen.crosspr.genotypes <<- eigen(GLOB.crosspr.genotypes)

	perm.info <- rep(1,14)

	
  # generate n.iter samples from the posterior
	for (iter in 1:n.iter) {
		if (iter == 1 || ( (iter %% 100) == 0) ) {
		  
			print(iter)
			print(paste('brr rank', model$brr$context$brr.rank))
		}
		
	  

			## Update confounder to phenotype coefficients. (Only if there are
			## any confounders in the model...)
			if (!(is.na(data$confounders)[1])) {
				print('updating confounders')
				tX <- proc.time()[3]
				
				
				residuals <- compute.residuals(type='confounder', model=model, data=data)
				tY <- proc.time()[3]
				t.res <- t.res + tY-tX
				conf.updated <- update.confounder.part(data=data, residuals=residuals, model=model)
				model <- conf.updated$model
			} else {
				conf.updated <- NULL
			}

			## Update independent-noise -part 

				tX <- proc.time()[3]
				
				
				residuals <- compute.residuals(type='fa', model=model, data=data)
				tY <- proc.time()[3]
				t.res <- t.res + tY-tX
			
				fa.res <- sparse.fa.gibbs(n.iter=1, context=model$fa$context, prior=model$fa$prior, Y=residuals, vars.to.update = fa.vars.to.update)
				
				model$fa$context <- fa.res$context	
				

			
			## Update the reduced-rank regression part of the model.
			tX <- proc.time()[3]
			

			residuals <- compute.residuals(type='brr', model=model, data=data)
			tY <- proc.time()[3]
			t.res <- t.res + tY-tX
			
			
			brr.updated <- update.brr.part(data=data, residuals=residuals, model=model, fixed.brr.rank=fixed.brr.rank, brr.vars.to.fix=brr.vars.to.fix, current.iter=iter, latent.noise = latent.noise)
			
		
			
			model <- brr.updated$model
			
			
			
			
			## Add current state to trace
			if (iter %% thin==0 || iter == 1) {
				
				
				traces <- add.to.trace(model=model, data=data, trace=traces, iter=iter, brr.vars.to.record=brr.vars.to.record, fa.vars.to.record=fa.vars.to.record)

				curr.length <- length(traces$tpve)
				
			}

			## Update times consumed by each different update
			cpu.times <- record.cpu.times(cpu.times=cpu.times, conf.res=conf.updated, fa.res=fa.res, brr.res=brr.updated)
			
		
	}


	return(list(model=model, traces=traces, cpu.times=cpu.times, t.res=t.res, perm.info = perm.info))
}
