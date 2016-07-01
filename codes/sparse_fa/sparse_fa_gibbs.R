# copyright by the authors
sparse.fa.gibbs <- function(n.iter, context, prior, Y,
vars.to.update=c('Lambda','variances','Eta','local.shrinkage','deltas','a1a2','rank'), thin = NULL, burnin=NULL, only.variances = FALSE) {
	#
	# Function for simulating from the posterior distribution of 
	# a sparse Bayesian factor analysis model, as defined in 
	# Bhattacharya and Dunson (2011).
	
	# Inputs:
	#	n.iter: The number of iterations to run
	#
	#	context: contains initial values for all variables that
	#            will be updated by the algorithm: "Lambda", 
	#            "variances", "Eta", "local.shrinkage", "deltas",
	#            "a1a2", "rank" (Initialization from the prior 
	#            can be done by function "initialize.fa.from.prior").
	#            Rank must always be >= 2
	#
	#	prior: contains hyperparameters of the model:
	#          "local.shrinkage.nu", "a.sigma", "b.sigma",
	#          "a1.shape", "a1.rate", "a1.lower.bound",
	#          "a2.shape", "a2.rate", "a2.lower.bound",
	#          "factor.relevance.cutoff", "alpha0", "alpha1"
	#
	#	vars.to.update: Only the variables specified here will be
	#                   updated by the algorithm (must represent
	#                   a subset of variables from the "context").
  #
  # only.variances: only update the residual noise variance parameters
  #
  # thin, burnin: mcmc parameters 
	#
	# Outputs, a list containing:
	#	context: updated context
	#
	#	trace: MCMC traces of the variables.

	# Unlist all values in context and prior for  easier use in the function.
	
	library(MASS)
	
	if (only.variances) {
		
		vars.to.update = c('variances')
		a.sigma <- prior$a.sigma
		b.sigma <- prior$b.sigma
		
		n.patients <- nrow(Y)
		n.pheno <-  ncol(Y)
		
	} else {
	
		Lambda <- context$Lambda
		rank <- ncol(Lambda)
		if (rank<2) {
			stop('Rank must be >= 2')
		}
		n.pheno <- nrow(Lambda)
		
		variances <- context$variances
		precisions <- 1/variances
		
		Eta <- context$Eta
		Lambda <- context$Lambda
		n.patients <- nrow(Eta)
		
		local.shrinkage <- context$local.shrinkage
		deltas <- context$deltas
		taus <- cumprod(deltas)
		a1 <- context$a1a2[1]  # Hyperparameters a1 and a2.
		a2 <- context$a1a2[2]
		
		local.shrinkage.nu <- prior$local.shrinkage.nu
		a.sigma <- prior$a.sigma
		b.sigma <- prior$b.sigma
		
		# Parameters of the gamma distributions for a1 and a2.
		a1.shape <- prior$a1.shape
		a1.rate <- prior$a1.rate
		a1.lower.bound <- prior$a1.lower.bound
		a2.shape <- prior$a2.shape
		a2.rate <- prior$a2.rate
		a2.lower.bound <- prior$a2.lower.bound
		
		# Parameters related to the adaptation of the rank
		alpha0 <- prior$alpha0
		alpha1 <- prior$alpha1
		factor.relevance.cutoff <- prior$factor.relevance.cutoff
		
		if (!is.null(prior$prior.var.Eta)) {
			
			prior.var.Eta <- prior$prior.var.Eta
		} else {
			prior.var.Eta <- 1
		}
	}
	
	

	# Data structure for the MCMC traces
	trace <- list()
	for (var.name in names(context)) {
		trace[[var.name]] <- list()
	}
	cpu.times <- rep(0, length(names(context)))
	names(cpu.times) <- names(context)

	if (!is.null(burnin) && !is.null(thin)) {
		samples.to.store <- seq(from=round(burnin*n.iter),to=n.iter, by=thin)
	}

	
	
	

	for (iter in 1:n.iter) {
		

		
		## Update Lambda
		if (is.element('Lambda',vars.to.update)) {
			t1 <- proc.time()
			
			Lambda <- matrix(0, nrow=n.pheno, ncol=rank)
			
			crossprod.Eta.Y <- crossprod(Eta,Y)
			crossprod.Eta <- crossprod(Eta)
			for (j in 1:n.pheno) {

				D_j_inv <- diag(local.shrinkage[j,]*taus)
				
				cov.matrix <- chol2inv(chol(D_j_inv + precisions[j]*crossprod.Eta))
				
				mean.vector <- precisions[j] * cov.matrix %*% crossprod.Eta.Y[,j]

				
				Lambda[j,] <- mvr.norm.own(mu=mean.vector, Sigma=cov.matrix)

			}
			
			trace$Lambda[[iter]] <- Lambda
			t2 <- proc.time()
			cpu.times['Lambda'] <- cpu.times['Lambda'] + (t2[3] - t1[3])
		}

		## Update Eta
		if (is.element('Eta', vars.to.update)) {
			t1 <- proc.time()

			
			cov.matrix <- chol2inv( chol( 1/prior.var.Eta * diag(rank) + crossprod( sqrt(precisions) * Lambda ) ) )

			
			mean.vectors <- cov.matrix %*% t(precisions * Lambda) %*% t(Y)

			
			mean.vectors <- t(mean.vectors)

			noise <- mvrnorm(n=n.patients, mu=rep(0,rank), Sigma=cov.matrix)

			Eta <- mean.vectors + noise
			

			trace$Eta[[iter]] <- Eta
			t2 <- proc.time()
			cpu.times['Eta'] <- cpu.times['Eta'] + (t2[3] - t1[3])
		}
		
		## Update variances (and precisions)
		if (is.element('variances', vars.to.update)) {
			t1 <- proc.time()
			
			shape <- a.sigma + n.patients/2
			if (only.variances) {
				rates <- b.sigma + 0.5 * apply((Y)^2,2,sum)
			} else {
				rates <- b.sigma + 0.5 * apply((Y-tcrossprod(Eta, Lambda))^2,2,sum)	
			}
			
			
			precisions <- rgamma(n=n.pheno, shape=shape, rate=rates)
			variances <- 1/precisions
			
			
			trace$variances[[iter]] <- variances
			t2 <- proc.time()
			cpu.times['variances'] <- cpu.times['variances'] + (t2[3] - t1[3])
		}

		## Update local.shrinkage
		if (is.element('local.shrinkage',vars.to.update)) {
			# local.shrinkage[j,h] is the same as \phi_{jh} in the article.
			t1 <- proc.time()

			shape <- (local.shrinkage.nu + 1) / 2
			rate.pars <- (local.shrinkage.nu + t(taus * t(Lambda^2))) / 2

			local.shrinkage <- matrix(rgamma(n=prod(dim(rate.pars)), shape=shape, rate=rate.pars), nrow=nrow(rate.pars), ncol=ncol(rate.pars))

			trace$local.shrinkage[[iter]] <- local.shrinkage
			t2 <- proc.time()
			cpu.times['local.shrinkage'] <- cpu.times['local.shrinkage'] + (t2[3] - t1[3])
		}

		## Update deltas (and taus)
		if (is.element('deltas',vars.to.update)) {
			# Sample deltas[1]
			t1 <- proc.time()

			shape <- a1 + n.pheno * rank / 2
			
			phi.lambda.vector <- colSums(Lambda^2 * local.shrinkage)
			
			deltas[1] <- 1
			tau.vector <- cumprod(deltas)
			
			rate <- 1 + 1/2 * tau.vector %*% phi.lambda.vector

			deltas[1] <- rgamma(n=1, shape=shape, rate=rate)

			# Sample the rest of the deltas

			for (h in 2:rank) {
				
				shape <- a2 + n.pheno/2 * (rank-h+1)
				
				deltas[h] <- 1
				tau.vector <- cumprod(deltas)[-seq(1,h-1,by=1)]
				
				rate <- 1 + 1/2 * tau.vector %*% phi.lambda.vector[-seq(1,h-1,by=1)]

				deltas[h] <- rgamma(n=1, shape=shape, rate=rate)
			}

			taus <- cumprod(deltas)

			trace$deltas[[iter]] <- deltas
			t2 <- proc.time()
			cpu.times['deltas'] <- cpu.times['deltas'] + (t2[3] - t1[3])
		}


		## Update a1 and a2
		if (is.element('a1a2',vars.to.update)) {
			t1 <- proc.time()
			# Sample a1 and a2

			# Move from (a1,a2) to (a1.star, a2.star) is propsed
			a1.proposal.std <- 1
			a2.proposal.std <- 0.2
			a1.star <- rnorm(n=1, mean=a1, sd=a1.proposal.std)
			a2.star <- rnorm(n=1, mean=a2, sd=a2.proposal.std)
			
			if ((a1.star>a1.lower.bound) & (a2.star>a2.lower.bound)) {
				# Otherwise the move is automatically rejected.
				proposal.prob <- dnorm(a1.star, mean=a1, sd=a1.proposal.std) * dnorm(a2.star, mean=a2, sd=a2.proposal.std)
				inverse.proposal.prob <- dnorm(a1, mean=a1.star, sd=a1.proposal.std) * dnorm(a2.star, mean=a2, sd=a2.proposal.std)
				# This would actually cancel...

				prob.current.a1 <- deltas[1]^(a1-1) * a1^(a1.shape-1) * exp(-1*a1.rate*a1) / gamma(a1)
				prob.current.a2 <- prod(deltas[-1]^(a2-1)) * a2^(a2.shape-1) * exp(-1*a2.rate*a2) / (gamma(a2)^(rank-1))
				prob.current <- prob.current.a1 * prob.current.a2

				prob.a1.star <- deltas[1]^(a1.star-1) * a1.star^(a1.shape-1) * exp(-1*a1.rate*a1.star) / gamma(a1.star)
				prob.a2.star <- prod(deltas[-1]^(a2.star-1)) * a2.star^(a2.shape-1) * exp(-1*a2.rate*a2.star) / (gamma(a2.star)^(rank-1))
				prob.proposed <- prob.a1.star * prob.a2.star

				acceptance.prob <- min(1, prob.proposed * inverse.proposal.prob / prob.current / proposal.prob)

				if (runif(n=1,min=0,max=1)<acceptance.prob) {
					# Proposal is accepted
					a1 <- a1.star
					a2 <- a2.star
				}
			}
			trace$a1a2[[iter]] <- c(a1,a2)
			t2 <- proc.time()
			cpu.times['a1a2'] <- cpu.times['a1a2'] + (t2[3] - t1[3])
		}

		## Adapt rank
		if (is.element('rank',vars.to.update)) {
			# Perform adaptation of the rank as described in the paper
			t1 <- proc.time()

			if (runif(n=1, min=0, max=1) < exp(alpha0 + alpha1*iter)) {
				
				relevance.scores <- apply(Lambda,2,function(x){crossprod(x)})
				relevance.scores <- relevance.scores / relevance.scores[1]
				col.relevant <- relevance.scores > factor.relevance.cutoff

				#col.relevant <- apply(Lambda,2,function(x){max(abs(x))}) > factor.relevance.cutoff

				if (all(col.relevant)) {
					# Add another column from the prior
 					values.to.add <- simulate.new.factor(a.sigma=a.sigma, b.sigma=b.sigma, deltas=deltas, local.shrinkage.nu=local.shrinkage.nu, a2=a2, n.pheno=n.pheno, n.patients=n.patients, prior.var.Eta)

					Lambda <- cbind(Lambda, values.to.add$new.Lambda.col)
					Eta <- cbind(Eta, values.to.add$new.Eta.column)
					local.shrinkage <- cbind(local.shrinkage, values.to.add$new.local.shrinkage.column)
					deltas <- c(deltas, values.to.add$new.delta)
					taus <- c(taus, values.to.add$new.tau)
					rank <- rank+1
				} else {
					# Remove the non-relevant columns; however, leave at least two columns.

					cols.to.remove <- which(!col.relevant)
					if (length(cols.to.remove) >= ncol(Lambda)-1) {
						# We are about to remove all columns, 
						# or all but one. However, rank must be
						# at least two.
						
						# Check which two columns are most significant
						# and retain these.
						relevance.scores <- apply(Lambda,2,function(x){max(abs(x))})
						retain.these <- order(relevance.scores, decreasing=T)[c(1,2)]
						cols.to.remove <- setdiff(cols.to.remove, retain.these)
					}
					
					if (length(cols.to.remove)>0) {
						Lambda <- Lambda[,-cols.to.remove, drop=FALSE]
						Eta <- Eta[, -cols.to.remove, drop=FALSE]
						local.shrinkage <- local.shrinkage[, -cols.to.remove, drop=FALSE]
						deltas <- deltas[-cols.to.remove]
						taus <- cumprod(deltas)
						rank <- rank - length(cols.to.remove)
					}
				}
			}
			trace$rank[[iter]] <- rank
			t2 <- proc.time()
			cpu.times['rank'] <- cpu.times['rank'] + (t2[3] - t1[3])
		}



	}
	
	updated.context <- list()
	updated.context$variances <- variances
	
	if (!only.variances) {

		
		updated.context$Lambda <- Lambda
		
		updated.context$Eta <- Eta
		updated.context$rank <- rank
		updated.context$local.shrinkage <- local.shrinkage
		updated.context$deltas <- deltas
		updated.context$a1a2 <- c(a1,a2)
	}
	
	
	return(list(context=updated.context, trace=trace, cpu.times=cpu.times))
	

}
