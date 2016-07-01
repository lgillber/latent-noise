# copyright by the authors
infinite.brr.gibbs <- function(n.iter=500, vars.to.update=c('Psi','Gamma','Psi.local.shrinkage','Gamma.local.shrinkage','star.deltas','a3a4','brr.rank'), context, prior, genotypes, phenotypes, crossprod.genotypes=NULL) {
	#
	# A function for running the updates for the
	# the infinite reduced rank regression model.
	#
	# Inputs:
	#   n.iter: number of MCMC iterations
	#   vars.to.update: list of variables to update
	#   context: contains the current values of the variables
	#   prior: contains the hyperparameter and the current
	#   	values of the diagonal elements of the phenotype
	#       covariance matrix
	#   genotypes: centered 0,1,2 coded genotypes (however,
	#       this assumption is not used anywhere, so,
	#       function can be used with other kinds of
	#       regressors as well.
	#   phenotypes: scaled and centered phenotypes
	#   
	#   crossprod.genotypes: pre-computed crossprod(genotypes)
	#
	#	new.parameterization: instead of shrinking columns/rows of Psi/Gamma, shrink twice Gamma
	#
	#
	# Outputs (a list with elements:
	#   updated.context: context with updated variables
	#   trace: MCMC traces of the updated variables
	#	time: total time taken by the updates
	#
	
	
	
	
	
	if (is.null(crossprod.genotypes)) {
		crossprod.genotypes = crossprod(genotypes)
	}
	
	# Unlist variables for easier access:
	Psi <- context$Psi
	Gamma <- context$Gamma
	Omega <- context$Omega
	Psi.local.shrinkage <- context$Psi.local.shrinkage
	Gamma.local.shrinkage <- context$Gamma.local.shrinkage
	star.deltas <- context$star.deltas
	star.taus <- cumprod(star.deltas)
	a3 <- context$a3a4[1]
	a4 <- context$a3a4[2]
	brr.rank <- context$brr.rank


	
	# Unlist hyperparameters and variances for easier access:
	local.shrinkage.nu <- prior$local.shrinkage.nu
	a3.shape <- prior$a3.shape
	a3.rate <- prior$a3.rate
	a3.lower.bound <- prior$a3.lower.bound
	a4.shape <- prior$a4.shape
	a4.rate <- prior$a4.rate
	a4.lower.bound <- prior$a4.lower.bound
	brr.factor.relevance.cutoff <- prior$brr.factor.relevance.cutoff
	alpha0 <- prior$alpha0
	alpha1 <- prior$alpha1
	
	
	
	
	
	if (is.null(prior$step.size)) {
		step.size <- 10
	} else {
		step.size <- prior$step.size
	}
	


	
	variances <- prior$variances
	precisions <- 1/variances

	
	if (!is.null(context$latent.noise.var)) {
		latent.noise.var <- context$latent.noise.var
	} else {
		latent.noise.var <- 1
	}
	
	
	

	
	n.snps <- ncol(genotypes)
	n.pheno <- ncol(phenotypes)
	n.patients <- nrow(genotypes)
	
	traces <- list()
	for (name in names(context)) {
		traces[[name]] <- list()
	}
	cpu.times <- rep(0, length(names(context)))
	names(cpu.times) <- names(context)
	


	for (iter in 1:n.iter) {
		

		if (any(vars.to.update=='Psi')) {
	
			
			t1 <- proc.time()
			Psi.variances <- rep(1, length(Psi))
			
			
			# see Stegle, NIPS 2011: 
			# "Efficient inference in matrix-variate Gaussian models
			# with iid observation noise"
			

			# first mean vector of the distribution
			GGT.inv <- chol2inv(chol(latent.noise.var*t(Gamma) %*% Gamma + diag(variances)))

		
			eig.gamma <- eigen((Gamma %*% GGT.inv)%*%t(Gamma), symmetric=T)
			eig.vals.part <- as.numeric(eig.gamma$values %x% GLOB.eigen.crosspr.genotypes$values + Psi.variances)
			
			
			post.mean <- as.numeric(  
				GLOB.eigen.crosspr.genotypes$vectors %*% 	 
					((t(GLOB.eigen.crosspr.genotypes$vectors) %*% 
					  	(crossprod(genotypes,phenotypes) %*% 
					  	GGT.inv)%*%t(Gamma) %*% eig.gamma$vectors) / 
					 	matrix(eig.vals.part, nrow=nrow(Psi)) ) %*% 
					t(eig.gamma$vectors) 
			)
			
			# the the variance
			vec.part <- rnorm(length(Psi)) / sqrt(eig.vals.part)
			post.variance.terms <- as.numeric( ( GLOB.eigen.crosspr.genotypes$vectors %*% matrix(vec.part, nrow=nrow(Psi)) %*% t(eig.gamma$vectors)) )
				
			Psi.vec <- post.mean + post.variance.terms
			Psi <- matrix(Psi.vec, nrow=n.snps, ncol=brr.rank)
			
			traces$Psi[[iter]] <- Psi
			t2 <- proc.time()
			cpu.times['Psi'] <- cpu.times['Psi'] + (t2[3] - t1[3])
			
		}
		
		
		
		if (any(vars.to.update=='Gamma')) {
		
			t1 <- proc.time()
				
			  # Gamma shrunk twice
				Gamma.element.precisions <- Gamma.local.shrinkage * (star.taus)^2	
				
				
				if (!is.null(context$Omega)) {
				
				  # latent-noise case	
					design.matrix <- (genotypes %*% Psi) + Omega
					
				} else {
					design.matrix <- genotypes %*% Psi
				}
				
				crossprod.X = crossprod(design.matrix)
				
				# Prior variances of the elements:
				Sigma.gamma.p <- 1/Gamma.element.precisions
				fit <- fit.bayes.lm.diag(design.matrix, phenotypes, variances, Sigma.gamma.p, crossprod.X=crossprod.X)
				
				
				for (p in 1:n.pheno) {
					# Update the pth column of Gamma
					
					# could be made faster by just computing the cholesky factors
					# and using them directly instead of computing covariance
					# matrices
					
					Gamma[,p] <- mvr.norm.own(mu=fit$posterior.mean[[p]], Sigma=fit$posterior.cov[[p]])
					
				}	
				
			

			
			traces$Gamma[[iter]] <- Gamma
			t2 <- proc.time()
			cpu.times['Gamma'] <- cpu.times['Gamma'] + (t2[3] - t1[3])
			
		}
		

		
		## Sample star.deltas (and star.taus)
		if (any(vars.to.update=='star.deltas')) {
			
			
			# Sample star.deltas[1]
			t1 <- proc.time()
			
		
			# effect of Gamma
			gamma.terms <- compute.gamma.terms(brr.rank = brr.rank, Gamma = Gamma, Gamma.local.shrinkage = Gamma.local.shrinkage)
			
			# first update the first star delta parameter
			current.val <- star.deltas[1]
			
			
			star.deltas[1] <- 1
			tau.vector <- (cumprod(star.deltas))^2
			exp.sum.term <- sum(gamma.terms * tau.vector)
			
			

			# approximate distribution with a grid (as it can be computed quickly)
			proposal <- seq(from=1e-6, to=5*current.val,by=min(0.1, current.val/3))
			ln.prop.to.pr.delta1 <- (n.pheno * brr.rank + a3 -1) * log(proposal) - proposal - 0.5 * (proposal^2) * exp.sum.term
			

			
			
			
			# gibbs approximation code:
			star.deltas[1] <- normalize.ln.prop.to.distr(ln.prop.to.pr = ln.prop.to.pr.delta1, proposal=proposal)
			
			
			
			# if the rank is higher than 1, sample the remaining star delta
			# -parameters
			if (brr.rank>1) {
				# Sample star.deltas[2:]
				for (h in 2:brr.rank) {
					# Gaussian proposal; symmetric -> no need to evaluate jump probability
					current.val <- star.deltas[h]
					
					
					# approximate gibbs grid:
					proposal <- seq(from=1e-6, to=5*current.val,by=min(0.3, current.val/3))
					
					
					star.deltas[h] <- 1
					selector <- -seq(1,h-1,by=1)
					
				
					tau.vector <- ( cumprod(star.deltas) )^2
				
					exp.sum.term <- sum( (gamma.terms * tau.vector)[selector] )
					ln.prop.to.pr.delta.h <- (n.pheno * (brr.rank - h + 1) + a4 -1) * log(proposal) - proposal - 0.5 *(proposal^2) * exp.sum.term
					
					# use approximate gibbs:
					star.deltas[h] <- normalize.ln.prop.to.distr(ln.prop.to.pr = ln.prop.to.pr.delta.h, proposal=proposal)
				

				}
			}				

			star.taus <- cumprod(star.deltas)
			traces$star.deltas[[iter]] <- star.deltas
			t2 <- proc.time()
			cpu.times['star.deltas'] <- cpu.times['star.deltas'] + (t2[3] - t1[3])
			

		}
		
		if (any(vars.to.update=='Gamma.local.shrinkage')) {
			
			
			# Gamma.local.shrinkage[j,h] is the same as \phi_{jh}^{\Gamma} in the article.
			t1 <- proc.time()
			
			
			num <- table(1:ncol(Gamma))
			num <- matrix(rep(num, brr.rank), nrow=brr.rank, byrow=TRUE)
			shape.pars <- (local.shrinkage.nu + num) / 2
			rate.pars <- shape.pars
			rate.pars[,] <- NA
			
			
			Gamma.local.shrinkage.elements <- matrix(NA, brr.rank, ncol(Gamma))
			for (cl in 1:ncol(Gamma)) {
				# update local shrinkages for all target variables
				
				gSg <- diag(Gamma[, cl,drop=F] %*% t(Gamma[, cl,drop=F]))
				rate.pars[, cl] <- ( local.shrinkage.nu + (star.taus)^2 * gSg ) / 2	
			} 
			
			Gamma.local.shrinkage.elements <- matrix(rgamma(n=length(rate.pars), shape=shape.pars, rate=rate.pars), nrow=nrow(rate.pars), ncol=ncol(rate.pars))
			
			Gamma.local.shrinkage <- Gamma.local.shrinkage.elements
			
			
			traces$Gamma.local.shrinkage[[iter]] <- Gamma.local.shrinkage
			t2 <- proc.time()
			cpu.times['Gamma.local.shrinkage'] <- cpu.times['Gamma.local.shrinkage'] + (t2[3] - t1[3])
		}
		
		
		## Update a3 and a4
		if (any(vars.to.update=='a3a4')) {
			# Sample a3 and a4
			t1 <- proc.time()
			
			
			if (!is.na(a3.rate)) {
				
				# Move from (a3,a4) to (a3.star, a4.star) is propsed
				
				a3.proposal.std <- log(a3)/step.size
				a4.proposal.std <- log(a4)/step.size
				log.a3.star <- rnorm(n=1, mean=log(a3), sd=a3.proposal.std)
				log.a4.star <- rnorm(n=1, mean=log(a4), sd=a4.proposal.std)
				a3.star <- exp(log.a3.star)
				a4.star <- exp(log.a4.star)
				
				
				a3.star.proposal.std <- log(a3.star)/step.size
				a4.star.proposal.std <- log(a4.star)/step.size
				
				
				if ((a3.star>a3.lower.bound) & (a4.star>a4.lower.bound)) {
					# Otherwise the move is automatically rejected.
				
					log.proposal.prob <- dnorm(log(a3.star), mean=log(a3), sd=a3.proposal.std, log=TRUE) + dnorm(log(a4.star), mean=log(a4), sd=a4.proposal.std, log=TRUE)
					log.inverse.proposal.prob <- dnorm(log(a3), mean=log(a3.star), sd=a3.star.proposal.std, log=TRUE) + dnorm(log(a4), mean=log(a4.star), sd=a4.star.proposal.std, log=TRUE)
					
					log.prob.current.a3 <- (a3-1) * log(star.deltas[1]) + (a3.shape-1) * log(a3) - a3.rate*a3 - lgamma(a3)
					
					if (brr.rank>1) {
						log.prob.current.a4 <- (a4-1) * sum(log(star.deltas[-1])) + (a4.shape-1)*log(a4) - a4.rate*a4 - (brr.rank-1) * lgamma(a4)
						
					} else {
						# value for a4 comes from the prior
						log.prob.current.a4 <- dgamma(a4, shape=a4.shape, rate=a4.rate, log=TRUE)
					}
					log.prob.current <- log.prob.current.a3 + log.prob.current.a4
					
					log.prob.a3.star <- (a3.star-1)*log(star.deltas[1]) + (a3.shape-1)* log(a3.star) - a3.rate*a3.star - lgamma(a3.star)
					
					if (brr.rank>1) {
						log.prob.a4.star <- (a4.star-1) * sum(log(star.deltas[-1])) + (a4.shape-1)*log(a4.star) - a4.rate*a4.star - (brr.rank-1) * lgamma(a4.star)
						
					} else {
						# Value for a4.star comes from the prior
						log.prob.a4.star <- dgamma(a4.star, shape=a4.shape, rate=a4.rate, log=TRUE)
					}
					log.prob.proposed <- log.prob.a3.star + log.prob.a4.star
					
					log.acceptance.prob <- min(0, log.prob.proposed + log.inverse.proposal.prob - log.prob.current - log.proposal.prob)
					acceptance.prob <- exp(log.acceptance.prob)
					
					if (runif(n=1,min=0,max=1)<acceptance.prob) {
						# Proposal is accepted
						a3 <- a3.star
						a4 <- a4.star
					}
				}
			}
			traces$a3a4[[iter]] <- c(a3,a4)
			t2 <- proc.time()
			cpu.times['a3a4'] <- cpu.times['a3a4'] + (t2[3] - t1[3])
		}
		
		
		## Adapt rank
		if (any(vars.to.update=='brr.rank')) {
			# Perform adaptation of the rank similarly to the infinite FA model
			
			# Omega does not need to be adapted, it is sampled after brr.rank and
			
			t1 <- proc.time()
			
			if (runif(n=1, min=0, max=1) < exp(alpha0 + alpha1*iter)) {
				
				relevance.score <- rep(NA, brr.rank)
				
				
				relevance.score <- rep(NA, brr.rank)
				for (i in 1:brr.rank) {
					relevance.score[i] <- max(abs(Psi[,i,drop=FALSE] %*% Gamma[i,,drop=FALSE]))
				}

				col.relevant <- (relevance.score > brr.factor.relevance.cutoff)
				
				if (all(col.relevant)) {
					## Add another column from the prior
					values.to.add <- simulate.new.brr.factor(star.deltas=star.deltas, local.shrinkage.nu=local.shrinkage.nu, a4=a4, n.pheno=n.pheno, n.snps=n.snps)
					
					
					
					Psi <- cbind(Psi, values.to.add$new.Psi.col)
					
					Gamma.local.shrinkage <- rbind(Gamma.local.shrinkage, values.to.add$new.Gamma.local.shrinkage.row)
					Gamma <- rbind(Gamma, values.to.add$new.Gamma.row)
					
					star.deltas <- c(star.deltas, values.to.add$new.star.delta)
					star.taus <- c(star.taus, values.to.add$new.star.tau)
					brr.rank <- brr.rank+1
					
				} else {
					# Remove the non-relevant columns; however, leave at least one column.
					if (length(col.relevant) > 1) {
						cols.to.remove <- which(!col.relevant)
						if (length(cols.to.remove)==ncol(Psi)) {
							# Check which column is the most relevant
							# and retain that column.
							retain.this <- which.max(relevance.score)
							cols.to.remove <- setdiff(cols.to.remove, retain.this)
						}
					
						Psi <- Psi[,-cols.to.remove, drop=FALSE]
						#Psi.local.shrinkage <- Psi.local.shrinkage[, -cols.to.remove, drop=FALSE]
					
						Gamma <- Gamma[-cols.to.remove, , drop=FALSE]
						Gamma.local.shrinkage <- Gamma.local.shrinkage[-cols.to.remove, , drop=FALSE]
						
						star.deltas <- star.deltas[-cols.to.remove]
						star.taus <- cumprod(star.deltas)
						brr.rank <- brr.rank - length(cols.to.remove)

					} 
					
					
				
				}
			}


			traces$brr.rank[[iter]] <- brr.rank
			t2 <- proc.time()
			cpu.times['brr.rank'] <- cpu.times['brr.rank'] + (t2[3] - t1[3])
		}

	if (any(vars.to.update=='Omega')) {
		
		t1 <- proc.time()
		precisions <- 1/variances		
		
		# fit Omega to residual
		phenotypes.tmp <- phenotypes - genotypes %*% Psi %*% Gamma
		
	
		
		cov.matrix <- chol2inv( chol( 1/latent.noise.var*diag(brr.rank) + crossprod( sqrt(precisions) * t(Gamma) ) ) )
		
		
		mean.vectors <- cov.matrix %*% t(precisions * t(Gamma)) %*% t(phenotypes.tmp)
		mean.vectors <- t(mean.vectors)
		
		noise <- mvrnorm(n=n.patients, mu=rep(0,brr.rank), Sigma=cov.matrix)
		
		Omega <- mean.vectors + noise
		
		
		traces$Omega[[iter]] <- Omega
		t2 <- proc.time()
		cpu.times['Omega'] <- cpu.times['Omega'] + (t2[3] - t1[3])
	}
	
	
	if (!exists('Gamma.times')) {
		Gamma.times <- rep(0,6)
	}
	
	if (!exists('Psi.times')) {
		Psi.times <- rep(0,5)
	}
	
	
	
	}
	
	for (name in names(context)){
		# Update current values of the variable to the context
		if (name=='a3a4') {
			context$a3a4 <- c(a3,a4)
		} else {
			eval(parse(text=paste('context$', name, '<-', name, sep='')))
		}
	}

	
	to.return <- list(updated.context=context, traces=traces, cpu.times=cpu.times)

	
	return(to.return)
}

