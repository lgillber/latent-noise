# copyright by the authors
simulate.model.and.data <- function(n.pheno=60, n.snps=30, n.confounders = NA, n.patients=1000, fa.rank=3, brr.rank=3, test.set.indexes, X.ptve, res.var.ptve, lambda, total.var = NA) {
  
  # This function simulates models and data. The parameters do the
  # following:
  #
  # n.pheno: the number of target variables (K)
  # n.snps: the number of covariates (P)
  # n.confounders: the the number of known confounders
  # n.patients: the number of observations (N)
  # fa.rank: the rank of the independent-noise model (H \Lambda)
  # brr.rank: the rank of the regression and latent-noise model
  #           (\Psi, \Gamma, \Omega)
  # test.set.indexes: the elements to be used as the test set 
  # X.ptve: the proportion of variance explained by the covariates
  # res.var.ptve: the proportion of variance explained by the   
  #               the unstructured residual noise
  #
  # the remaining variance is assumed to be structured, either independent
  # or latent-noise. 
  #
  # lambda: the relative proportion of independent structured noise,
  #         \in [0,1], the relative proportion of latent-noise is 1-lambda
  # total.var: the total variation in the target variable data,
  #            sum(trace(var(Y))), where Y is the phenotypes -matrix
  #
	
  
  
	library(MASS)
	
		
	
	max.cors <- 1
	PTVE.ACCEPT <- FALSE
	
	while (max.cors > 0.96 || (PTVE.ACCEPT == FALSE)) {
		
		true.model <- generate.simulation.study.model(n.pheno=n.pheno, n.snps=n.snps, n.patients=n.patients, n.confounders = n.confounders,fa.rank=fa.rank, brr.rank=brr.rank, a3.shape=3, a3.rate=0.28, a3.lower.bound=2, a4.shape=3.1, a4.rate=0.31, a4.lower.bound=3, local.shrinkage.nu = 10)
		
		
		
		
		data <- simulate.from.full.low.rank.brr(true.model, mafs = true.model$mafs)


		# avoid numerical problems by requiring data to be not too collinear
		cors <- cor(data$phenotypes)
		max.cors <- max(sort(cors, decreasing=T)[n.pheno+1:20])
		
		
		
		
		# compute correction coefficients to normalize the variances from 
		# the different sources (covariates, structured noise, 
		# residual noise )to the values given as X.ptve, res.var.ptve, ...
		correction.coefs <- compute.ptve.corr.coef.with.Omega(test.set.indexes, true.model, data, simu.data.ptve = X.ptve, res.var.ptve = res.var.ptve, tot.var = total.var, lambda = lambda)
		tmp.var.comps <- correction.coefs$orig.var.comps / sum(correction.coefs$orig.var.comps)
		ptve <-  tmp.var.comps['var.X.Psi.Gamma']
		
		# check that the data is acceptable
		if ( ptve > 0 && (max.cors < 0.96)) {
			
			
			if (abs(log10(ptve) - log10(X.ptve)) > 0.5 ) {
				print('ptve')
				print(ptve)
				
				PTVE.ACCEPT <- TRUE
				
			}
			
		} else {
			
			PTVE.ACCEPT <- FALSE
		}
		
		
		print('sum of total phenotype variance')
		print(sum(apply(data$phenotypes,2,var)))
			
		
		
		# next fine tune effects to get exact effect sizes
		if (PTVE.ACCEPT) {
			
			for (coef.iter in 1:20) {
				
				print('fine tuning')
				print(coef.iter)
				#print(paste('max cors', max.cors))
				
				
				correction.coefs <- compute.ptve.corr.coef.with.Omega(test.set.indexes, true.model, data, simu.data.ptve = X.ptve, res.var.ptve = res.var.ptve, tot.var = total.var, lambda = lambda)

				
				if (coef.iter == 1) {
					print('sum of original variance')
					orig.sum.var <- sum(correction.coefs$orig.var.comps)
					print(orig.sum.var)
					tot.var <- orig.sum.var
				} 
				
				

				
				true.model$brr$context$Gamma <- true.model$brr$context$Gamma * correction.coefs$Gamma
				true.model$fa$context$Lambda <- true.model$fa$context$Lambda * correction.coefs$Lambda
				true.model$fa$context$variances <- true.model$fa$context$variances * correction.coefs$Sigma
				true.model$brr$context$Omega <- true.model$brr$context$Omega * correction.coefs$Omega
				
				
				
				
				
				
				data <- simulate.from.full.low.rank.brr(true.model, genotypes = data$genotypes)
				
				
				
				correction.coefs <- compute.ptve.corr.coef.with.Omega(test.set.indexes, true.model, data, simu.data.ptve = X.ptve, res.var.ptve = res.var.ptve, tot.var = total.var, lambda = lambda)
				
				print('sum of current variances')
				print(sum(correction.coefs$orig.var.comps))
				
				
				# recompute scores
				tmp.var.comps <- correction.coefs$orig.var.comps / sum(correction.coefs$orig.var.comps)
				ptve <-  tmp.var.comps['var.X.Psi.Gamma']
				ptve.recomputed <- compute.ptve(data, true.model, test.set.indexes)
				print('sum of total pheno variance')
				print(sum(apply(data$phenotypes,2,var)))
				
				
				PTVE.ACCEPT <- TRUE
				
				print(sum(correction.coefs$orig.var.comps))
				if ( ((ptve/X.ptve) > 1.05) || ((ptve/X.ptve) < 0.95) ) {
					PTVE.ACCEPT <- FALSE
				}
				

				if (  (  (tmp.var.comps['var.Sigma'] / res.var.ptve) > 1.05) ||  (  (tmp.var.comps['var.Sigma'] / res.var.ptve) < 0.95) ) {
					
					PTVE.ACCEPT <- FALSE
				}
				
				if (PTVE.ACCEPT) {
					print(coef.iter)
					break;
				}
				
				
				
			}
		
		}
		
	}
	
	print('ptve aim:')
	print(X.ptve)
	print('exact ptve in toy data')
	print(ptve)
	print('final portions for variance')
	print(correction.coefs$orig.var.comps / sum(correction.coefs$orig.var.comps))
	true.model$portions.of.variance <- correction.coefs$orig.var.comps / sum(correction.coefs$orig.var.comps)
	true.model$total.var <- sum(apply(data$phenotypes,2,var))

	# current ptve is slighly theoretical, compute actual value
	ptve.recomputed <- compute.ptve(data, true.model, test.set.indexes)
	ptve <- ptve.recomputed
	
	# remove confounders so that they will not be used
	data$confounders <- NA	
	
	
	
	return(list(true.model = true.model, full.data = data, ptve = ptve))
	
	
}

