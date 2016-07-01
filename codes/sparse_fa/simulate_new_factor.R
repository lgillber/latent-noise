# copyright by the authors
simulate.new.factor <- function(a.sigma, b.sigma, deltas, local.shrinkage.nu, a2, n.pheno, n.patients, prior.var.Eta) {
	#
	# This is a helper function of "sparse.fa.gibbs". The
	# function is used when adding columns during the 
	# rank adaptation
	#

	# a2<-10; n.pheno<-8; n.patients<-100; local.shrinkage.nu<-3; a.sigma<-1; b.sigma<-0.3; deltas<-c(4,5)

	new.delta <- rgamma(n=1, shape=a2, rate=1)  # Column to be added is never the first.
	new.tau <- prod(deltas) * new.delta
	
	new.local.shrinkage.column <- rgamma(n=n.pheno, shape=local.shrinkage.nu/2, rate=local.shrinkage.nu/2)
	
	new.Lambda.col <- rnorm(n=n.pheno, mean=0, sd=1/sqrt(new.local.shrinkage.column * new.tau))
	
	#new.precision <- rgamma(n=1, shape=a.sigma, rate=b.sigma)
	#new.variance <- 1/new.precision

	new.Eta.column <- rnorm(n=n.patients, mean=0, sd=sqrt(prior.var.Eta))

	#return(list(new.delta=new.delta, new.tau=new.tau, new.local.shrinkage.column=new.local.shrinkage.column, new.Lambda.col=new.Lambda.col, new.precision=new.precision, new.variance=new.variance, new.Eta.column=new.Eta.column))

	return(list(new.delta=new.delta, new.tau=new.tau, new.local.shrinkage.column=new.local.shrinkage.column, new.Lambda.col=new.Lambda.col, new.Eta.column=new.Eta.column))
}

