# copyright by the authors
simulate.new.brr.factor <- function(star.deltas, local.shrinkage.nu, a4, n.pheno, n.snps) {
	
  # Column to be added is never the first.
	new.star.delta <- rgamma(n=1, shape=a4, rate=1)  
	new.star.tau <- prod(star.deltas) * new.star.delta
	
	# Psi
	
	
	# Psi is not shrunk
	new.Psi.col <- rnorm(n=n.snps, mean=0, sd=1)
		
	
	
	
	#n.pheno.clusters <- length(unique(output.clustering))
	
	new.Gamma.local.shrinkage.parameters <- rgamma(n=n.pheno, shape=local.shrinkage.nu/2, rate=local.shrinkage.nu/2)
	
	new.Gamma.local.shrinkage.row <- new.Gamma.local.shrinkage.parameters#[output.clustering]
	

		# Gamma is doubly shrunk
		new.Gamma.row <- rnorm(n=n.pheno, mean=0, sd=1/sqrt(new.Gamma.local.shrinkage.row * (new.star.tau)^2))		
		
	
	
	return(list(new.star.delta=new.star.delta, new.star.tau=new.star.tau, new.Psi.col=new.Psi.col, new.Gamma.local.shrinkage.row=new.Gamma.local.shrinkage.row, new.Gamma.row=new.Gamma.row))
	
}
