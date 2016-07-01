# copyright by the authors
compute.gamma.terms <- function(brr.rank, Gamma, Gamma.local.shrinkage) {
  # helper function for sampling star.taus from their posterior
	
	gamma.term.matr <- array(NA, dim= c(brr.rank, ncol(Gamma)))
	
	
		
	for (cl in 1:ncol(Gamma)) {
		
		Gamma.tmp.term <- Gamma[, cl, drop=F] %*% t(Gamma[, cl, drop=F])
		Gamma.tmp.term <- diag(Gamma.tmp.term)
		gamma.term.matr[, cl] <- Gamma.local.shrinkage[, cl] * Gamma.tmp.term
	}
		
	gamma.terms <- rowSums(gamma.term.matr)	
	return(gamma.terms)
	
	
}
