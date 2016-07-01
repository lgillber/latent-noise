# copyright by the authors
compute.coef.matrix.mcmc.estimate <- function(mcmc.output) {
  
  # this function computes the estimate of the regression 
  # coefficient matrix from the posterior samples
  
  n.iter <- length(mcmc.output$traces[['Gamma']])
  n.pheno <- ncol(mcmc.output$traces$Gamma[[1]])
  n.snps <- nrow(mcmc.output$traces$Psi[[1]])	
  
  coef.matrix.samples <- array(0, dim=c(n.snps, n.pheno, n.iter))
  
  for (i in 1:n.iter) {
    coef.matrix.samples[,,i] <- mcmc.output$traces$Psi[[i]] %*% mcmc.output$traces$Gamma[[i]]
  }
  coef.matrix.mean <- apply(coef.matrix.samples, 1:2, mean)
  coef.matrix.std <- apply(coef.matrix.samples, 1:2, sd)
  
  return(list(coef.matrix.mean=coef.matrix.mean, coef.matrix.std=coef.matrix.std))
}



