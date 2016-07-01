# copyright by the authors

check.mcmc.result <- function(context, mcmc.output, name, plot.path = NULL, plot.title = NULL) {
  # a function to check a MCMC result
  #
  # "name" is the name of the variable
  # "context" contains the true variable values
  # mcmc.output contains the "traces" variable
  
  if (name=='coefMat') {
    
    # Compare the coefficient matrices (obtainable from Psi and Gamma)
    
    res <- compute.coef.matrix.mcmc.estimate(mcmc.output)
    coef.matrix.mean <- res$coef.matrix.mean
    coef.matrix.std <- res$coef.matrix.std
    
    correct <- context$Psi %*% context$Gamma
    comp <- compare.matrices(correct=correct, est.mean=coef.matrix.mean, est.std=coef.matrix.std)
    
    
    plot.matrix.comparison(correct=correct, est.mean=coef.matrix.mean, plot.path = plot.path, plot.title)
    
    
  } 
  return(comp)
}