# copyright by the authors
compare.matrices <- function(correct, est.mean, est.std) {
  # a function to compare a mcmc estimates for a matrix of
  # parameters with the values used to generate the data
  comp <- list()
  comp$mse <- mean((correct-est.mean)^2)
  lower.bound <- est.mean-2*est.std
  upper.bound <- est.mean+2*est.std
  num.in.interval <- length(which(correct>lower.bound & correct<upper.bound))
  comp$prop.in.interval <- num.in.interval / prod(dim(correct))
  comp$correct <- correct
  comp$est.mean <- est.mean
  comp$est.std <- est.std
  # the correlation between the correct and the estimated parameter
  # matrices after vectorization
  comp$correlation <- cor(c(correct), c(est.mean))
  return(comp)
}