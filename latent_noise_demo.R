#
# copyright by the authors
#
#################################################################
#                                                               #
# This is a demo about latent-noise and the method presented in # 
# "Multiple Output Regression with Latent Noise", accepted for  #
# publication at the Journal of Machine Learning Research by    #
# Gillberg et al.                                               #
#                                                               #
#################################################################
rm(list=ls())

# the following function source_directory is by 
# Mehmet Gonen (mehmet.gonen@gmail.com)
source_directory <- function(path) {
  
  files <- sort(dir(path, "\\.[rR]$", full.names = TRUE))
  lapply(files, source, chdir = TRUE)
}


################################
# load functions and libraries #
################################

current.path <- '/home/lgillber/latent_noise_2015_no_highlights/latent_noise_code'
#current.path <- 'addMyPathHere'
if (current.path == 'addMyPathHere') {
  stop('Add the path to the folder containing this file as current.path')
}

#
code.path <- paste0(current.path, '/codes')


# load codes
# data simulation etc
source_directory(paste0(code.path, '/simulation_study/'))
# combining all parts of the model (latent-noise, low-rank noise ...)
source_directory(paste0(code.path, '/full_low_rank_brr/'))
# code for the infinite BRRR
source_directory(paste0(code.path, '/infinite_brr/'))
# code for the independent structured noise (H \Lambda) 
# and residual variances 
source_directory(paste0(code.path, '/sparse_fa/'))
# general
source_directory(paste0(code.path, '/common/'))



library(coda)
library(MASS)




#####################
# Generate toy data #
#####################

set.seed(100)
n.pheno <- 60
n.snps <- 30
n.train <- 500
n.test <- 5000

training.set.indexes <- 1:n.train  
test.set.indexes <- (1:n.test+n.train)

# The variance of Y is partitioned according to the following:
#
# 3 % is explained by the covariates (X)
# 20 % is explained by unstructured residual variance explains
# the remaining (77 %) is explained by the structured noise
#
# of the remaining 77 %,
# 100 % is latent noise (1-lambda) and 0 % (lambda) is independent noise
#
# These values are theoretical, the actual values can differ slightly

generative.model.and.data <- simulate.model.and.data(n.pheno=n.pheno, n.snps=n.snps, n.patients = (n.train+n.test), test.set.indexes = test.set.indexes, X.ptve = 0.03, res.var.ptve = 0.2, lambda = 0)


# split generated data into training and test set
data <- generative.model.and.data$full.data
data$genotypes <- data$genotypes[training.set.indexes,]
data$phenotypes <- data$phenotypes[training.set.indexes,]
data$crossprod.genotypes <- NULL
  
test.data <- generative.model.and.data$full.data
test.data$genotypes <- test.data$genotypes[test.set.indexes,]
test.data$phenotypes <- test.data$phenotypes[test.set.indexes,]
test.data$crossprod.genotypes <- NULL


############################
# Learn latent-noise BRRR  #
############################


# rank is sampled
brr.rank <- NULL

# latent signal-to-noise-ratio (latent SNR)
# Omega.coef: 1/(latent SNR)
# we are dealing with weak effects, this parameterization seems 
# convenient
Omega.coef <- c(7.5)





learnt.params <- list()

# sampling parameters
n.iter <- 1000
burnin <- 0.5
thin<-10

# BRRR rank used for inititialization
if (is.null(brr.rank)) brr.rank.init <- 3 else brr.rank.init <- brr.rank


######################
# initilialize model #
######################

set.seed(300)
init.model <- initialize.from.prior(n.pheno= n.pheno, n.snps=n.snps, n.patients=n.train, fa.rank=3, brr.rank=brr.rank.init,  a.sigma = 2.2, b.sigma = 0.5)

# parameters related to sampling the rank: how long should it be updated
init.model$brr$prior$alpha0 <- -2
init.model$brr$prior$alpha1 <- (log(0.1) - init.model$brr$prior$alpha0) / (n.iter*burnin)


# set BRRR shrinkage parameters. Note that in the notation of the
# paper these are a1 and a2 (and vice versa)
init.model$brr$context$a3a4 <- c(10, 4.1)

#########################
# latent-noise variance #
#########################

# Set the variance of the latent noise to match the prior assumption
# about the a priori latent signal-to-noise ratio

# first compute prior variance of X\Psi, then set
# prior variance of \Omega to match it when multiplied with
# a factor of "Omega.coef"
# for how to compute var(X %*% Psi), see e.g.
# http://en.wikipedia.org/wiki/Variance#Basic_properties



# each component has same variance
# do following var calculation 1 comp at a time and sum
prior.var.X.Psi <- 0
prior.vars.Psi <- diag(ncol(data$genotypes))

for (rank.tmp in 1:init.model$brr$context$brr.rank) {

  prior.var.X.Psi <- prior.var.X.Psi + sum(diag(data$genotypes %*% prior.vars.Psi %*% t(data$genotypes)))
}


# all elements of Omega assumed to have the same variance
init.model$brr$context$latent.noise.var <- Omega.coef * prior.var.X.Psi / (n.train * brr.rank.init)
init.model$brr$context$Omega <- matrix(rnorm(brr.rank.init*n.train, sd = sqrt(init.model$brr$context$latent.noise.var)), nrow=n.train, ncol=brr.rank.init)

# random initialization with excess variance
init.model$brr$context$Psi[,] <- rnorm(n.snps*brr.rank.init, sd=1)
init.model$brr$context$Gamma[,] <- rnorm(brr.rank.init*n.pheno, sd=1)


  
  

  
set.seed(125)
mcmc.output <- gibbs.full.low.rank.brr(model=init.model, data=data, n.iter=n.iter, thin=thin, fixed.brr.rank=brr.rank, fa.vars.to.record=c('variances','local.shrinkage','rank','a1a2','deltas'))



# to use the independent-noise BRRR model
#
#  in the initialization, set:
#    init.model$brr$context$Omega <- NULL
#
#  give the following arguments to gibbs.full.low.rank.brr:
#    ind.struct.noise = TRUE
#    latent.noise = FALSE





print('model learnt')	

  
# study MCMC chain: with large number of samples and data points, 
# this will converge to the model used to generate the data
tmp.mcmc.output <- remove.burnin(mcmc.output=mcmc.output, burnin=burnin*round(n.iter/thin))

# compare learnt parameters with the values used to generate toy data
# 
comp <- check.mcmc.result(context=generative.model.and.data$true.model$brr$context, mcmc.output=tmp.mcmc.output, name='coefMat', plot.path = '/home/lgillber/latent_noise_2015_no_highlights/latent_noise_code/', plot.title = 'coefficient matrix')

# test set performance
test.scores <- compute.prediction.error(test.data, tmp.mcmc.output, burnin = burnin)
print(test.scores)

