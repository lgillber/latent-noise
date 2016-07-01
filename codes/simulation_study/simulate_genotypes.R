# copyright by the authors
simulate.genotypes <- function(n.patients, n.snps, mafs = NULL) {
  #
  # A function for simulating independent SNP genotypes.
  #
  # Inputs:
  #   n.patients: number of patients
  #   n.snps: number of genotypes
  #
  # Outputs:
  #   genotypes (unnormalized)
  
  if (is.null(mafs)) {
    mafs <- runif(n=n.snps, min=0.05, 0.5)
  }
  
  prob.minor.homozygote <- mafs^2
  prob.heterozygote <- 2* mafs * (1-mafs)
  prob.major.homozygote <- (1-mafs)^2
  
  table <- matrix(runif(n=n.patients*n.snps, min=0, max=1),nrow=n.patients, ncol=n.snps)
  
  genotypes <- matrix(2, nrow=n.patients, ncol=n.snps)
  heterozygote <- t(t(table)>prob.minor.homozygote & t(table)<prob.minor.homozygote+prob.heterozygote)
  genotypes[heterozygote]=1
  major.homozygote <- t(t(table)>prob.minor.homozygote+prob.heterozygote)
  genotypes[major.homozygote]=0
  
  #return(list(genotypes=genotypes, mafs=mafs))
  return(genotypes)
}
