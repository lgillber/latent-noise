# copyright by the authors
compute.amount.total.variance.explained <- function(genotypes, Psi, Gamma) {
	#
	# Function for computing the amount of total variation
	# explained by the reduced rank part of the model, given
	# fixed parameters Psi and Gamma.
	#
	# Inputs:
	#	genotypes: matrix of centered genotypes
	#
	#	Psi, Gamma: low rank representation of the 
	#		coefficient matrix
	#
	#
	# Outputs:
	#	amount.total.var.explained: a scalar specifying the amount
	#		of variance explained by the BRR part of the model.	
	#

	aux <- genotypes %*% Psi %*% Gamma
	amount.total.var.explained <- sum(apply(aux,2,var))

}
