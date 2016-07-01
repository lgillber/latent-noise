# copyright by the authors
compute.residuals <- function(type, model, data) {
	#
	# Helper function for the Gibbs sampling. The function
	# computes the residuals 
	#
	# Inputs:
	#	type: One of 'brr', 'fa', or 'confounder', depending on
	#         which part of the model we are going to update
	#         using the obtained residuals.
	#
	#	model: The current model
	#	
	#	data: List with fields: 'genotypes', 'phenotypes', and
	#		  'confounders'
	#
	#
	# Outputs, a list containing:
	#	residuals: matrix of residuals

	phenotypes <- data$phenotypes
	if (type=='brr') {
		
		residuals <- phenotypes
		
		confounders <- data$confounders		
		if (!is.na(confounders)[1]) {
			A <- model$A
			residuals <- residuals - confounders %*% A
		} 
		
		if (!is.null(model$fa$context$Eta)) {
			Eta <- model$fa$context$Eta
			Lambda <- model$fa$context$Lambda
			# if using noise model
			residuals <- residuals - tcrossprod(Eta,Lambda)
		}

		
	} else if (type=='fa') {

 		A <- model$A
 		confounders <- data$confounders
		genotypes <- data$genotypes
 		Psi <- model$brr$context$Psi
 		Gamma <- model$brr$context$Gamma
		
		
		
		residuals <- phenotypes
		
		if (!is.na(confounders)[1]) {
			residuals <- residuals - confounders %*% A 
		} 
		
		if (is.null(model$brr$context$Omega)) {
			residuals <- residuals - genotypes %*% Psi %*% Gamma
		} else {
			Omega <- model$brr$context$Omega
			residuals <- residuals - (genotypes %*% Psi + Omega)  %*% Gamma 
		}

	} else if (type=='confounder') {
		
		print('warning: confounder -part of the code not tested!')
		# Omega modifications not done yet, see type='fa' case above
		genotypes <- data$genotypes
 		Psi <- model$brr$context$Psi
 		Gamma <- model$brr$context$Gamma
		Omega <- model$brr$context$Omega
		
		residuals <- phenotypes - (genotypes %*% Psi + Omega) %*% Gamma 
		
		if (!is.null(model$fa$Eta)) {
			
			Eta <- model$fa$context$Eta
			Lambda <- model$fa$context$Lambda
			residuals <- residuals - tcrossprod(Eta,Lambda)
			
		
		}

	} else {
		stop('Unknown residual type')
	}

	return(residuals)
}
