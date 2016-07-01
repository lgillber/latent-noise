# copyright by the authors
normalize.ln.prop.to.distr <- function(ln.prop.to.pr, proposal) {
	# assuming evenly spaced grid!
	
	grid.density <- proposal[3] - proposal[2]
	
	ln.prop.to.pr <- ln.prop.to.pr - max(ln.prop.to.pr)
	ln.prop.to.pr <- exp(ln.prop.to.pr)
	ln.prop.to.pr <- ln.prop.to.pr / sum(ln.prop.to.pr)
	ln.prop.to.pr <- cumsum(ln.prop.to.pr)
	ind.to.select <- which(runif(1) < ln.prop.to.pr)[1]
	
	# make sure not to give negative values!
	if (ind.to.select == 1) {
		to.return <- proposal[ind.to.select] + grid.density*(runif(1)*0.5)
	} else {
		# the regular approximative gibbs case: add uniformly distributed
		# jitter of grid density
		to.return <- proposal[ind.to.select] + grid.density*(runif(1)-0.5)	
	}
	
	return(to.return)
}
