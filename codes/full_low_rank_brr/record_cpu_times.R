# copyright by the authors
record.cpu.times <- function(cpu.times=NULL, conf.res=NULL, fa.res=NULL, brr.res=NULL) {
  # The function updates a vector specifying CPU-times taken by 
  # updates of each variable.
  #
  # Inputs:
  #	cpu.times: a vector where the total time spent at each update
  #              is written
  #
  #	conf.res, fa.res, brr.res: lists, each having cpu.times as one 
  #              of the components. The cpu.times is a list where the
  #              elements tell the cpu time taken by updates of FA,
  #              BRR, or confounder parts of the models, respectively.
	
	if (is.null(cpu.times)) {
		# Initialize trace

		var.names <- NULL
		var.names <- names(brr.res$cpu.times)
		var.names <- c(var.names, names(fa.res$cpu.times))
		var.names <- c(var.names, names(conf.res$cpu.times))
		cpu.times <- rep(0, length(var.names))
		names(cpu.times) <- var.names

	}


	for (res in c('conf.res','fa.res','brr.res')) {
		
		if (!is.null(res)) {
			cpu.times.to.add <- eval(parse(text=paste(res, '$cpu.times', sep='')))
			var.names <- names(cpu.times.to.add)
			
			for (name in var.names) {
				cpu.times[name] <- cpu.times[name] + cpu.times.to.add[name]
			}
		}
	}

	return(cpu.times)
}
