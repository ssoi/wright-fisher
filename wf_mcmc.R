
wf_mh_step <- function(data, prop, mcmc) {
	list(loglik=, )
	return()	
}

# initialize parameters for MCMC
numIters <- 1e5
mcmc <- list() # proposal, log-likelihood, and other deets for MCMC
mcmc[["B"]] <- 100 # number of replicates per simulation
mcmc[["sigma.m"]] <- 0.001 # s.d. for migration proposal
mcmc[["sigma.Ne"]] <- 0.001 # s.d. for migration proposal
mcmc[["sigma.Tm"]] <- 1 # s.d. for time of migration proposal
mcmc[["sigma.Tdiv"]] <- 1 # s.d. for time of divergence proposal
pars <- list() # samples from posterior distribution of samples
pars[["Ne"]] <- matrix(nrow=numIters, ncol=3)
pars[["m"]] <- matrix(nrow=numIters, ncol=2)
pars[["alpha"]] <- matrix(nrow=numIters, ncol=4)
pars[["Tdiv"]] <- pars[["Tm"]] <- numeric(numIters)
pars$m[1,] <- rep(0,2)
pars$alpha[1,] <- rep(0.1, 4)

for(iter in 2:numIters) {
		
}
