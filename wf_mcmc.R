library(truncnorm)

# propose new parameters for Metropolis-Hastings
wf_mh_propose <- function(param, cntl) {
	proposal <- param
	if(cntl$currIter %% 2 == 1) {
		proposal$M1 <- rtruncnorm(n=1, mean=param$M1, sd=cntl$sigma.m, 
			a=0, b=1)
		#proposal$M1 <- rtruncnorm(n=1, mu=param$M1, sigma=cntl$sigma.m, 
		#	lower=0, upper=1)
	} else {
		proposal$M2 <- rtruncnorm(n=1, mean=param$M2, sd=cntl$sigma.m, 
			a=0, b=1)
		#proposal$M2 <- rtruncnorm(n=1, mu=param$M2, sigma=cntl$sigma.m, 
		#	lower=0, upper=1)
	}
	return(proposal)
}

# Metropolis-Hastings ABC/BIL step
wf_mh_step <- function(obs, param, cntl) {
	t <- c(param$Tmig, param$Tdiv)
	b <- cntl$B
	ne <- c(param$Ne1, param$Ne2)
	m <- c(param$M1, param$M2)
	a <- c(param$alpha1, param$alpha2, param$alpha3, param$alpha4)
	sim <- wf_sim(param=list(T=t, B=b, Ne=ne, A=a, M=m))
	if(cntl$currIter %% 2 == 1) {			
		model <- nls(A ~ C0*exp(-C1*r-C2*r^2), data=data.frame(A=sim$A1, r=rec),
			start=list(C0=1, C1=1, C2=1))
		loglik <- dnorm(x=sum(obs$A1 - fitted(model)), mean=0,
			sd=sd(sim$A1 - fitted(model)), log=TRUE)
	} else {
		model <- nls(A ~ C0*exp(-C1*r), data=data.frame(A=sim$A2, r=rec),
			start=list(C0=1, C1=1))
		loglik <- dnorm(x=sum(obs$A2 - fitted(model)), mean=0,
			sd=sd(sim$A2 - fitted(model)), log=TRUE)
	}
	alpha <- min(0, loglik - param$loglik
	if(runif(1) < loglik) {
		return(list(acc=TRUE, param=))
	} else {
		return(list(acc=FALSE, param=loglik))
	}
}

# generate "pseudo-observed dataset"
# alternatively read in data set (future)
rec <- seq(1e-5, 1e-2, length=100)
pod <- wf_sim(param=list(T=c(1L, 150L, 200L), B=100L, Ne=c(1e4L, 1e4L, 3e3L), 
	A=c(2, 1, 1, 2), M=c(0.01, 0.0)))

# define control for control MCMC
cntl <- list() # proposal, log-likelihood, and other deets for MCMC
cntl$currIter <- 1 # current iteration of MCMC chain
cntl$B <- 100 # number of replicates per simulation
cntl$sigma.m <- 0.001 # s.d. for migration proposal
cntl$sigma.Ne <- 0.001 # s.d. for migration proposal
cntl$sigma.Tm <- 1 # s.d. for time of migration proposal
cntl$sigma.Tdiv <- 1 # s.d. for time of divergence proposal

# define parameters for posterior sampling 
numIters <- 1e4 # length of chain
param <- data.frame(
	Ne1=rep(NA, numIters), 
	Ne2=rep(NA, numIters), 
	M1=rep(NA, numIters), 
	M2=rep(NA, numIters), 
	alpha1=rep(NA, numIters), 
	alpha2=rep(NA, numIters), 
	alpha3=rep(NA, numIters), 
	alpha4=rep(NA, numIters), 
	Tmig=rep(NA, numIters), 
	Tdiv=rep(NA, numIters),
	loglik=rep(NA, numIters)
)

# initialize MCMC chain
param$Ne1[1] <- 1e4L
param$Ne2[1] <- 3e3L
param$M1[1] <- param$M2[1] <- 0
param$alpha1[1] <- param$alpha4[1] <- 2
param$alpha2[1] <- param$alpha3[1] <- 1
param$Tmig[1] <- 150L 
param$Tdiv[1] <- 200L

# calculate log-likelihood of initial state of chain
mh <- wf_mh_step(obs=pod, param=param[1,], cntl=cntl)
mh$loglik[1] <- mh$loglik

# run chain
for(iter in 2:numIters) {
	cntl$currIter <- iter
	repeat {
		proposal <- try(wf_mh_propose(param[iter-1,], cntl))
		if(class(proposal) != "try-error") break
	}	
	mh <- wf_mh_step(obs=pod, param=proposal, cntl=cntl)
	if(mh$acc) {
		param[iter,] <- proposal
	} else {
		param[iter,] <- param[iter-1,] 
	}
}
