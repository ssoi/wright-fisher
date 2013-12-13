library(truncnorm)
source("wf_sim.R")

# density of truncated normal
dtnorm <- function(x, mean, sd, a, b) 
	dnorm(x, mean, sd)/(sd*(pnorm(b, mean, sd) - pnorm(a, mean, sd)))

# propose new parameters for Metropolis-Hastings
wf_mh_propose <- function(param, cntl) {	
	proposal <- param
	if(cntl$currIter %% 4 == 1) {
		proposal$M1 <- rtruncnorm(n=1, mean=param$M1, sd=cntl$sigma.M, 
			a=0, b=1)
		#proposal$M1 <- rtruncnorm(n=1, mu=param$M1, sigma=cntl$sigma.M, 
		#	lower=0, upper=1)
	} else if(cntl$currIter %% 4 == 2) {
		proposal$M2 <- rtruncnorm(n=1, mean=param$M2, sd=cntl$sigma.M, 
			a=0, b=1)
		#proposal$M2 <- rtruncnorm(n=1, mu=param$M2, sigma=cntl$sigma.M, 
		#	lower=0, upper=1)
	} else if(cntl$currIter %% 4 == 2) {
		proposal$M1 <- rtruncnorm(n=1, mean=param$M1, sd=cntl$sigma.M, 
			a=0, b=1)
	} else if(cntl$currIter %% 4 == 3) {
		proposal$Ne1 <- rtruncnorm(n=1, mean=param$Ne1, sd=cntl$sigma.M, 
			a=1e2, b=1e5)
	} else if(cntl$currIter %% 4 == 0) {
		proposal$Ne2 <- rtruncnorm(n=1, mean=param$Ne2, sd=cntl$sigma.M, 
			a=1e2, b=1e5)
	}
	return(proposal)
}

# Metropolis-Hastings ABC/BIL step
wf_mh_step <- function(obs, curr, prev, cntl, init=F) {
	t <- c(1L, curr$Tmig, curr$Tdiv)
	b <- cntl$B
	ne <- c(curr$Ne1, curr$Ne1, curr$Ne2)
	m <- c(curr$M1, curr$M2)
	a <- c(curr$alpha1, curr$alpha2, curr$alpha3, curr$alpha4)
	sim <- wf_sim(param=list(T=t, B=b, Ne=ne, A=a, M=m))
	if(cntl$currIter %% 2 == 1) {			
		model <- nls(A ~ C0*exp(-C1*r-C2*r^2), data=data.frame(A=sim$A1, r=rec),
			start=list(C0=1, C1=1, C2=1))
		curr$loglik <- dnorm(x=sum(obs$A1 - fitted(model)), mean=0,
			sd=sd(sim$A1 - fitted(model)), log=TRUE)
		q01 <- log(dtnorm(x=curr$M1, mean=prev$M1, sd=cntl$sigma.M, a=0, b=1))
		q10 <- log(dtnorm(x=prev$M1, mean=curr$M1, sd=cntl$sigma.M, a=0, b=1))
	} else if(cntl$currIter %% 4 == 2) {
		model <- nls(A ~ C0*exp(-C1*r-C2*r^2), data=data.frame(A=sim$A2, r=rec),
			start=list(C0=1, C1=1, C2=1))
		curr$loglik <- dnorm(x=sum(obs$A2 - fitted(model)), mean=0,
			sd=sd(sim$A2 - fitted(model)), log=TRUE)
		q01 <- log(dtnorm(x=curr$M2, mean=prev$M2, sd=cntl$sigma.M, a=0, b=1))
		q10 <- log(dtnorm(x=prev$M2, mean=curr$M2, sd=cntl$sigma.M, a=0, b=1))
	} else if(cntl$currIter %% 4 == 3) {	
		model <- nls(D ~ C0*exp(-C1*r-C2*r^2), data=data.frame(A=sim$D1^2, r=rec),
			start=list(C0=1, C1=1, C2=1))
		curr$loglik <- dnorm(x=sum(obs$D1^2 - fitted(model)), mean=0,
			sd=sd(sim$D1^2 - fitted(model)), log=TRUE)
		q01 <- log(dtnorm(x=curr$Ne1, mean=prev$Ne1, sd=cntl$sigma.Ne, a=0, b=1))
		q10 <- log(dtnorm(x=prev$Ne1, mean=curr$Ne1, sd=cntl$sigma.Ne, a=0, b=1))	
	} else if(cntl$currIter %% 4 == 0) {	
		model <- nls(D ~ C0*exp(-C1*r-C2*r^2), data=data.frame(A=sim$D2^2, r=rec),
			start=list(C0=1, C1=1, C2=1))
		curr$loglik <- dnorm(x=sum(obs$D2^2 - fitted(model)), mean=0,
			sd=sd(sim$D2^2 - fitted(model)), log=TRUE)
		q01 <- log(dtnorm(x=curr$Ne2, mean=prev$Ne2, sd=cntl$sigma.Ne, a=0, b=1))
		q10 <- log(dtnorm(x=prev$Ne2, mean=curr$Ne2, sd=cntl$sigma.Ne, a=0, b=1))
	}
	if(init) return(curr)
	alpha <- min(0, curr$loglik + q10 - prev$loglik - q01)
	if(runif(1) < curr$loglik) {
		return(list(acc=TRUE, param=curr))
	} else {
		return(list(acc=FALSE, param=prev))
	}
}

# generate "pseudo-observed dataset"
# alternatively read in data set (future)
rec <- seq(1e-5, 1e-2, length=50)
pod <- wf_sim(param=list(T=c(1L, 150L, 200L), B=20L, Ne=c(1e4L, 1e4L, 3e3L), 
	A=c(2, 1, 1, 2), M=c(0.01, 0.0)))

# define control for control MCMC	
cntl <- list() # proposal, log-likelihood, and other deets for MCMC
cntl$currIter <- 1 # current iteration of MCMC chain
cntl$B <- 20L # number of replicates per simulation
cntl$sigma.M <- 0.01 # s.d. for migration proposal
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
param[1,] <- wf_mh_step(obs=pod, curr=param[1,], prev=param[1,], cntl=cntl, init=T)

# run chain
for(iter in 3556:numIters) {
	cntl$currIter <- iter
	proposal <- try(wf_mh_propose(param[iter-1,], cntl))
	repeat {
		mh <- wf_mh_step(obs=pod, curr=proposal, prev=param[iter-1,], cntl=cntl)
		if(class(mh) != "try-error") break
	}	
	if(mh$acc) {
		param[iter,] <- mh$param
	} else {
		param[iter,] <- param[iter-1,] 
	}
	if(iter %% 100 == 0)
		print(paste(iter, ":", proposal$M1, proposal$M2, "|", 
			paste(param[iter,], collapse=" ")))
}
	if(init) return(curr)
	alpha <- min(0, curr$loglik + q10 - prev$loglik - q01)
	if(runif(1) < curr$loglik) {
		return(list(acc=TRUE, param=curr))
	} else {
		return(list(acc=FALSE, param=prev))
	}
}

# generate "pseudo-observed dataset"
# alternatively read in data set (future)
rec <- seq(1e-5, 1e-2, length=50)
pod <- wf_sim(param=list(T=c(1L, 150L, 200L), B=20L, Ne=c(1e4L, 1e4L, 3e3L), 
	A=c(2, 1, 1, 2), M=c(0.01, 0.0)))

# define control for control MCMC	
cntl <- list() # proposal, log-likelihood, and other deets for MCMC
cntl$currIter <- 1 # current iteration of MCMC chain
cntl$B <- 20L # number of replicates per simulation
cntl$sigma.M <- 0.01 # s.d. for migration proposal
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
param[1,] <- wf_mh_step(obs=pod, curr=param[1,], prev=param[1,], cntl=cntl, init=T)

# run chain
for(iter in 3556:numIters) {
	cntl$currIter <- iter
	proposal <- try(wf_mh_propose(param[iter-1,], cntl))
	repeat {
		mh <- wf_mh_step(obs=pod, curr=proposal, prev=param[iter-1,], cntl=cntl)
		if(class(mh) != "try-error") break
	}	
	if(mh$acc) {
		param[iter,] <- mh$param
	} else {
		param[iter,] <- param[iter-1,] 
	}
	if(iter %% 100 == 0)
		print(paste(iter, ":", proposal$M1, proposal$M2, "|", 
			paste(param[iter,], collapse=" ")))
}
