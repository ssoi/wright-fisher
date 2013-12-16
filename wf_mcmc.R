library(splines, quietly=TRUE)
library(tmvtnorm, quietly=TRUE)
source("wf_sim.R")

# density of truncated normal
dtnorm <- function(x, mean, sd, a, b) 
	dnorm(x, mean, sd)/(sd*(pnorm(b, mean, sd) - pnorm(a, mean, sd)))

# return coefficients from linear model
auxillary <- function(y, x) {
	fit <- lm(y ~ bs(x, df = 3))
	return(list(Y=fitted(fit), res=y-fitted(fit), coef=coef(fit)))
}

# propose new parameters for Metropolis-Hastings
wf_mh_propose <- function(param, cntl) {
	proposal <- param
	if(cntl$currIter %% 4 == 1) {
		proposal$M1 <- rtmvnorm(n=1, mean=param$M1, sigma=cntl$sigma2.M,
			upper=1, lower=0)
	} else if(cntl$currIter %% 4 == 2) {
		proposal$M2 <- rtmvnorm(n=1, mean=param$M2, sigma=cntl$sigma2.M,
			upper=1, lower=0)
	} else if(cntl$currIter %% 4 == 3) {
		proposal$Ne1 <- rtmvnorm(n=1, mean=param$Ne1, sigma=cntl$sigma2.Ne, 
			upper=1e5, lower=1e2)
	} else if(cntl$currIter %% 4 == 0) {
		proposal$Ne2 <- rtmvnorm(n=1, mean=param$Ne2, sigma=cntl$sigma2.Ne, 
			upper=1e5, lower=1e2)
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
	# simulate data with proposed parameters
	sim <- wf_sim(param=list(T=t, B=b, Ne=ne, A=a, M=m))
	# fit curve and calculate likelihood based on multivariate normal
	mu <- colSums(simplify2array(sim) - obs$aux)
	curr$loglik <- dmvnorm(x=mu, sigma=obs$sigma, log=T)
	#	calculate appropriate transition densities
	if(cntl$currIter %% 4 == 1) {
		q01 <- log(dtnorm(x=curr$M1, mean=prev$M1, sd=cntl$sigma.M, a=0, b=1))
		q10 <- log(dtnorm(x=prev$M1, mean=curr$M1, sd=cntl$sigma.M, a=0, b=1))
	} else if(cntl$currIter %% 4 == 2) {
		q01 <- log(dtnorm(x=curr$M2, mean=prev$M2, sd=cntl$sigma.M, a=0, b=1))
		q10 <- log(dtnorm(x=prev$M2, mean=curr$M2, sd=cntl$sigma.M, a=0, b=1))
	} else if(cntl$currIter %% 4 == 3) {	
		q01 <- log(dtnorm(x=curr$Ne1, mean=prev$Ne1, sd=cntl$sigma.Ne, a=1e2, b=1e5))
		q10 <- log(dtnorm(x=prev$Ne1, mean=curr$Ne1, sd=cntl$sigma.Ne, a=1e2, b=1e5))	
	} else if(cntl$currIter %% 4 == 0) {	
		q01 <- log(dtnorm(x=curr$Ne2, mean=prev$Ne2, sd=cntl$sigma.Ne, a=1e2, b=1e5))
		q10 <- log(dtnorm(x=prev$Ne2, mean=curr$Ne2, sd=cntl$sigma.Ne, a=1e2, b=1e5))
	}
	if(init) return(curr)
	alpha <- min(0, curr$loglik + q10 - prev$loglik - q01)
	if(log(runif(1)) < alpha) {
		return(list(acc=TRUE, param=curr))
	} else {
		return(list(acc=FALSE, param=prev))
	}
}

# generate "pseudo-observed dataset" and fit auxillary model
# alternatively read in data set (future)
rec <- seq(1e-5, 1e-2, length=50)
pod <- wf_sim(param=list(T=c(1L, 300L, 400L), B=50L, Ne=c(1e4L, 1e4L, 3e3L), 
	A=c(2, 1, 1, 2), M=c(0.01, 0.0)))
aux <- lapply(pod, auxillary, x=rec)
pod$aux <- simplify2array(sapply(aux, "[[", "Y"))
pod$sigma <- cov(simplify2array(sapply(aux, "[[", "res")))

# define control for control MCMC	
cntl <- list() # proposal, log-likelihood, and other deets for MCMC
cntl$numIters <- 1e4 # length of chain
cntl$currIter <- 1 # current iteration of MCMC chain
cntl$B <- 50L # number of replicates per simulation
cntl$sigma.M <- 0.01 # variance for migration proposal
cntl$sigma.Ne <- 1e3 # variance for migration proposal
cntl$sigma2.M <- cntl$sigma.M^2 # variance for migration proposal
cntl$sigma2.Ne <- cntl$sigma.Ne^2 # variance for migration proposal
cntl$sigma.Tm <- 1 # variance for time of migration proposal
cntl$sigma.Tdiv <- 1 # variance for time of divergence proposal

# define parameters for posterior sampling 

param <- data.frame(
	Ne1=rep(NA, cntl$numIters), 
	Ne2=rep(NA, cntl$numIters), 
	M1=rep(NA, cntl$numIters), 
	M2=rep(NA, cntl$numIters), 
	alpha1=rep(NA, cntl$numIters), 
	alpha2=rep(NA, cntl$numIters), 
	alpha3=rep(NA, cntl$numIters), 
	alpha4=rep(NA, cntl$numIters), 
	Tmig=rep(NA, cntl$numIters), 
	Tdiv=rep(NA, cntl$numIters),
	loglik=rep(NA, cntl$numIters)
)

# initialize MCMC chain
param$Ne1[1] <- 1e3L
param$Ne2[1] <- 1e3L
param$M1[1] <- param$M2[1] <- 0
param$alpha1[1] <- param$alpha4[1] <- 2
param$alpha2[1] <- param$alpha3[1] <- 1
param$Tmig[1] <- 300L 
param$Tdiv[1] <- 400L

# calculate log-likelihood of initial state of chain
param[1,] <- wf_mh_step(obs=pod, curr=param[1,], prev=param[1,], cntl=cntl, init=T)

# run chain
for(iter in 2:cntl$numIters) {
	cntl$currIter <- iter
	repeat {
		proposal <- wf_mh_propose(param[iter-1,], cntl)
		mh <- try(wf_mh_step(obs=pod, curr=proposal, prev=param[iter-1,], cntl=cntl))
		if(class(mh) != "try-error") break
	}	
	if(mh$acc) {
		param[iter,] <- mh$param
	} else {
		param[iter,] <- param[iter-1,] 
	}
	cntl$acc[iter-1] <- mh$acc
	if(iter %% 10 < 4)
		print(paste(iter, ":", paste(round(proposal[1:4],4), collapse=" "), "|", 
			paste(round(param[iter,c(1:4, 11)],4), collapse=" ")))
}
