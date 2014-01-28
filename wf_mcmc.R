library(multicore, quietly=TRUE)
library(splines, quietly=TRUE)
library(tmvtnorm, quietly=TRUE)
source("wf_sim.R")

# density of truncated normal
dtnorm <- function(x, mean, sd, a, b) 
	dnorm(x, mean, sd)/(sd*(pnorm(b, mean, sd) - pnorm(a, mean, sd)))

# return coefficients from linear model
auxillary <- function(y, x) {
	fit <- nls(Y ~ c0*exp(c1*X+c2*X^2), data=list(Y=y, X=x), start=c(c0=1, c1=1, c2=2))
	return(list(Y=fitted(fit), res=y-fitted(fit), coef=coef(fit)))
}

loglik_sample <- function(m, s2, mu, sigma2, N) {
	T <- (m - mu)/sqrt(s2/N)
	a <- (N-1)/2
	b <- (N-3)/2
	s.loglik <- a*(log(N) - log(2*sigma2)) + b*log(s2) - N*s2/(2*sigma2) -lgamma(a)
	m.loglik <- dt(x=T, df=N-1, log=TRUE)
	return(sum(s.loglik + m.loglik))
}

# propose new parameters for Metropolis-Hastings
wf_mh_propose <- function(param, cntl) {
	proposal <- param
	mod <- cntl$currIter %% cntl$P
	if(mod == 0) {
		proposal$M <- rtmvnorm(n=1, mean=param$M, sigma=cntl$sigma2.M,
			upper=1, lower=0)
		} else if(mod == 1) { 
			proposal$Ne1 <- rtmvnorm(n=1, mean=param$Ne1, sigma=cntl$sigma2.Ne,
				upper=1e5, lower=1e2)
	} else if(mod == 2) {
		proposal$Ne2 <- rtmvnorm(n=1, mean=param$Ne2, sigma=cntl$sigma2.Ne, 
			upper=1e5, lower=1e2)	
	} else if(mod == 3) {
		proposal$Tdiv <- rtmvnorm(n=1, mean=param$Tdiv, sigma=cntl$sigma2.Tdiv, 
			upper=Inf, lower=1e1)
	} else if(mod == 4) {
		proposal$Tmig <- rtmvnorm(n=1, mean=param$Tmig, sigma=cntl$sigma2.Tmig, 
			upper=param$Tdiv-1, lower=1)	
	} else if(mod == 5) {
		alpha <- rtmvnorm(n=1, 
			mean=c(proposal$alpha1, proposal$alpha2, proposal$alpha3, proposal$alpha4),
			sigma=diag(4)*cntl$sigma2.alpha, lower=rep(0, 4))	
		proposal$alpha1 <- alpha[1]
		proposal$alpha2 <- alpha[2]
		proposal$alpha3 <- alpha[3]
		proposal$alpha4 <- alpha[4]
	}
	return(proposal)
}

# Metropolis-Hastings ABC/BIL step
wf_mh_step <- function(obs, curr, prev, cntl, plot=F, init=F) {
	t <- c(curr$Tmig, curr$Tdiv)
	b <- cntl$B
	ne <- c(curr$Ne1, curr$Ne2) 
	m <- curr$M
	a <- c(curr$alpha1, curr$alpha2, curr$alpha3, curr$alpha4)
	# simulate data with proposed parameters
	sim <- sapply(psv_sim(par=list(T=t, B=b, Ne=ne, A=a, M=m)), unlist)
	if(plot) {
		plot(rec, sim[,1], col="mediumturquoise", main=paste("D (Pop. 1)\n", cntl$currIter),
			xlab="rec", ylab="D")
		lines(rec, obs$aux[,1], lwd=3)
		points(rec, obs$D1, pch=20, col="mediumturquoise")
		plot(rec, obs$D2, col="mediumturquoise", main=paste("D (Pop. 2)\n", cntl$currIter),
			xlab="rec", ylab="D")
		lines(rec, obs$aux[,2], lwd=3)
		points(rec, sim[,2], pch=20, col="mediumturquoise")
		plot(rec, obs$A1, col="mediumturquoise", main=paste("ALDER (Pop. 1)\n", cntl$currIter), 
			xlab="rec", ylab="ALDER")
		lines(rec, obs$aux[,3], lwd=3)
		points(rec, sim[,3], pch=20, col="mediumturquoise")
		plot(rec, obs$A2, col="mediumturquoise", main=paste("ALDER (Pop. 2)\n", cntl$currIter),
			xlab="rec", ylab="ALDER")
		lines(rec, obs$aux[,4], lwd=3)
		points(rec, sim[,4], pch=20, col="mediumturquoise")
	}
	# fit curve and calculate likelihood based on multivariate normal
	xhat <- colSums(sim	- obs$aux)
	s2 <- apply(sim	- obs$aux, 2, var)
	#curr$loglik <- dmvnorm(x=mu, sigma=obs$sigma, log=T)
	curr$loglik <- loglik_sample(m=xhat, s2=s2, mu=rep(0, 4), sigma2=diag(obs$sigma), 
		N=length(rec))
	#	calculate appropriate transition densities
	mod <- cntl$currIter %% cntl$P
	if(mod == 0) {
		q01 <- log(dtnorm(x=curr$M, mean=prev$M, sd=cntl$sigma.M, a=0, b=1))
		q10 <- log(dtnorm(x=prev$M, mean=curr$M, sd=cntl$sigma.M, a=0, b=1))
	} else if(mod == 1) {
		q01 <- log(dtnorm(x=curr$Ne1, mean=prev$Ne1, sd=cntl$sigma.Ne, a=1e2, b=1e5))
		q10 <- log(dtnorm(x=prev$Ne1, mean=curr$Ne1, sd=cntl$sigma.Ne, a=1e2, b=1e5))	
	} else if(mod == 2) {	
		q01 <- log(dtnorm(x=curr$Ne2, mean=prev$Ne2, sd=cntl$sigma.Ne, a=1e2, b=1e5))
		q10 <- log(dtnorm(x=prev$Ne2, mean=curr$Ne2, sd=cntl$sigma.Ne, a=1e2, b=1e5))	
	} else if(mod == 3) {	
		q01 <- log(dtnorm(x=curr$Tdiv, mean=prev$Tdiv, sd=cntl$sigma.Tdiv, a=1e1, b=Inf))
		q10 <- log(dtnorm(x=prev$Tdiv, mean=curr$Tdiv, sd=cntl$sigma.Tdiv, a=1e1, b=Inf))	
	} else if(mod == 4) {	
		q01 <- log(dtnorm(x=curr$Tmig, mean=prev$Tmig, sd=cntl$sigma.Tmig, a=1, b=prev$Tmig-1))
		q10 <- log(dtnorm(x=prev$Tmig, mean=curr$Tmig, sd=cntl$sigma.Tmig, a=1, b=curr$Tmig-1))	
	} else {
		alpha.prev <- c(prev$alpha1, prev$alpha2, prev$alpha3, prev$alpha4)
		alpha.curr <- c(curr$alpha1, curr$alpha2, curr$alpha3, curr$alpha4)
		q01 <- log(dtnorm(x=alpha.curr, mean=alpha.prev, sd=cntl$sigma.alpha, a=0, b=Inf))
		q10 <- log(dtnorm(x=alpha.prev, mean=alpha.curr, sd=cntl$sigma.alpha, a=0, b=Inf))
	}
	if(init) return(curr)
	alpha <- min(0, curr$loglik + q10 - prev$loglik - q01)
	if(log(runif(1)) < alpha) {
		return(list(acc=TRUE, param=curr))
	} else {
		return(list(acc=FALSE, param=prev))
	}
}

# define control for control MCMC	
cntl <- list() # proposal, log-likelihood, and other deets for MCMC
cntl$numIters <- 6e4L
cntl$currIter <- 1
cntl$B <- 200L 
cntl$P <- 6 # number of parameters being sampled
cntl$sigma.M <- 0.01 # variance for migration proposal
cntl$sigma.Ne <- 4e3 # variance for migration proposal
cntl$sigma2.M <- cntl$sigma.M^2 # variance for migration proposal
cntl$sigma2.Ne <- cntl$sigma.Ne^2 # variance for migration proposal
cntl$sigma.Tmig <- 50 # variance for time of migration proposal
cntl$sigma.Tdiv <- 50 # variance for time of divergence proposal
cntl$sigma2.Tmig <- cntl$sigma.Tmig^2 # variance for time of migration proposal
cntl$sigma2.Tdiv <- cntl$sigma.Tdiv^2 # variance for time of divergence proposal
cntl$sigma.alpha <- 0.5 # variance for migration proposal
cntl$sigma2.alpha <- cntl$sigma.alpha^2 # variance for migration proposal
cntl$acc <- numeric(cntl$numIters-1) # store acceptance rate

# define parameters for posterior sampling 
param <- data.frame(
	Ne1=rep(NA, cntl$numIters), 
	Ne2=rep(NA, cntl$numIters), 
	M=rep(NA, cntl$numIters), 
	alpha1=rep(NA, cntl$numIters), 
	alpha2=rep(NA, cntl$numIters), 
	alpha3=rep(NA, cntl$numIters), 
	alpha4=rep(NA, cntl$numIters), 
	Tmig=rep(NA, cntl$numIters), 
	Tdiv=rep(NA, cntl$numIters),
	loglik=rep(NA, cntl$numIters)
)

# initialize MCMC chain
param$Ne1[1] <- 5e3L
param$Ne2[1] <- 5e3L
param$M <- 0
param$alpha1[1] <- param$alpha4[1] <- 1
param$alpha2[1] <- param$alpha3[1] <- 1
param$Tmig[1] <- 500L 
param$Tdiv[1] <- 1000L

# generate "pseudo-observed dataset" and fit auxillary model
# alternatively read in data set (future)
rec <- seq(1e-5, 1e-2, length=30)
pod <- psv_sim(par=list(T=c(300L, 400L), B=cntl$B, Ne=c(1e4L, 3e3L), 
	A=c(0.05, 0.02, 0.02, 0.05), M=0.01))
pod <- lapply(pod, unlist)
aux <- lapply(pod, auxillary, x=rec)
pod$aux <- simplify2array(lapply(aux, "[[", "Y"))
pod$sigma <- cov(simplify2array(lapply(aux, "[[", "res")))
#load("pod.RData")

# calculate log-likelihood of initial state of chain
param[1,] <- wf_mh_step(obs=pod, curr=param[1,], prev=param[1,], cntl=cntl, init=T)

# run chain
par(mfrow=c(2, 2))
for(iter in 2:cntl$numIters) {
	cntl$currIter <- iter
	repeat {
		PLOT <- ifelse(iter %% (cntl$P*10) == 0, TRUE, FALSE)
		proposal <- wf_mh_propose(param[iter-1,], cntl)
		mh <- try(wf_mh_step(obs=pod, curr=proposal, prev=param[iter-1,], cntl=cntl, plot=PLOT))
		if(class(mh) != "try-error") {
			break
		} else {
			print(paste("ERROR:", paste(proposal, collapse=" ")))
			error <- proposal
		}
	}	
	if(mh$acc) {
		param[iter,] <- mh$param
	} else {
		param[iter,] <- param[iter-1,] 
	}
	cntl$acc[iter-1] <- mh$acc
	if(iter %% (cntl$P*10) == 0) {
		print(paste(iter, ":", paste(round(proposal[-10],4), collapse=" "), "|", 
			paste(round(param[iter,],4), collapse=" ")))
	}
}
