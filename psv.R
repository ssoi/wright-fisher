library(mvtnorm)

wf <- function(N, p, t, r) {
	# simulate WF exactly
	for(i in 2:t) {
		D <- (p[1]*p[4] - p[2]*p[3])
		p[1] <- p[1] - r*D
		p[2] <- p[2] + r*D
		p[3] <- p[3] + r*D
		p[4] <- p[4] - r*D
		p <- as.numeric(rmultinom(n=1, size=N, prob=p)/N)
	}
	return(as.numeric(p))
}

wf_mig <- function(N, p0, t, r, tm, m) {
	# simulate WF for 2 populations with migration,
	# which starts after tm generations and is from pop 2->1 
	p2 <- p1 <- p1
	for(i in 2:t) {
		if(i <= tm) p1 <- m*p2 + (1-m)*(p1)
		D <- (p1[1]*p1[4] - p1[2]*p1[3])
		p1[1] <- p1[1] - r*D
		p1[2] <- p1[2] + r*D
		p1[3] <- p1[3] + r*D
		p1[4] <- p1[4] - r*D
		p2[1] <- p2[1] - r*D
		p2[2] <- p2[2] + r*D
		p2[3] <- p2[3] + r*D
		p2[4] <- p2[4] - r*D
		if(any(p1 < 0)) p1 <- pmax(0, p1)/sum(pmax(0, p1))
		if(any(p2 < 0)) p2 <- pmax(0, p2)/sum(pmax(0, p2))
		p1 <- as.numeric(rmultinom(n=1, size=N, prob=p1)/N)
		p2 <- as.numeric(rmultinom(n=1, size=N, prob=p2)/N)
	}
	return(c(p1, p2))
}

diffuse <- function(N, p, t, r) {
	# simulate WF diffusion using Gaussian
	# with same mean and variance as multinomial probability
	res <- matrix(nrow=t, ncol=4)
	for(i in 1:t) {
		D <- (p[1]*p[4] - p[2]*p[3])
		p[1] <- p[1] - r*D
		p[2] <- p[2] + r*D
		p[3] <- p[3] + r*D
		p[4] <- p[4] - r*D
		S <- outer(p, p)/N
		diag(S) <- p*(1-p)/N
		p <- p + as.numeric(rmvnorm(n=1, sigma=S))
		p <- res[i,] <- p/sum(p)
	}
	return(p)
}

psv <- function(N, p, t, r) {
	# approximate WF diffusion using U() variables 
	# with same mean and variance as Gaussian noise 
	a <- sqrt(3)
	u <- matrix(nrow=t, ncol=3, runif(n=t*3, min=-a, max=a))
	for(i in 1:t) {
		D <- (p[1]*p[4] - p[2]*p[3])
		p[1] <- p[1] - r*D
		p[2] <- p[2] + r*D
		p[3] <- p[3] + r*D
		p[4] <- p[4] - r*D
		s <- sqrt(p*(1-p)/N)
		c21 <- -sqrt(p[1]*p[2]/((1-p[1])*(1-p[2])))
		c22 <- sqrt(1-c21^2)
		c31 <- -sqrt(p[1]*p[3]/((1-p[1])*(1-p[3])))
		c32 <- -sqrt(p[2]*p[3]/((1-p[2])*(1-p[3])))/(1-p[1])
		c33 <- sqrt(1 - c31^2 - c32^2)
		e1 <- s[1]*u[i-1,1]
		e2 <- s[2]*(c21*u[i-1,1] + c22*u[i-1,2])
		e3 <- s[3]*(c31*u[i-1,1] + c32*u[i-1,2] + c33*u[i-1,3])
		p <- p + c(e1, e2, e3, 0)
		p[4] <- 1 - sum(p[1:3])
	}
	return(p)
}

psv2 <- function(N, p, t, r, i) {
	# accelerate approximation using Euler method
	# sum drift within interval [T,T+t/i)
	# diffusion within interval [T,T+t/i) is proportial to sqrt(t/i)
	delta <- t/i
	a <- sqrt(3)
	b <- sqrt(delta)
	u <- matrix(nrow=t-1, ncol=3, runif(n=(t-1)*3, min=-a, max=a))
	curr <- 0
	while(curr < t) {
		D <- (p[1]*p[4] - p[2]*p[3])
		p[1] <- p[1] - r*D*delta
		p[2] <- p[2] + r*D*delta
		p[3] <- p[3] + r*D*delta
		p[4] <- p[4] - r*D*delta
		s <- b*sqrt(p*(1-p)/N)
		c21 <- -sqrt(p[1]*p[2]/((1-p[1])*(1-p[2])))
		c22 <- sqrt(1-c21^2)
		c31 <- -sqrt(p[1]*p[3]/((1-p[1])*(1-p[3])))
		c32 <- -sqrt(p[2]*p[3]/((1-p[2])*(1-p[3])))/(1-p[1])
		c33 <- sqrt(1 - c31^2 - c32^2)
		e1 <- s[1]*u[i-1,1]
		e2 <- s[2]*(c21*u[i-1,1] + c22*u[i-1,2])
		e3 <- s[3]*(c31*u[i-1,1] + c32*u[i-1,2] + c33*u[i-1,3])
		p[1] <- max(min(1, p[1] + e1), 0)
		p[2] <- max(min(1, p[2] + e2), 0)
		p[3] <- max(min(1, p[3] + e3), 0)
		p[4] <- 1 - sum(p[-4])
		if(any(p == 0) | any(p == 1)) return(p)
		curr <- curr + delta
		print(paste(curr, paste(p, collapse=" "), collapse=" "))
	}
	return(p)
}

# simulate a given scenario under all 3 schemes
sim.wf <- replicate(1000, wf(N=1e4, p=c(0.4,0.1,0.1,0.4), t=100, r=0))
sim.diffuse <- replicate(1000, diffuse(N=1e4, p=c(0.4,0.1,0.1,0.4), t=100, r=0))
sim.psv <- replicate(1000, psv(N=1e4, p=c(0.4,0.1,0.1,0.4), t=100, r=0))

# calculate simulated statistic to evaluate similarity w.r.t to WF (reference)
D.wf <- apply(approx.wf, 2, function(x) x[1]*x[4] - x[2]*x[3])
D.diffuse <- apply(sim.diffuse, 2, function(x) x[1]*x[4] - x[2]*x[3])
D.psv <- apply(sim.psv, 2, function(x) x[1]*x[4] - x[2]*x[3])

# visualize simulations of D
plot(density(D.diffuse, bw=0.01), lwd=2, lty=2, 
	main="D simulated using different schemes")
lines(density(D.wf, bw=0.01), lwd=2, col="indianred2")
lines(density(D.psv, bw=0.01), lwd=2, lty=3)
legend("topright", legend=c("WF", "Gaussian", "PSV"), 	
	col=c("indianred2", "black", "black"), lwd=2, lty=1:3, box.col=NA)

# test difference between D simulated under PSV/Gaussian and WF
t.test(D.wf, D.diffuse)
# t = -2.4333, df = 1903.515, p-value = 0.01506
# alternative hypothesis: true difference in means is not equal to 0
t.test(D.wf, D.psv)
# t = -0.8436, df = 1992.616, p-value = 0.399
# alternative hypothesis: true difference in means is not equal to 0
