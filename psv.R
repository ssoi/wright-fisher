psv <- function(N, p, t, r) {
	a <- sqrt(3)
	u <- matrix(nrow=t-1, ncol=3, runif(n=(t-1)*3, min=-a, max=a))
	for(i in 2:t) {
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
		p[1] <- max(min(1, p[1] + e1), 0)
		p[2] <- max(min(1, p[2] + e2), 0)
		p[3] <- max(min(1, p[3] + e3), 0)
		p[4] <- 1 - sum(p[-4])
		if(any(p == 0) | any(p == 1)) return(p)
	}
	return(p)
}

psv2 <- function(N, p, t, r, i) {
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

wf <- function(N, p, t, r) {
	for(i in 2:t) {
		D <- (p[1]*p[4] - p[2]*p[3])
		p[1] <- p[1] - r*D
		p[2] <- p[2] + r*D
		p[3] <- p1[3] + r*D
		p[4] <- p[4] - r*D
		p1 <- rmultinom(n=1, size=N, prob=p1)/N
		p2 <- rmultinom(n=1, size=N, prob=p2)/N
	}
	return(c(p1, p2))
}

wf_mig <- function(N, p0, t, r, tm, m) {
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
	for(i in 2:t) {
		D <- (p[1]*p[4] - p[2]*p[3])
		p[1] <- p[1] - r*D
		p[2] <- p[2] + r*D
		p[3] <- p1[3] + r*D
		p[4] <- p[4] - r*D
		S <- outer(p, p)/N
		diag(S) <- p*(1-p)/N
		p <- p + as.numeric(rmvnorm(n=1, sigma=S))
	}
	return(p)
}
