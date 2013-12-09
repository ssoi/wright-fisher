library(inline)
library(hexbin)
library(multicore)

# load C++ source code for WF simulation of 2 loci in 1 pop
wf_2loc_1pop_src <- paste(readLines("wf_2loc_1pop.cpp"), collapse="\n")
wf_2loc_1pop_sig <- signature(T="integer", B="integer", N="numeric",
	R="numeric", A="numeric")
wf_2loc_1pop <- cfunction(wf_2loc_1pop_sig, body=wf_2loc_1pop_src, Rcpp=TRUE, 
	includes=c("#include <gsl/gsl_randist.h>", "#include <gsl/gsl_rng.h>"), 
	libargs="-lgsl -lgslcblas")

# load C++ source code for WF simulation of 2 loci in 2 pop w/ migration
wf_2loc_2pop_src <- paste(readLines("wf_2loc_2pop.cpp"), collapse="\n")
wf_2loc_2pop_sig <- signature(T="integer", B="integer", N="numeric", R="numeric", 
	A="numeric",  M="numeric")
wf_2loc_2pop <- cfunction(wf_2loc_2pop_sig, body=wf_2loc_2pop_src, Rcpp=TRUE, 
	includes=c("#include <gsl/gsl_randist.h>", "#include <gsl/gsl_rng.h>"), 
	libargs="-lgsl -lgslcblas")

# define simulation parameters
b <- 200L # number of loci T <- 1000L # number of generations
t <- c(100L, 800L, 1000L) # number of loci T <- 1000L # number of generations
n <- c(10000, 10000, 3000) # effective population sizes
m <- c(0.05, 0.000) # migration rates
theta <- 2 
rec <- seq(1e-6, 1e-4, length=50)
res <- simplify2array(mclapply(rec, function(r) {
	s <- wf_2loc_2pop(T=t, B=b, N=n, R=r, A=rep(theta, 4), M=m)
	sim <- matrix(ncol=8, nrow=length(s), data=unlist(s), byrow=TRUE)
	d1 <- sim[,1]*sim[,4] - sim[,2]*sim[,3]
	d2 <- sim[,5]*sim[,8] - sim[,6]*sim[,7]
	a <- d1 * (sim[,1] + sim[,2] - sim[,5] - sim[,6]) * 
		(sim[,3] + sim[,4] - sim[,7] - sim[,8])
	colMeans(cbind(d1, d2, a))
}))
