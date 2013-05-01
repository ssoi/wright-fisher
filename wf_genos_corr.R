library(inline)
library(hexbin)
library(multicore)

# load C++ source code for WF simulation
wf_2loc_1pop_src <- paste(readLines("wf_2loc_1pop.cpp"), collapse="\n")
wf_2loc_2pop_src <- paste(readLines("wf_2loc_2pop.cpp"), collapse="\n")
# define Rcpp function signatures
wf_2loc_1pop_sig <- signature(T="integer", B="integer", N="numeric", U="numeric", 
	R="numeric", G="numeric")
wf_2loc_2pop_sig <- signature(T="integer", B="integer", N="numeric", U="numeric", 
	R="numeric", G="numeric", H="numeric", M="numeric")
# compile source code
wf_2loc_1pop <- cfunction(wf_2loc_1pop_sig, body=wf_2loc_1pop_src, Rcpp=TRUE, 
	includes=c("#include <gsl/gsl_randist.h>", "#include <gsl/gsl_rng.h>"), 
	libargs="-lgsl -lgslcblas")
wf_2loc_2pop <- cfunction(wf_2loc_2pop_sig, body=wf_2loc_2pop_src, Rcpp=TRUE, 
	includes=c("#include <gsl/gsl_randist.h>", "#include <gsl/gsl_rng.h>"), 
	libargs="-lgsl -lgslcblas")
# define simulation parameters
S <- 1e5
B <- 1 # number of loci 
T <- 1e4 # number of generatnois
N <- rep(1.0e3, 2) # effective population sizes
M <- c(0.0001, 0.0001) # migration rates
theta <- 0.04
p <- rbeta(S, theta, theta)
q <- rbeta(S, theta, theta)
G <- apply(cbind(p, q), 1, function(f) 
	return(c(prod(f), (1-f[1])*f[2], f[1]*(1-f[2]), prod(1-f)))
)
fixed <- apply(G, 2, function(g) any(g==1))
s <- apply(G[,which(!fixed)], 4, function(g) 
	wf_2loc_1pop(T=as.integer(T), B=as.integer(B), N=10000.0, U=1e-8, R=0.001, G=g)
)
sim <- matrix(ncol=length(s), nrow=4, data=unlist(s))
