library(Rcpp)
library(inline)

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

wf_sim <- function(pars) {
	sim <- sapply(rec, function(r) {
		s <- wf_2loc_2pop(T=pars$T, B=pars$B, N=pars$Ne, R=r, A=pars$A, M=pars$M)
		mat <- matrix(ncol=8, nrow=length(s), data=unlist(s), byrow=TRUE)
		d1 <- mat[,1]*mat[,4] - mat[,2]*mat[,3] # LD in pop 1
		d2 <- mat[,5]*mat[,8] - mat[,6]*mat[,7] # LD in pop 2
		a <- d1 * (mat[,5]+mat[,6]-mat[,1]-mat[,2]) * 
			(mat[,7]+mat[,8]-mat[,3]-mat[,4]) # calculate ALDER w.r.t. pop 1
		colMeans(cbind(d1, d2, a))
	})
	return(sim)
}

# generate "pseudo-observed dataset"
rec <- seq(1e-6, 1e-3, length=50)
pod <- wf_sim(pars=list(T=c(1L, 150L, 200L), B=100L, Ne=c(1e4L, 1e4L, 1e4L),
	A=c(2, 1, 1, 2), M=c(0.01, 0.0)))
