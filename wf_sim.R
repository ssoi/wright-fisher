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
wf_2loc_2pop_sig <- signature(T="integer", B="integer", N="integer", R="numeric", 
	A="numeric",  M="numeric")
wf_2loc_2pop <- cfunction(wf_2loc_2pop_sig, body=wf_2loc_2pop_src, Rcpp=TRUE, 
	includes=c("#include <gsl/gsl_randist.h>", "#include <gsl/gsl_rng.h>"), 
	libargs="-lgsl -lgslcblas")

wf_sim <- function(par) {
	par$Ne <- as.integer(round(par$Ne))
	sim <- sapply(rec, function(r) {
		s <- wf_2loc_2pop(T=par$T, B=par$B, N=par$Ne, R=r, A=par$A, 
			M=par$M)
		mat <- matrix(ncol=8, nrow=length(s), data=unlist(s), byrow=TRUE)
		d1 <- mat[,1]*mat[,4] - mat[,2]*mat[,3] # LD in pop 1
		d2 <- mat[,5]*mat[,8] - mat[,6]*mat[,7] # LD in pop 2
		delta <- (mat[,5]+mat[,6]-mat[,1]-mat[,2]) * 
			(mat[,7]+mat[,8]-mat[,3]-mat[,4]) 
		a1 <- d1*delta
		a2 <- d2*delta
		colMeans(cbind(d1, d2, a1, a2))
	})
	return(list(D1=sim[1,], D2=sim[2,], A1=sim[3,]))
}
