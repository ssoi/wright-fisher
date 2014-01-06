library(Rcpp)
library(inline)

# load C++ source code for WF simulation of 2 loci in 1 pop
wf_2loc_1pop_src <- paste(readLines("wf_2loc_1pop.cpp"), collapse="\n")
wf_2loc_1pop_sig <- signature(T="integer", B="integer", N="integer",
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

# load C++ source code for PSV simulation of 2 loci in 1 pop
psv_2loc_1pop_src <- paste(readLines("psv_2loc_1pop.cpp"), collapse="\n")
psv_2loc_1pop_sig <- signature(T="integer", B="integer", N="integer",
	R="numeric", A="numeric")
psv_2loc_1pop <- cfunction(psv_2loc_1pop_sig, body=psv_2loc_1pop_src, Rcpp=TRUE, 
	includes=c("#include <gsl/gsl_randist.h>", "#include <gsl/gsl_rng.h>"), 
	libargs="-lgsl -lgslcblas")

s <- psv_2loc_1pop(T=300L, B=50L, N=1e4L, R=1.0e-3, A=c(0.5, 0.1, 0.1, 0.5))

# load C++ source code for PSV simulation of 2 loci in 1 pop
psv_2loc_2pop_src <- paste(readLines("psv_2loc_2pop.cpp"), collapse="\n")
psv_2loc_2pop_sig <- signature(T="integer", B="integer", N="integer", R="numeric", 
	A="numeric",  M="numeric")
psv_2loc_2pop <- cfunction(psv_2loc_2pop_sig, body=psv_2loc_2pop_src, Rcpp=TRUE, 
	includes=c("#include <gsl/gsl_randist.h>", "#include <gsl/gsl_rng.h>"), 
	libargs="-lgsl -lgslcblas")

wf_sim <- function(par) {
	par$T <- as.integer(round(par$T))
	par$Ne <- as.integer(round(par$Ne))
	sim <- sapply(rec, function(r) {
		s <- wf_2loc_2pop(T=par$T, B=par$B, N=par$Ne, R=r, A=par$A, M=par$M)
		mat <- matrix(ncol=8, nrow=length(s), data=unlist(s), byrow=TRUE)
		d1 <- mat[,1]*mat[,4] - mat[,2]*mat[,3] # LD in pop 1
		d2 <- mat[,5]*mat[,8] - mat[,6]*mat[,7] # LD in pop 2
		delta <- (mat[,1]+mat[,2]-mat[,5]-mat[,6]) * 
			(mat[,1]+mat[,3]-mat[,5]-mat[,7]) 
		a1 <- d1*delta
		a2 <- d2*delta
		colMeans(cbind(d1, d2, a1, a2))
	})
	return(list(D1=sim[1,], D2=sim[2,], A1=sim[3,], A2=sim[4,]))
}


psv_sim <- function(par) {
	par$T <- as.integer(round(par$T))
	par$Ne <- as.integer(round(par$Ne))
	mat <- mclapply(rec, function(r)
		simplify2array(psv_2loc_2pop(T=par$T, B=par$B, N=par$N, R=r, A=par$A,
			M=par$M)))
	d1 <- mclapply(mat, function(r) 
		(apply(r, 2, function(x) x[1]*x[4]-x[2]*x[3])))
	d2 <- mclapply(mat, function(r) 
		mean(apply(r, 2, function(x) x[5]*x[8]-x[6]*x[7])))
	a1 <- mclapply(mat, function(r) {
		d <- apply(r, 2, function(x) x[1]*x[4]-x[2]*x[3])
		f1 <- apply(r, 2, function(x) x[1]+x[2]-x[5]-x[6])
		f2 <- apply(r, 2, function(x) x[1]+x[3]-x[5]-x[7])
		d*f1*f2
	})
	a2 <- mclapply(mat, function(r) {
		d <- apply(r, 2, function(x) x[5]*x[8]-x[6]*x[7])
		f1 <- apply(r, 2, function(x) x[1]+x[2]-x[5]-x[6])
		f2 <- apply(r, 2, function(x) x[1]+x[3]-x[5]-x[7])
		mean(d*f1*f2)
	})
	return(list(D1=d1, D2=d2, A1=a1, A2=a2))
}
