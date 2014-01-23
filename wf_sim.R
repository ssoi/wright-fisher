library(Rcpp)
library(inline)

# load C++ source code for PSV simulation of 2 loci in 1 pop
psv_2loc_1pop_src <- paste(readLines("psv_2loc_1pop.cpp"), collapse="\n")
psv_2loc_1pop_sig <- signature(T="integer", B="integer", N="integer",
	R="numeric", A="numeric")
psv_2loc_1pop <- cfunction(psv_2loc_1pop_sig, body=psv_2loc_1pop_src, Rcpp=TRUE, 
	includes=c("#include <gsl/gsl_randist.h>", "#include <gsl/gsl_rng.h>"), 
	libargs="-lgsl -lgslcblas")

# load C++ source code for PSV simulation of 2 loci in 1 pop
psv_2loc_2pop_src <- paste(readLines("psv_2loc_2pop.cpp"), collapse="\n")
psv_2loc_2pop_sig <- signature(T="integer", B="integer", N="integer", R="numeric", 
	A="numeric",  M="numeric")
psv_2loc_2pop <- cfunction(psv_2loc_2pop_sig, body=psv_2loc_2pop_src, Rcpp=TRUE, 
	includes=c("#include <gsl/gsl_randist.h>", "#include <gsl/gsl_rng.h>"), 
	libargs="-lgsl -lgslcblas")

psv_sim <- function(par) {
	par$T <- as.integer(round(par$T))
	par$Ne <- as.integer(round(par$N))
	mat <- mclapply(rec, function(r)
		simplify2array(psv_2loc_2pop(T=par$T, B=par$B, N=par$N, R=r, A=par$A,
			M=par$M)))
	d1 <- mclapply(mat, function(r) 
		mean(apply(r, 2, function(x) x[1]*x[4]-x[2]*x[3]), na.rm=TRUE))
	d2 <- mclapply(mat, function(r) 
		mean(apply(r, 2, function(x) x[5]*x[8]-x[6]*x[7]), na.rm=TRUE))
	a1 <- mclapply(mat, function(r) {
		d <- apply(r, 2, function(x) x[1]*x[4]-x[2]*x[3])
		f1 <- apply(r, 2, function(x) x[1]+x[2]-x[5]-x[6])
		f2 <- apply(r, 2, function(x) x[1]+x[3]-x[5]-x[7])
		mean(d*f1*f2, na.rm=TRUE)
	})
	a2 <- mclapply(mat, function(r) {
		d <- apply(r, 2, function(x) x[5]*x[8]-x[6]*x[7])
		f1 <- apply(r, 2, function(x) x[1]+x[2]-x[5]-x[6])
		f2 <- apply(r, 2, function(x) x[1]+x[3]-x[5]-x[7])
		mean(d*f1*f2, na.rm=TRUE)
	})
	return(list(D1=d1, D2=d2, A1=a1, A2=a2))
}
