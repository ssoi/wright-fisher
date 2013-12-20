/************************************************************
*	Approximate Wright-Fisher model using pseudosampling variable
*	approach introduced by Kimura
*
*	author: Sameer Soi
************************************************************/
#include <Rcpp.h>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace Rcpp ;
using std::vector ;

void gen_psv_var(double *p, unsigned int ne, double *u) {
	unsigned int i, j ;
	double tot, tmp1, tmp2, tmp3, *psv = new double[4],  *eta = new double[3], 
		*sd = new double[4], *c2 = new double[2], *c3 = new double[3] ;

	for(i = 0 ; i < 4 ; i++) sd[i] = sqrt(p[i]*(1-p[i])/((double) ne)) ;
	tmp1 = p[0]*p[1]/((1-p[0])*(1-p[1])) ;
	c2[0] = -sqrt(tmp1) ;
	c2[1] = sqrt(1-tmp1) ;
	tmp1 = p[0]*p[2]/((1-p[0])*(1-p[2])) ;
	tmp2 = p[1]*p[2]/((1-p[1])*(1-p[2])) ;
	tmp3 = 1-p[0] ;
	c3[0] = -sqrt(tmp1) ;
	c3[1] = -sqrt(tmp2)/tmp3 ;  
	c3[2] = sqrt(1 - tmp1 - tmp2/(tmp3*tmp3)) ;
	eta[0] = sd[0] * u[0] ;
	eta[1] = sd[1] * (c2[0] * u[0] + c2[1]*u[1]) ;
	eta[2] = sd[2] * (c3[0] * u[0] + c3[1]*u[1] + c3[2]*u[2]) ;
	for(tot=0.0, i = 0 ; i < 3 ; i++) {
		p[i] += eta[i] ;
		tot += p[i] ;
	}
	p[3] = 1 - tot ;
	return ;	
}

// [[Rcpp::export]]
NumericMatrix psv_2loc_1pop(NumericVector N, NumericVector T, NumericVector A, 
	NumericVector R, NumericVector B) {
	unsigned int  i = 0, j = 0, k = 0,
		b = as<unsigned int>(B), // # of 2-locus pairs
		t = as<unsigned int>(T), // # of generations
		n = as<unsigned int>(N), // effective pop. size
		*H = new unsigned int[4] ;
	double p = 0.0, q = 0.0, d = 0.0, maxU=sqrt(3.0), minU=-maxU,
		r = as<double>(R), // recombination rate
		*u = new double[3], // uniform variables
		*g = new double[4], 
		*prob = new double[4], 
		*alpha = new double[4], 
		*buff = new double[4] ;
	vector<double> a = as< vector<double> >(A) ; // Dirichlet parameters
	vector< vector<double> > res ; // 2d vector returned to R
	
	// Set up GSL RNG
	gsl_rng_env_setup() ;
	gsl_rng * rng ;
	const gsl_rng_type * rngType = gsl_rng_default ;
	rng = gsl_rng_alloc(rngType) ;
	
	// Begin simulation
	for(i = 0 ; i < 4 ; i++) alpha[i] = a[i] ;
	res.resize(b) ;
	for(i = 0; i < b ; i++) {
		res[i].resize(4) ;
		gsl_ran_dirichlet(rng, 4, alpha, g) ;
		for(j = 0 ; j < 4 ; j++) buff[j] = g[j] ;
		for(j = 0; j < t; j++) {
			d = (g[0] * g[3]) - (g[1] * g[2]) ; // LD
			g[0] -= r*d ;
			g[1] += r*d ;
			g[2] += r*d ;
			g[3] -= r*d ;
			for(k = 0 ; k < 3 ; k++) u[k] = gsl_ran_flat(rng, minU, maxU) ;
			gen_psv_var(g, n, u) ;
		}
		// load results and re-initialize buffer
		for(j = 0 ; j < 4 ; j++) {
			res[i][j] = g[j] ;
			g[j] = buff[j] ;
		}
	}
	
	// clean up memory
	delete[] H ;
	delete[] buff ;
	delete[] prob ; 
	gsl_rng_free(rng) ;
	
	// return simulation results
	return(wrap(res)) ;
}
