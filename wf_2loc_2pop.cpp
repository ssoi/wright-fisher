/************************************************************
*	Simulate two Wright-Fisher populations using
*	multinomial RVs for B pairs of bi-allelic
*	loci with recombination rate r from gen. 0
*	to gen. T. Initial configuration of haplotypes G
*	with constant effective pop. sizes N={n1, n2}, 
*	mutation rate U and migration rates M={m12, m21}.
*
*	author: Sameer Soi
************************************************************/


using namespace Rcpp ;
using std::vector ;

unsigned int  i = 0, j = 0, k = 0,
	b = as<unsigned int>(B), // # of 2-locus pairs
	t = as<unsigned int>(T), // # of generations
	*H2 = new unsigned int[4], 
	*H1 = new unsigned int[4] ;
double d1 = 0.0, d2 = 0.0, r = as<double>(R),
	u = as<double>(U),
	rho1 = new double[4],
	rho2 = new double[4],
	*prob1 = new double[4], 
	*prob2 = new double[4], 
	*buff1 = new double[4], 
	*buff2 = new double[4] ;
vector<double>	m = as< vector<double> >(M), 
	n = as< vector<double> >(N), // effective pop.sizes for 2 pops
	g = as< vector<double> >(G), // initial haplotype configuration
	h = as< vector<double> >(H) ; // initial haplotype configuration
vector< vector<double> > res ; // 2d vector returned to R

// Set up GSL RNG
gsl_rng * rng ;
const gsl_rng_type * rngType ;
gsl_rng_env_setup() ;
rngType = gsl_rng_ranlxd1 ;
rng = gsl_rng_alloc(rngType) ;

// Begin simulation
for(i = 0 ; i < 4 ; i++) buff1[i] = g[i] ;
for(i = 0 ; i < 4 ; i++) buff2[i] = h[i] ;
res.resize(b) ;
for(i = 0; i < b ; i++) {
	res[i].resize(8) ;	
	for(j = 0; j < t; j++) {
		d1 = (g[0] * g[3]) - (g[1] * g[2]) ; //  LD in pop. 1
		d2 = (h[0] * h[3]) - (h[1] * h[2]) ; //  LD in pop. 2
		// migration
		for(k = 0 ; k < 4 ; k++) {
			rho1[k] = n[1]*m[1]*h[k] ; // # migrants from pop. 2 -> 1
			rho2[k] = n[0]*m[0]*g[k] ; // # migrants from pop. 1 -> 2 
			prob1[k] = (n[0]*g[k] + rho1[k])/(n[0] + rho1[k]) ;
			prob2[k] = (n[1]*h[k] + rho2[k])/(n[1] + rho2[k]) ;
		}
		// mutation and recombination
		prob1[0] = (1-u)*(g[0] - r*d1) - 2*u ;
		prob1[1] = (1-u)*(g[1] + r*d1) + u*g[0] ;
		prob1[2] = (1-u)*(g[2] + r*d1) + u*g[0] ;
		prob1[3] = (1-u)*(g[3] - r*d1) + u*(g[1] + g[2]) ;
		prob2[0] = (1-u)*(h[0] - r*d2) - 2*u ;
		prob2[1] = (1-u)*(h[1] + r*d2) + u*h[0] ;
		prob2[2] = (1-u)*(h[2] + r*d2) + u*h[0] ;
		prob2[3] = (1-u)*(h[3] - r*d2) + u*(h[1] + h[2]) ;
		// G(t+1) | G(t) ~ Multi(N,p0(t)-rD,p1(t)+rD,p2+rD,p3-rD)
		gsl_ran_multinomial(rng, 4, n[0], prob1, H1) ;
		gsl_ran_multinomial(rng, 4, n[1], prob2, H2) ;
		for(k = 0 ; k < 4 ; k++) g[k] = ((double) H1[k])/n[0] ;
		for(k = 0 ; k < 4 ; k++) h[k] = ((double) H2[k])/n[1] ;
	}
	//  store results of simulation and reset for next round
	for(k = 0 ; k < 4 ; k++) {
		res[i][k] = g[k] ;
		res[i][k + 4] = h[k] ;
		g[k] = buff1[k] ;
		h[k] = buff2[k] ;
	}
}

// free memory before returning results
delete[] buff1 ;
delete[] buff2 ;
delete[] prob1 ; 
delete[] prob2 ; 
delete[] H1 ;
delete[] H2 ;
gsl_rng_free(rng) ;

// return results
return(wrap(res)) ;
