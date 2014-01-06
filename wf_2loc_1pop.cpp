/************************************************************
*	Simulate Wright-Fisher population using
*	multinomial RVs for B pairs of bi-allelic
*	loci with recombination rate r from gen. 0
*	to gen. T. Initial configuration of haplotypes G
*	with constant effective pop. size N, mutation
*	rate U (no 1->0 mutations).
*
*	author: Sameer Soi
************************************************************/

using namespace Rcpp ;
using std::vector ;

unsigned int  i = 0, j = 0, k = 0,
	b = as<unsigned int>(B), // # of 2-locus pairs
	t = as<unsigned int>(T), // # of generations
	*H = new unsigned int[4] ;
double p = 0.0, q = 0.0, d = 0.0, 
	n = as<double>(N), // effective pop. size
	r = as<double>(R), // recombination rate
	//u = as<double>(U), // mutation rate
	*g = new double[4], 
	*prob = new double[4], 
	*alpha = new double[4] ; 
vector<double> a = as< vector<double> >(A) ; // Dirichlet parameters
vector< vector<double> > res ; // 2d vector returned to R

// Set up GSL RNG
gsl_rng_env_setup() ;
gsl_rng * rng ;
const gsl_rng_type * rngType ;
rngType = gsl_rng_ranlxd1 ;
rng = gsl_rng_alloc(rngType) ;

// Begin simulation
for(i = 0 ; i < 4 ; i++) alpha[i] = a[i] ;
res.resize(b) ;
for(i = 0; i < b ; i++) {
	// res[i].resize(4) ;
	res[i].resize(t) ;
	gsl_ran_dirichlet(rng, 4, alpha, g) ;
	for(j = 0; j < t; j++) {
		res[i][j] = d = (g[0] * g[3]) - (g[1] * g[2]) ; // LD
		// mutation and recombination
		prob[0] = g[0] - r*d ;
		prob[1] = g[1] + r*d ;
		prob[2] = g[2] + r*d ;
		prob[3] = g[3] - r*d ;
		// G(t+1) | G(t) ~ Multi(N,p0(t)-rD,p1(t)+rD,p2+rD,p3-rD)
		gsl_ran_multinomial(rng, 4, n, prob, H) ;
		for(k = 0 ; k < 4 ; k++) g[k] = ((double)H[k])/n ;
	}
	// load results 
	// for(j = 0 ; j < 4 ; j++) res[i][j] = g[j] ;
}

// clean up memory
delete[] H ;
delete[] prob ; 
gsl_rng_free(rng) ;

// return simulation results
return(wrap(res)) ;
