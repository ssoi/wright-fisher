/************************************************************
*	Simulate Wright-Fisher population using
*	multinomial RVs for B pairs of bi-allelic
*	loci with recombination rate r from gen. 0
*	to gen. T. Initial configuration of haplotypes G
*	with constant effective pop. size N
*
*	author: Sameer Soi
************************************************************/

using namespace Rcpp ;
using std::vector ;

unsigned int i = 0, j = 0, k = 0, l = 0,
	b = as<unsigned int>(B), // # of 2-locus pairs
	*H1 = new unsigned int[4],
	*H2 = new unsigned int[4] ;
double p = 0.0, q = 0.0, d1 = 0.0, d2 = 0.0,
	g1Sum, g2Sum,
	r = as<double>(R), // recombination rate
	*ne = new double[2], 
	*mig = new double[2], 
	*tmp = new double[4], 
	*g1 = new double[4], 
	*g2 = new double[4], 
	*prob1 = new double[4], 
	*prob2 = new double[4], 
	*alpha = new double[4], 
	*buff = new double[4] ;
vector<unsigned int> n = as< vector<unsigned int> >(N) ; // effective pop. sizes
vector< vector<double> > res ; // 2d vector returned to R
vector<double> a = as< vector<double> >(A) ; // Dirichlet parameters
vector<double> m = as< vector<double> >(M) ; // migration rates
vector<double> t = as< vector<double> >(T) ; // # of generations

// Set up GSL RNG
gsl_rng_env_setup() ;
gsl_rng * rng ;
const gsl_rng_type * rngType ;
rngType = gsl_rng_ranlxd1 ;
rng = gsl_rng_alloc(rngType) ;

// Begin simulation
ne[0] = (double) n[0] ;
mig[0] = mig[1] = 0.0 ;
for(i = 0 ; i < 4 ; i++) {
	prob1[i] = prob2[i] = g1[i] = g2[i] = 0 ;
	alpha[i] = a[i] ;
}
gsl_ran_dirichlet(rng, 4, alpha, g1) ;
for(i = 0 ; i < 4 ; i++) buff[i] = g1[i] ;
res.resize(b) ;
for(i = 0; i < b ; i++) {
	res[i].resize(8) ;
	for(j = 0; j < t[2]; j++) {
		if(j == t[0]) { // change population sizes after divergence
			ne[0] = (double) n[1] ;
			ne[1] = (double) n[2] ;
			for(k = 0 ; k < 4 ; k++) g2[k] = g1[k] ; // initialize second population
		}
		if(j == t[1]) { // set migration sizes after divergence
			mig[0] = m[0] ;
			mig[1] = m[1] ;
		}
		d1 = (g1[0] * g1[3]) - (g1[1] * g1[2]) ; // LD
		prob1[0] = (g1[0] + g2[0]*mig[0]) - r*d1 ;
		prob1[1] = (g1[1] + g2[1]*mig[0]) + r*d1 ;
		prob1[2] = (g1[2] + g2[2]*mig[0]) + r*d1 ;
		prob1[3] = (g1[3] + g2[3]*mig[0]) - r*d1 ;
		for(g1Sum = 0, k = 0 ; k < 4 ; k++) g1Sum += g1[k] ;
		for(k = 0 ; k < 4 ; k++) g1[k] = g1[k]/g1Sum ;
		for(k = 0 ; k < 4 ; k++) tmp[k] = g1[k] ;
		gsl_ran_multinomial(rng, 4, (unsigned int) ne[0], prob1, H1) ;
		for(k = 0 ; k < 4 ; k++) g1[k] = ((double)H1[k])/ne[0] ;
		if(j >= t[0]) {
			d2 = (g2[0] * g2[3]) - (g2[1] * g2[2]) ; // LD
			prob2[0] = (g2[0] + g1[0]*mig[1]) - r*d2 ;
			prob2[1] = (g2[1] + g1[1]*mig[1]) + r*d2 ;
			prob2[2] = (g2[2] + g1[2]*mig[1]) + r*d2 ;
			prob2[3] = (g2[3] + g1[3]*mig[1]) - r*d2 ;
			for(g2Sum = 0, k = 0 ; k < 4 ; k++) g2Sum += g2[k] ;
			for(k = 0 ; k < 4 ; k++) g2[k] = g2[k]/g2Sum ;
			gsl_ran_multinomial(rng, 4, (unsigned int) ne[1], prob2, H2) ;
			for(k = 0 ; k < 4 ; k++) g2[k] = ((double) H2[k])/ne[1] ;
		}
	}
	// load results and re-initialize buffer
	for(j = 0 ; j < 4 ; j++) {
		res[i][j] = g1[j] ;
		res[i][j+4] = g2[j] ;
		g1[j] = buff[j] ;
	}
}

// clean up memory
delete[] g1 ;
delete[] g2 ;
delete[] H1 ;
delete[] H2 ;
delete[] ne ;
delete[] mig ;
delete[] buff ;
delete[] alpha ;
delete[] prob1 ; 
delete[] prob2 ; 
gsl_rng_free(rng) ;

// return simulation results
return(wrap(res)) ;
