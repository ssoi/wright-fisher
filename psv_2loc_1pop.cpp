/*****************************************************************
 *	Approximate Wright-Fisher model using pseudosampling variable
 *	approach introduced by Kimura
 *
 *	author: Sameer Soi
 *****************************************************************/

using namespace Rcpp ;
using std::cout ;
using std::endl ;
using std::vector ;

unsigned int  i = 0, j = 0, k = 0, l = 0, 
b = as<unsigned int>(B), // # of 2-locus pairs
t = as<unsigned int>(T), // # of generations
n = as<unsigned int>(N) ; // effective pop. size
double p = 0.0, q = 0.0, d = 0.0, 
	maxU=sqrt(3.0), minU=-sqrt(3.0),
	tot, tmp1, tmp2, tmp3, 
	r = as<double>(R), // recombination rate
	*g = new double[4],  // allele frequencies
	*alpha = new double[4], // Dirichlet parameters
	*u = new double[3], // uniform variables
	*eta = new double[3], // noise for PSV
	*sd = new double[4], // desired standard deviation
	*c2 = new double[2],
	*c3 = new double[3] ;
vector<double> a = as< vector<double> >(A) ; 
vector< vector<double> > res ; // 2d vector returned to R


// Set up GSL RNG
gsl_rng_env_setup() ;
gsl_rng * rng ;
const gsl_rng_type * rngType = gsl_rng_ranlxd1 ;
rng = gsl_rng_alloc(rngType) ;

// Begin simulation
for(i = 0 ; i < 4 ; i++) alpha[i] = a[i] ;
res.resize(b) ;
for(i = 0; i < b ; i++) {
	res[i].resize(4) ;
	gsl_ran_dirichlet(rng, 4, alpha, g) ;
	for(j = 0; j < t; j++) {
		d = (g[0] * g[3]) - (g[1] * g[2]) ; // LD
		if(isnan(d)) {
			cout << i << " " ;
			for(l = 0 ; l < 4 ; l++) cout << g[l] << " " ;
			cout << "| " ;
			for(l = 0 ; l < 4 ; l++) cout << res[i-1][l] << " " ;
			cout << endl ;
		}
		g[0] -= r*d ;
		g[1] += r*d ;
		g[2] += r*d ;
		g[3] -= r*d ;
		for(k = 0 ; k < 3 ; k++) g[k] = g[k] > 0 ? g[k] : 0 ;
		for(k = 0 ; k < 3 ; k++) u[k] = gsl_ran_flat(rng, minU, maxU) ;
		for(k = 0 ; k < 4 ; k++) sd[k] = sqrt(g[k]*(1-g[k])/((double) n)) ;
		tmp1 = g[0]*g[1]/((1-g[0])*(1-g[1])) ;
		c2[0] = -sqrt(tmp1) ;
		c2[1] = sqrt(1-tmp1) ;
		tmp1 = g[0]*g[2]/((1-g[0])*(1-g[2])) ;
		tmp2 = g[1]*g[2]/((1-g[1])*(1-g[2])) ;
		tmp3 = 1-g[0] ;
		c3[0] = -sqrt(tmp1) ;
		c3[1] = -sqrt(tmp2)/tmp3 ;  
		c3[2] = sqrt(1 - tmp1 - tmp2/(tmp3*tmp3)) ;
		eta[0] = sd[0] * u[0] ;
		eta[1] = sd[1] * (c2[0] * u[0] + c2[1]*u[1]) ;
		eta[2] = sd[2] * (c3[0] * u[0] + c3[1]*u[1] + c3[2]*u[2]) ;
		for(tot=0.0, k = 0 ; k < 3 ; k++) {
			g[k] += eta[k] ;
			tot += g[k] ;
		}
		g[3] = 1 - tot ;
	}
	// load results 
	for(j = 0 ; j < 4 ; j++) res[i][j] = g[j] ;
	// g[j] = buff[j] ;	
}

// clean up memory
delete[] g ;
delete[] u ;
delete[] sd ;
delete[] alpha ;
delete[] c2 ; 
delete[] c3 ; 
delete[] eta ; 
gsl_rng_free(rng) ;

// return simulation results
return(wrap(res)) ;
