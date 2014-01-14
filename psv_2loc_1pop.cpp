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

bool DOMULTI = false ;
unsigned int  i = 0, j = 0, k = 0, l = 0, 
b = as<unsigned int>(B), // # of 2-locus pairs
t = as<unsigned int>(T), // # of generations
n = as<unsigned int>(N),
*H = new unsigned int[4] ; // effective pop. size
double p = 0.0, q = 0.0, d = 0.0, dn = (double) n, 
	maxU=sqrt(3.0), minU=-sqrt(3.0),
	tot, tmp1, rho12, rho13, rho23, 
	c21, c22, c31, c32, c33, 
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
		g[0] -= r*d ;
		g[1] += r*d ;
		g[2] += r*d ;
		g[3] -= r*d ;
		for(k = 0 ; k < 3 ; k++) u[k] = gsl_ran_flat(rng, minU, maxU) ;
		for(k = 0 ; k < 4 ; k++) sd[k] = sqrt(g[k]*(1-g[k])/dn) ;
		rho12 = (g[0]*g[1])/((1-g[0])*(1-g[1])) ;
		rho13 = (g[0]*g[2])/((1-g[0])*(1-g[2])) ;
		rho23 = (g[1]*g[2])/((1-g[1])*(1-g[2])) ;
		c21 = -sqrt(rho12) ;
		c22 = sqrt(1 - rho12) ;
		c31 = -sqrt(rho13) ;
		c32 = (-sqrt(rho23) - sqrt(rho12*rho13))/c22 ;
		c33 = sqrt(1 - (c31*c31) - (c32*c32)) ;
		eta[0] = sd[0] * u[0] ;
		eta[1] = sd[1] * (c21 * u[0] + c22*u[1]) ;
		eta[2] = sd[2] * (c31 * u[0] + c32*u[1] + c33*u[2]) ;
		for(tot=0.0, k = 0 ; k < 3 ; k++) {
			tot += (g[k]+eta[k]) ;
			if(g[k] + eta[k] < 0) DOMULTI = true ;
		}
		if(1 - tot < 0 ) DOMULTI = true ;
		// if((isnan(c21) || isnan(c22) || isnan(c31) || isnan(c32) || isnan(c33)) || 
		//	 isnan(g[0]) || isnan(g[1]) || isnan(g[2]) || isnan(g[3])) {
		// }
		if(!DOMULTI) {
			for(k = 0 ; k < 3 ; k++) g[k] += eta[k] ;
			g[3] = 1 - tot ;
		} else {
			gsl_ran_multinomial(rng, 4, dn, g, H) ;
			for(k = 0 ; k < 3 ; k++) g[k] = H[k]/((double) n) ;
		}
		DOMULTI = false ;
	}
	// load results 
	for(j = 0 ; j < 4 ; j++) res[i][j] = g[j] ;
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
