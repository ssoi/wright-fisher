/*****************************************************************
 *	Approximate Wright-Fisher model for 2 populations using 
 *  pseudosampling variable approach introduced by Kimura
 *
 *	author: Sameer Soi
 *****************************************************************/

using namespace Rcpp ;
using std::cout ;
using std::endl ;
using std::vector ;

bool DOMULTI[2] ;
unsigned int C, i = 0, j = 0, k = 0, l = 0, 
	b = as<unsigned int>(B), 
	*H = new unsigned int[8] ; // # of alleles in pop. 2
double p = 0.0, q = 0.0, d, d1, d2, I, 
	maxU=sqrt(3.0), minU=-sqrt(3.0),
	tot, cov12, cov13, cov23, 
	rho12, rho13, rho23, 
	c21, c22, c31, c32, c33, 
	m = as<double>(M), // recombination rate
	r = as<double>(R), // recombination rate
	*dn = new double[2], 
	*g0 = new double[4], // initial allele frequencies
	*prob = new double[4], 
	*g = new double[8],  // allele frequencies
	*alpha = new double[4], // Dirichlet parameters
	*u = new double[3], // uniform variables
	*eta = new double[3], // noise for PSV
	*sd = new double[3] ; // desired standard deviation
vector<unsigned int> n = as< vector<unsigned int> >(N),
	 t = as< vector<unsigned int> >(T) ; 
vector<double> a = as< vector<double> >(A) ; 
vector< vector<double> > res ; // 2d vector returned to R

// Set up GSL RNG
gsl_rng_env_setup() ;
gsl_rng * rng ;
const gsl_rng_type * rngType = gsl_rng_ranlxd1 ;
rng = gsl_rng_alloc(rngType) ;

// Initialize variables
DOMULTI[0] = DOMULTI[1] = false ;
dn[0] = (double) n[0] ;
dn[1] = (double) n[1] ;
for(i = 0 ; i < 4 ; i++) alpha[i] = a[i] ;
res.resize(b) ;

// Run b simulations
for(i = 0; i < b ; i++) {
	res[i].resize(8) ;
	gsl_ran_dirichlet(rng, 4, alpha, g0) ;
	for(j = 0 ; j < 8 ; j++) g[j] = g0[j % 4] ;
	// Run PSV for t[1] generations
	for(j = 0; j < t[1]; j++) {
		// Simulate migration from pop 2 to 1 after T[0] generations
		if(j > t[0]) for(k = 0; k < 4; k++) g[k] = m*g[k+4] + (1 - m)*g[k] ;
		d1 = (g[0] * g[3]) - (g[1] * g[2]) ; // LD
		d2 = (g[4] * g[7]) - (g[5] * g[6]) ; // LD
		for(k = 0; k < 8; k++) {
			I = k % 4 == 1 || k % 4 == 2 ? 1.0 : -1.0 ; // add recombinants if hets
			g[k] = k < 4 ? g[k] + I*r*d1 : g[k] + I*r*d2 ; // recombination
			C = k < 4 ? 0 : 1 ;
			g[k] = g[k] > 0.0 ? g[k] : 0.0 ; // make sure freqs > 0
			if((g[k]*dn[C] < 3.0 || (1-g[k])*dn[C] < 3.0) && k < 4)  DOMULTI[0] = true ;
			if((g[k]*dn[C] < 3.0 || (1-g[k])*dn[C] < 3.0) && k >= 4) DOMULTI[1] = true ;
		}
		for(C = 0, k = 0 ; k < 2 ; k++, C+=4) {	
			if(!DOMULTI[k]) {
				// If there aren't low or very high allele frequencies (i.e. > 3)
				// simulate change in allele frequencies using PSV
				for(l = 0 ; l < 3 ; l++) {
					u[l] = gsl_ran_flat(rng, minU, maxU) ;
					sd[l] = sqrt(g[l+C]*(1-g[l+C])/(dn[k])) ;
				}
				cov12 = -1*g[0+C]*g[1+C]/dn[k] ;
				cov13 = -1*g[0+C]*g[2+C]/dn[k] ;
				cov23 = -1*g[1+C]*g[2+C]/dn[k] ;
				c21 = rho12 =	cov12/(sd[0]*sd[1]) ; 
				c31 = rho13 =	cov13/(sd[0]*sd[2]) ; 
				rho23 =	cov23/(sd[1]*sd[2]) ; 
				c22 = sqrt(1 - rho12*rho12) ;
				c32 = (rho23 - rho12*rho13)/c22 ;
				c33 = sqrt(1 - (c31*c31) - (c32*c32)) ;
				eta[0] = sd[0] * u[0] ;
				eta[1] = sd[1] * (c21 * u[0] + c22*u[1]) ;
				eta[2] = sd[2] * (c31 * u[0] + c32*u[1] + c33*u[2]) ;
				for(tot=0.0, l = 0 ; l < 3 ; l++) {
					g[l+C] += eta[l] ;
					tot += g[l+C] ;
				}
				g[3+C] = 1 - tot ;
			} else {
				for(l = 0 ; l < 4 ; l++) prob[l] = g[l+C] ;
				gsl_ran_multinomial(rng, 4, n[k], prob, H) ;
				for(l = 0 ; l < 4 ; l++) 
					g[k+C] = H[k]/dn[k] ;
			}
		}
		DOMULTI[0] = DOMULTI[1] = false ;
	}
	// load results 
	for(j = 0 ; j < 8 ; j++) res[i][j] = g[j] ;
}

// clean up memory
delete[] g ;
delete[] u ;
delete[] sd ;
delete[] alpha ;
delete[] eta ; 
gsl_rng_free(rng) ;

// return simulation results
return(wrap(res)) ;
