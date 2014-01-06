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

// [[Rcpp::export]]
unsigned int C, i = 0, j = 0, k = 0, l = 0, 
b = as<unsigned int>(B) ; // # of 2-locus pairs
double p = 0.0, q = 0.0, d, d1, d2, I, 
	maxU=sqrt(3.0), minU=-sqrt(3.0),
	tot, tmp1, tmp2, tmp3, 
	m = as<double>(M), // recombination rate
	r = as<double>(R), // recombination rate
	*g0 = new double[4], // initial allele frequencies
	*g = new double[8],  // allele frequencies
	*alpha = new double[4], // Dirichlet parameters
	*u = new double[3], // uniform variables
	*eta = new double[3], // noise for PSV
	*sd = new double[3], // desired standard deviation
	*c2 = new double[2],
	*c3 = new double[3] ;
vector<unsigned int> n = as< vector<unsigned int> >(N),
	 t = as< vector<unsigned int> >(T) ; 
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
	res[i].resize(8) ;
	gsl_ran_dirichlet(rng, 4, alpha, g0) ;
	for(j = 0 ; j < 8 ; j++) g[j] = g0[j % 4] ;
	// Run PSV for T[1] generations
	for(j = 0; j < t[1]; j++) {
		// Simulate migration from pop 2 to 1 after T[0] generations
		if(j > t[0]) for(k = 0; k < 4; k++) g[k] = m*g[k+4] + (1 - m)*g[k] ;
		d1 = (g[0] * g[3]) - (g[1] * g[2]) ; // LD
		d2 = (g[4] * g[7]) - (g[5] * g[6]) ; // LD
		for(k = 0; k < 8; k++) {
			I = k % 4 == 1 || k % 4 == 2 ? 1.0 : -1.0 ; // add recombinants if hets
			g[k] = k < 4 ? g[k] + I*r*d1 : g[k] + I*r*d2 ; // recombination
			g[k] = g[k] > 0.0 ? g[k] : 0.0 ; // make sure freqs > 0
		}
		for(C = 0, k = 0 ; k < 2 ; k++, C+=4) {	
			for(l = 0 ; l < 3 ; l++) {
				u[l] = gsl_ran_flat(rng, minU, maxU) ;
				sd[l] = sqrt(g[l+C]*(1-g[l+C])/((double) n[k])) ;
			}
			tmp1 = g[0+C]*g[1+C]/((1-g[0+C])*(1-g[1+C])) ;
			c2[0] = -sqrt(tmp1) ;
			c2[1] = sqrt(1-tmp1) ;
			tmp1 = g[0+C]*g[2+C]/((1-g[0+C])*(1-g[2+C])) ;
			tmp2 = g[1+C]*g[2+C]/((1-g[1+C])*(1-g[2+C])) ;
			tmp3 = 1-g[0+C] ;
			c3[0] = -sqrt(tmp1) ;
			c3[1] = -sqrt(tmp2)/tmp3 ;  
			c3[2] = sqrt(1 - tmp1 - tmp2/(tmp3*tmp3)) ;
			eta[0] = sd[0] * u[0] ;
			eta[1] = sd[1] * (c2[0] * u[0] + c2[1]*u[1]) ;
			eta[2] = sd[2] * (c3[0] * u[0] + c3[1]*u[1] + c3[2]*u[2]) ;
			for(tot=0.0, l = 0 ; l < 3 ; l++) {
				g[l+C] += eta[l] ;
				tot += g[l+C] ;
			}
			g[3+C] = 1 - tot ;
		}
	}
	// load results 
	for(j = 0 ; j < 8 ; j++) res[i][j] = g[j] ;
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
