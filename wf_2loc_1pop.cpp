using namespace Rcpp ;
using std::vector ;

unsigned int  i = 0, j = 0, k = 0,
							b = as<unsigned int>(B),
							t = as<unsigned int>(T), 
						 *H = new unsigned int[4] ;
double p = 0.0, q = 0.0, d = 0.0, 
			 n = as<double>(N), 
			 r = as<double>(R),
			 u = as<double>(U),
			 *prob = new double[4], 
			 *buff = new double[4] ;
vector<double> g = as< vector<double> >(G) ;
vector< vector<double> > res ; 
  gsl_rng_env_setup() ;
  gsl_rng * rng ;
const gsl_rng_type * rngType ;
  
rngType = gsl_rng_ranlxd1 ;
rng = gsl_rng_alloc(rngType) ;
for(i = 0 ; i < 4 ; i++) buff[i] = g[i] ;
res.resize(b) ;
for(i = 0; i < b ; i++) {
	res[i].resize(4) ;
	for(j = 0; j < t; j++) {
		d = (g[0] * g[3]) - (g[1] * g[2]) ;
		prob[0] = (1-u)*(g[0] - r*d) - 2*u*g[0] ;
		prob[1] = (1-u)*(g[1] + r*d) + u*(g[0] - g[1]) ;
		prob[2] = (1-u)*(g[2] + r*d) + u*(g[0] - g[2]) ;
		prob[3] = (1-u)*(g[3] - r*d) + u*(g[1] + g[2]) ;
		gsl_ran_multinomial(rng, 4, n, prob, H) ;
		for(k = 0 ; k < 4 ; k++) g[k] = ((double)H[k])/n ;
	}
	for(j = 0 ; j < 4 ; j++) {
		res[i][j] = g[j] ;
		g[j] = buff[j] ;
	}
}
delete[] H ;
delete[] buff ;
delete[] prob ; 
gsl_rng_free(rng) ;

return(wrap(res)) ;
