/***************************************************************************************
 *  Multivariate Normal density function and random number generator
 *  Multivariate Student t density function and random number generator
 *  Wishart random number generator
 *  Using GSL -> www.gnu.org/software/gsl
 *
 *  Copyright (C) 2006  Ralph dos Santos Silva
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  AUTHOR
 *     Ralph dos Santos Silva,  [EMAIL PROTECTED]
 *     March, 2006      
 ***************************************************************************************/

int rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *r){
  /* multivariate normal distribution random number generator */
  /*
   *    n       dimension of the random vetor
   *    mean    vector of means of size n
   *    var     variance matrix of dimension n x n
   *    result  output variable with a sigle random vector normal distribution generation
   */
  int k;
  gsl_matrix *work = gsl_matrix_alloc(n,n);
 
  gsl_matrix_memcpy(work,var);
  gsl_linalg_cholesky_decomp(work);
 
  for(k=0; k<n; k++)
    gsl_vector_set( result, k, gsl_ran_ugaussian(r) );
 
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
  gsl_vector_add(result,mean);
 
  gsl_matrix_free(work);
 
  return 0;
}
