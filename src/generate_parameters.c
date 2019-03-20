#ifdef HAVE_GSL

#include <gsl/gsl_randist.h>
#include <petscmat.h>

PetscErrorCode GeneratePermeability(PetscReal *parameters ,const PetscInt n, const double mu, const double sigma, const double scale){
  /*
   * Using GSL for generating array of values which will be interpreted as pemeability parametrs on finite elements.
   * The code here presumes that the size of the array parameters is n**2, there is no checking.
   * Mu and sigma are parameters of the lognormal distribution with pd
   * p(x) dx = {1 \over x \sqrt{2 \pi \sigma^2} } \exp(-(\ln(x) - \zeta)^2/2 \sigma^2) dx
   * copied from
   * http://www.gnu.org/software/gsl/manual/html_node/The-Lognormal-Distribution.html#The-Lognormal-Distribution
   * The parameter scale effectively multiplies this distribution.
   * */
  gsl_rng* r;
  double p;
  PetscInt i;
  PetscFunctionBeginUser;
  /* setting type of random generator, taken from
  * http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html#Random-number-generator-algorithms */
  r = gsl_rng_alloc(gsl_rng_mt19937);
  /* seed setting, 0 means default seed see
  * http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-initialization.html#Random-number-generator-initialization */
  gsl_rng_set(r,1);
  for( i = 0; i < n*n; i++){
    p = gsl_ran_lognormal(r,mu,sigma);
    parameters[i] = scale*p;
  }
  gsl_rng_free(r);
  PetscFunctionReturn(0);
}
#endif

