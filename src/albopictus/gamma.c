#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include "uthash.h"
#include "gamma.h"

#define MAKE_INT(x) (floor(1e8*(x)))

typedef struct {
  int n;
  int k;
  int theta;
} gamma_key_t;

typedef struct {
  gamma_key_t key;
  double value;
  UT_hash_handle hh;
} gamma_t;

gamma_t *gammas = NULL;

void gamma_destroy(void) {
  gamma_t *p, *tmp;
  HASH_ITER(hh, gammas, p, tmp) {
    HASH_DEL(gammas, p);
    free(p);
  }
  // printf("Gamma hash is cleared\n");
}

double gamma_prob(double n, double k, double theta) {
  gamma_t *elm;
  gamma_t gam;
  
  memset(&gam, 0, sizeof(gamma_t));
  gam.key.n = MAKE_INT(n);
  gam.key.k = MAKE_INT(k);
  gam.key.theta = MAKE_INT(theta);
  HASH_FIND(hh, gammas, &gam.key, sizeof(gamma_key_t), elm);

  // if found return the value
  if (elm) return elm->value;
  // if not found calculate the value and add it to the hash

  // calculate the value
  double val;
  //
  // Q(x) = \int_{x}^{+\infty} dx' p(x')
  // 1.0 - (GammaInc(k,(1.0+x)/theta) / GammaInc(k,x/theta))
  //
  if ( k > GSL_SF_GAMMA_XMAX ) {
    // Solve the problem of overflow
    double Q = gsl_cdf_gamma_Q(n,k,theta);
    if ( Q < EPS ) {
      /*
        Gamma function overflowed and the alternative formulation failed!
        Use the limit 1-e^{-1/theta}
      */
      val = 1.0 - exp(-1.0/theta);
    } else {
      val = ( ( gsl_cdf_gamma_P(n+1.0,k,theta) + Q - 1.0 ) / Q );
    }
  } else {
    val = 1.0 - (gsl_sf_gamma_inc(k,(1.0+n)/theta) / gsl_sf_gamma_inc(k,n/theta));
  }
  if (isnan(val)) val = 1.0;

  // add it to the hash
  elm = (gamma_t*)malloc( sizeof(gamma_t) );
  memset(elm, 0, sizeof(gamma_t));
  elm->key.n = MAKE_INT(n);
  elm->key.k = MAKE_INT(k);
  elm->key.theta = MAKE_INT(theta);
  elm->value = val;
  HASH_ADD(hh, gammas, key, sizeof(gamma_key_t), elm);
  
  // return the value
  return val;
}

// This is to be implemented with a hash table!
double nbinom_prob(unsigned int k, double p, double n) {
  double val;
  // p(k) = {\Gamma(n + k) \over \Gamma(k+1) \Gamma(n) } p^n (1-p)^k
  //
  // double Q = gsl_cdf_negative_binomial_Q(k,p,n);
  // val = ( ( gsl_cdf_negative_binomial_P(k+1,p,n) + Q - 1.0 ) / Q );
  //
  val = gsl_sf_beta_inc(n, 1 + k, p);
  val = (val - gsl_sf_beta_inc(n, 2 + k, p)) / (val - 1.0);
  //
  return isnan(val) ? 1.0 : val;
}

