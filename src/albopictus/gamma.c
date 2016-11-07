#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include "uthash.h"
#include "gamma.h"

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

#define MAKE_INT(x) (round(1e3*(x)))
#define MAKE_DOUBLE(x) (1e-3*(x))

#define MAX_MEM   1000000000

typedef struct {
  double n;
  double value;
  UT_hash_handle hh;
} gamma_n_t;

typedef struct {
  double sd;
  gamma_n_t *n_hash;
  UT_hash_handle hh;
} gamma_sd_t;

typedef struct {
  double mean;
  gamma_sd_t *sd_hash;
  UT_hash_handle hh;
} gamma_mean_t;

gamma_mean_t *gamma_mean = NULL;
size_t gamma_mean_sd_mem = 0;

void gamma_dist_destroy(void) {
  gamma_mean_t *p, *tmp;
  gamma_sd_t *p2, *tmp2;
  gamma_n_t *p3, *tmp3;
  HASH_ITER(hh, gamma_mean, p, tmp) {
    HASH_ITER(hh, p->sd_hash, p2, tmp2) {
      HASH_ITER(hh, p2->n_hash, p3, tmp3) {
        HASH_DEL(p2->n_hash, p3);
        free(p3);
      }
      HASH_DEL(p->sd_hash, p2);
      free(p2);
    }
    HASH_DEL(gamma_mean, p);
    free(p);
  }
  gamma_mean_sd_mem = 0;
  gamma_mean = NULL;
  //
  printf("Gamma hash is successfully cleared\n");
}

void gamma_dist_check(void) {
  if (gamma_mean_sd_mem > MAX_MEM)
    gamma_dist_destroy();
}

double gamma_pdf(double n, double mean, double sd) {
  double theta = sd * sd / mean;
  double k = mean / theta;
  return gsl_ran_gamma_pdf(n,k,theta);
}

double gamma_dist_prob(double mean, double sd, double n) {
  gsl_set_error_handler_off();
  double tmp;
  //
  double theta = sd * sd / mean;
  double k = mean / theta;
  //
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
    tmp = gsl_sf_gamma_inc(k,(1.0+n)/theta);
    if (tmp == GSL_EMAXITER) {
      printf("WARNING: gsl_sf_gamma_inc returned GSK_EMAXITER\n");
      gsl_set_error_handler(NULL);
      return 1.0;
    }
    val = tmp;
    tmp = gsl_sf_gamma_inc(k,n/theta);
    if (tmp == GSL_EMAXITER) {
      printf("WARNING: gsl_sf_gamma_inc returned GSK_EMAXITER\n");
      gsl_set_error_handler(NULL);
      return 1.0;
    }
    val = 1.0 - (val / tmp);
  }
  if (isnan(val)) val = 1.0;

  // return the value
  gsl_set_error_handler(NULL);
  return val;
}

char gamma_dist_hash(double mean_d, double sd_d, double n_d, double *value) {
  (*value) = 0;
  double mean = MAKE_INT(mean_d);
  double sd = MAKE_INT(sd_d);
  double n = MAKE_INT(n_d);
  //
  gamma_mean_t *g_mean;
  HASH_FIND(hh, gamma_mean, &mean, sizeof(double), g_mean);
  if (!g_mean) {
    g_mean = (gamma_mean_t*)malloc(sizeof(gamma_mean_t));
    if (g_mean == NULL) return 0;
    g_mean->mean = mean;
    g_mean->sd_hash = NULL;
    HASH_ADD(hh, gamma_mean, mean, sizeof(double), g_mean);
    gamma_mean_sd_mem += sizeof(gamma_mean_t);
  }
  //
  gamma_sd_t *g_sd;
  HASH_FIND(hh, g_mean->sd_hash, &sd, sizeof(double), g_sd);
  if (!g_sd) {
    g_sd = (gamma_sd_t*)malloc(sizeof(gamma_sd_t));
    if (g_sd == NULL) return 0;
    g_sd->sd = sd;
    g_sd->n_hash = NULL;
    HASH_ADD(hh, g_mean->sd_hash, sd, sizeof(double), g_sd);
    gamma_mean_sd_mem += sizeof(gamma_sd_t);
  }
  //
  gamma_n_t *g_n;
  HASH_FIND(hh, g_sd->n_hash, &n, sizeof(double), g_n);
  if (!g_n) {
    g_n = (gamma_n_t*)malloc(sizeof(gamma_n_t));
    if (g_n == NULL) return 0;
    g_n->n = n;
    g_n->value = gamma_dist_prob(MAKE_DOUBLE(mean),MAKE_DOUBLE(sd),MAKE_DOUBLE(n)); // precision is already lost with the keys
    HASH_ADD(hh, g_sd->n_hash, n, sizeof(double), g_n);
    gamma_mean_sd_mem += sizeof(gamma_n_t);
  }
  //
  (*value) = g_n->value;
  return 1;
}

/*
  ---------------------------------------------------------------------------------
*/

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

/*
  ---------------------------------------------------------------------------------
*/

double gamma_matrix[gamma_SIZE];
double gamma_matrix_done = 0;

void prepare_gamma_matrix(void) {
  int i;
  int n_i;
  int mean_i;
  double n;
  double mean;
  double sd;
  
  for (mean_i=0; mean_i<mean_i_MAX; mean_i++) {
    for (n_i=0; n_i<n_i_MAX; n_i++) {
      i = mean_i * row_SIZE + n_i;
      mean = mean_i*mean_STEP;
      if (mean==0) {
        gamma_matrix[i] = 1.0;
        continue;
      }
      n = n_i*n_STEP;
      sd = gamma_matrix_sd * mean;
      gamma_matrix[i] = gamma_dist_prob(mean,sd,n);
      // printf("mean=%g, n=%g, i=%d, gamma=%g\n",mean,n,i,gamma_matrix[i]);
    }
  }
  
  gamma_matrix_done = 1;
  printf("Gamma matrix is ready to use!\n");
}

double gamma_dist_matrix(double mean, double n) {
  if (!gamma_matrix_done)
    prepare_gamma_matrix();
  //
  double m1 = floor(mean / mean_STEP);
  double m2 = ceil(mean / mean_STEP);
  int i1 = m1 * row_SIZE + floor(n / n_STEP);
  int i2 = m2 * row_SIZE + floor(n / n_STEP);
  if (i1>=gamma_SIZE || i2>=gamma_SIZE) {
    double tmp;
    if (gamma_dist_hash(mean,gamma_matrix_sd*mean,n,&tmp)) {
      printf("WARNING: Gamma hash is used instead of the matrix!\n");
      return tmp;
    } else {
      printf("ERROR: Gamma distribution failed!\n");
      exit(1);
    }
  }
  return (i1==i2) ? gamma_matrix[i1] : gamma_matrix[i1] + (mean/mean_STEP - m1)*(gamma_matrix[i2]-gamma_matrix[i1])/(m2-m1); // Linear interpolation
}

