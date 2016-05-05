#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include "uthash.h"
#include "gamma.h"

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

#define MAKE_INT(x) (round(1e3*(x)))
#define MAKE_DOUBLE(x) (1e-3*(x))

#define MAX_COUNT INT_MAX
#define MIN_COUNT 1000
#define MAX_MEM   100000000

typedef struct {
  double n;
  double value;
  UT_hash_handle hh;
} gamma_n_t;

typedef struct {
  double sd;
  size_t count;
  gamma_n_t *n_hash;
  UT_hash_handle hh;
} gamma_sd_t;

typedef struct {
  double mean;
  size_t count;
  gamma_sd_t *sd_hash;
  UT_hash_handle hh;
} gamma_mean_t;

gamma_mean_t *gamma_mean = NULL;
size_t gamma_mean_sd_mem = 0;

void gamma_mean_destroy(void) {
  gamma_mean_t *p, *tmp;
  gamma_sd_t *p2, *tmp2;
  gamma_n_t *p3, *tmp3;
  HASH_ITER(hh, gamma_mean, p, tmp) {
    printf("gamma_mean = %g, count = %zu\n",p->mean,p->count);
    HASH_ITER(hh, p->sd_hash, p2, tmp2) {
      printf("gamma_mean = %g, count = %zu, gamma_sd = %g, count = %zu\n",p->mean,p->count,p2->sd,p2->count);
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
  printf("Gamma hash is successfully cleared\n");
}

void gamma_mean_trim(void) {
  printf("Memory used by the Gamma hash: %zu\n",gamma_mean_sd_mem);
  if (gamma_mean_sd_mem < MAX_MEM) return;
  //
  gamma_mean_t *p, *tmp;
  gamma_sd_t *p2, *tmp2;
  gamma_n_t *p3, *tmp3;
  HASH_ITER(hh, gamma_mean, p, tmp) {
    if (p->count < MIN_COUNT) {
      HASH_ITER(hh, p->sd_hash, p2, tmp2) {
        if (p2->count < MIN_COUNT) {
          if ((p->count)-- == 0) { // gamma_mean is checked in the next cal to this function
            printf("ERROR trimming the gamma hash!\n");
            exit(1);
          }
          HASH_ITER(hh, p2->n_hash, p3, tmp3) {
            HASH_DEL(p2->n_hash, p3);
            free(p3);
            if (gamma_mean_sd_mem < sizeof(gamma_n_t)) {
              printf("ERROR trimming the gamma hash (n)!\n");
              exit(1);
            }
            gamma_mean_sd_mem -= sizeof(gamma_n_t);
          }
          HASH_DEL(p->sd_hash, p2);
          free(p2);
          if (gamma_mean_sd_mem < sizeof(gamma_sd_t)) {
            printf("ERROR trimming the gamma hash (sd)!\n");
            exit(1);
          }
          gamma_mean_sd_mem -= sizeof(gamma_sd_t);
        }
      }
      HASH_DEL(gamma_mean, p);
      free(p);
      if (gamma_mean_sd_mem < sizeof(gamma_mean_t)) {
        printf("ERROR trimming the gamma hash (mean)!\n");
        exit(1);
      }
      gamma_mean_sd_mem -= sizeof(gamma_mean_t);
    }
  }
  printf("Gamma hash is trimmed down = %zu\n",gamma_mean_sd_mem);
}

double gamma_mean_prob(double mean, double sd, double n) {
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

char gamma_mean_hash(double mean_d, double sd_d, double n_d, double *value) {
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
    g_mean->count = 1;
    g_mean->sd_hash = NULL;
    HASH_ADD(hh, gamma_mean, mean, sizeof(double), g_mean);
    gamma_mean_sd_mem += sizeof(gamma_mean_t);
  } else {
    g_mean->count = min( g_mean->count+1, MAX_COUNT );
  }
  //
  gamma_sd_t *g_sd;
  HASH_FIND(hh, g_mean->sd_hash, &sd, sizeof(double), g_sd);
  if (!g_sd) {
    g_sd = (gamma_sd_t*)malloc(sizeof(gamma_sd_t));
    if (g_sd == NULL) return 0;
    g_sd->sd = sd;
    g_sd->count = 1;
    g_sd->n_hash = NULL;
    HASH_ADD(hh, g_mean->sd_hash, sd, sizeof(double), g_sd);
    gamma_mean_sd_mem += sizeof(gamma_sd_t);
  } else {
    g_sd->count = min( g_sd->count+1, MAX_COUNT );
  }
  //
  gamma_n_t *g_n;
  HASH_FIND(hh, g_sd->n_hash, &n, sizeof(double), g_n);
  if (!g_n) {
    g_n = (gamma_n_t*)malloc(sizeof(gamma_n_t));
    if (g_n == NULL) return 0;
    g_n->n = n;
    //g_n->value = gamma_mean_prob(mean,sd,n);
    g_n->value = gamma_mean_prob(MAKE_DOUBLE(mean),MAKE_DOUBLE(sd),MAKE_DOUBLE(n)); // precision is already lost with the keys
    HASH_ADD(hh, g_sd->n_hash, n, sizeof(double), g_n);
    gamma_mean_sd_mem += sizeof(gamma_n_t);
  }
  //
  (*value) = g_n->value;
  return 1;
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

