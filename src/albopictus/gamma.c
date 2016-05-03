#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include "gamma.h"

double gamma_matrix[gamma_SIZE];
double gamma_matrix_done = 0;

void test_gamma_matrix(void) {
  double eps = 1e-6;
  int i;
  double n, mean, sd;
  double theta, k;
  double r1, r2;
  if (!gamma_matrix_done) {
    prepare_gamma();
  }
  srandom((unsigned int)time(NULL));
  for (i=0; i<1000; i++) {
    mean = (double)rand() / (double)RAND_MAX * mean_MAX;
    sd = gamma_sd * mean;
    n = random() % n_i_MAX;
    theta = sd * sd / mean;
    k = mean / theta;
    printf("Testing mean=%g, sd=%g, n=%g\n",mean,sd,n);
    r1 = prob_gamma(n,k,theta);
    r2 = prob_gamma_matrix(n,mean);
    if (fabs(r1-r2)>eps) {
      printf("Test failed: %g - %g = %g\n",r1,r2,r1-r2);
      return;
    }
    printf("Success!\n");
  }
  printf("Accuracy at epsilon = %g\n",eps);
}

void prepare_gamma(void) {
  int i;
  int n_i;
  int mean_i;
  double n;
  double mean;
  double sd;
  double theta;
  double k;
  
  for (mean_i=0; mean_i<mean_i_MAX; mean_i++) {
    for (n_i=0; n_i<n_i_MAX; n_i++) {
      i = mean_i * row_SIZE + n_i;
      mean = mean_i*mean_STEP;
      if (mean==0) {
        gamma_matrix[i] = 1.0;
        continue;
      }
      n = n_i*n_STEP;
      sd = gamma_sd * mean;
      theta = sd * sd / mean;
      k = mean / theta;
      gamma_matrix[i] = prob_gamma(n,k,theta);
      // printf("mean=%g, n=%g, i=%d, gamma=%g\n",mean,n,i,gamma_matrix[i]);
    }
  }
  
  gamma_matrix_done = 1;
  printf("Gamma matrix is ready to use!\n");
}

double prob_gamma_matrix(double n, double mean) {
  double m1 = floor(mean / mean_STEP);
  double m2 = ceil(mean / mean_STEP);
  int i1 = m1 * row_SIZE + floor(n / n_STEP);
  int i2 = m2 * row_SIZE + floor(n / n_STEP);
  if (i1>=gamma_SIZE || i2>=gamma_SIZE) {
    printf("ERROR: Not in matrix (mean=%g, n=%g)\n",mean,n);
    exit(1);
  }
  return (i1==i2) ? gamma_matrix[i1] : gamma_matrix[i1] + (mean/mean_STEP - m1)*(gamma_matrix[i2]-gamma_matrix[i1])/(m2-m1); // Linear interpolation
}

double prob_gamma(double n, double k, double theta) {
  double val;
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
  return isnan(val) ? 1.0 : val;
}

double prob_nbinom(unsigned int k, double p, double n) {
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

