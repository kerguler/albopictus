#include "Python.h"
#include <time.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include "incubator.h"

double gamma_matrix[gamma_SIZE];
double gamma_matrix_done = 0;

volatile clock_t istart = 0, idiff;
double itime2here(void) {
  idiff = clock() - istart;
  istart = clock();
  return idiff;
}

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

void incubator_empty(incubator *s) {
  incubator tmp;
  while ((*s)) {
    tmp = (*s)->next;
    free((*s));
    (*s) = tmp;
  }
}

void incubator_add(incubator *s, double popsize, double popdev) {
  if (popsize<=0) return;
  incubator tmp = (incubator)malloc(sizeof(struct incubator_st));
  tmp->data.popsize = popsize;
  tmp->data.popdev = popdev;
  tmp->next  = *s;
  *s = tmp;
}

void incubator_update(incubator s, void (*fun)(double *size, double *dev, double *par), double *par) {
  while (s) {
    fun(&(s->data.popsize),&(s->data.popdev),par);
    s = s->next;
  }
}

void incubator_remove(incubator *s, double *val) {
  *val = 0;
  incubator dummy = (incubator)malloc(sizeof(struct incubator_st));
  dummy->next = *s;
  incubator prev = dummy;
  while (*s) {
    if ((*s)->data.popdev >= 1.0-EPS) {
      *val += (*s)->data.popsize;
      (*s) = (*s)->next;
      free(prev->next);
      prev->next = (*s);
    } else {
      prev = (*s);
      (*s) = (*s)->next;
    }
  }
  (*s) = dummy->next;
  free(dummy);
}

void incubator_develop_survive(incubator *s,
                               double d_mean, // development time (ignore if < 0)
                               double d_sd, // 
                               double p_mean, // survival time (ignore if < 0)
                               double p_sd, // 
                               double juvthr, // days to become mature
                               double *juv, // juveniles
                               double *avec, // matures
                               double *dev, // developed
                               char mode) { // option for gamma distribution
  /*
    Calculate daily survival and development according to 
    1. expected life span (d_mean),
    2. expected development time (p_mean) and
    3. current age ((*s)->data.popdev)
   */
  if (mode == 2 && !gamma_matrix_done) {
    prepare_gamma();
  }
  //
  double d_theta, d_k, d_p, d_r;
  double p_theta, p_k, p_p, p_r;
  if (mode == 1) {
    d_p = p_p = 0.5;
    d_r = d_mean;
    p_r = p_mean;
  } else if (mode == 0) {
    d_theta = d_sd * d_sd / d_mean;
    d_k = d_mean / d_theta;
    p_theta = p_sd * p_sd / p_mean;
    p_k = p_mean / p_theta;
  }
  //
  *avec = 0;
  *juv = 0;
  if (d_mean >= 0)
    *dev = 0;
  incubator dummy = (incubator)malloc(sizeof(struct incubator_st));
  dummy->next = *s;
  incubator prev = dummy;
  double prob = 0.0;
  while (*s) {
    // Probability of dying or developing between day d and d+1
    //
    // survival
    if (p_mean >= 0) {
      if (mode == 1)
        prob = 1.0 - prob_nbinom((unsigned int)floor((*s)->data.popdev),p_p,p_r);
      else if (mode == 0)
        prob = 1.0 - prob_gamma((*s)->data.popdev,p_k,p_theta);
      else if (mode == 2)
        prob = 1.0 - prob_gamma_matrix((*s)->data.popdev,p_mean);
      (*s)->data.popsize *= prob;
      //
      if ((*s)->data.popsize < EPS) { // all dead, remove this batch
        (*s) = (*s)->next;
        free(prev->next);
        prev->next = (*s);
        continue;
      }
    }
    // development
    if (d_mean >= 0) {
      if (mode == 1)
        prob = prob_nbinom((unsigned int)floor((*s)->data.popdev),d_p,d_r);
      else if (mode == 0)
        prob = prob_gamma((*s)->data.popdev,d_k,d_theta);
      else if (mode == 2)
        prob = prob_gamma_matrix((*s)->data.popdev,d_mean);
      (*dev) = (*s)->data.popsize * prob;
      (*s)->data.popsize *= 1.0 - prob;
      //
      if ((*s)->data.popsize < EPS) { // all developed, remove this batch
        (*s) = (*s)->next;
        free(prev->next);
        prev->next = (*s);
        continue;
      }
    }
    // move on to the next batch
    (*s)->data.popdev += 1.0;
    if ((*s)->data.popdev < juvthr) // juveniles
      (*juv) += (*s)->data.popsize;
    else // matures
      (*avec) += (*s)->data.popsize;
    prev = (*s);
    (*s) = (*s)->next;
  }
  (*s) = dummy->next;
  free(dummy);
}

double incubator_sum(incubator s) {
  double sum = 0;
  while (s) {
    sum += s->data.popsize;
    s = s->next;
  }
  return sum;
}

void incubator_print(incubator s) {
  printf("/------------------>\n");
  while (s) {
    printf("%g\t%g\n",s->data.popsize,s->data.popdev);
    s = s->next;
  }
  printf("\\------------------>\n");
}
