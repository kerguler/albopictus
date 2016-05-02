#include "Python.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include "incubator.h"

#define EPS 1e-14

#define n_MAX 400.0
#define n_STEP 1.0
#define n_i_MAX 400
#define mean_MAX 401.0
#define mean_STEP 0.01
#define mean_i_MAX 40100
#define gamma_SIZE 16040000
#define row_SIZE 400

double gamma_matrix[gamma_SIZE];
double gamma_matrix_done = 0;

volatile clock_t istart = 0, idiff;
double itime2here() {
  idiff = clock() - istart;
  istart = clock();
  return idiff;
}

void test() {
  double n=50, mean = 137.365, sd = 0.375*137.365;
  double theta, k, p, r;
  p = 0.5;
  r = mean;
  theta = sd * sd / mean;
  k = mean / theta;
  /*
  printf("Test:\nmean:%g,sd:%g,n:%g\np:%g,r:%g\nk:%g,theta:%g\n%g,%g,%g\n",
         mean,
         sd,
         n,
         p,
         r,
         k,
         theta,
         prob_gamma(n,k,theta),
         prob_nbinom((unsigned int)floor(n),p,r),
         prob_gamma_matrix(n,mean));
  */
  for (n=0;n<mean;n++) {
    printf("%g,%g,%g,%g\n",
           n,
           prob_gamma(n,k,theta),
           prob_nbinom((unsigned int)floor(n),p,r),
           prob_gamma_matrix(n,mean));
  }
}

void prepare_gamma() {
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
      sd = 0.375 * mean;
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
  int i = floor(mean / mean_STEP) * row_SIZE + floor(n / n_STEP);
  if (i>=gamma_SIZE) {
    printf("ERROR: Not in matrix (mean=%g, n=%g)\n",mean,n);
    exit(1);
  }
  return gamma_matrix[i];
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

void incubator_develop(incubator *s, double mean, double sd, double *val, double thr, double *valthr, char mode) { // This is custom-made for improved performance
  //double timeall = 0;
  //double timeof = 0;
  /*
    Calculate daily survival according to 
    1. expected life span (mean) and
    2. current age ((*s)->data.popdev)
   */
  if (mode == 2 && !gamma_matrix_done) {
    prepare_gamma();
  }
  //
  double theta, k, p, r;
  if (mode == 1) {
    //p = mean / sd / sd;
    //r = mean * p / (1.0 - p);
    p = 0.5;
    r = mean;
    //printf("NBIN: %g,%g,%g,%g\n",mean,sd,p,r);
  } else if (mode == 0) {
    //sd = sqrt(2.0*mean);
    theta = sd * sd / mean;
    k = mean / theta;
  }
  //
  *val = 0; // total population size
  *valthr = 0; // total (juvenile) population size
  incubator dummy = (incubator)malloc(sizeof(struct incubator_st));
  dummy->next = *s;
  incubator prev = dummy;
  double prob;
  while (*s) {
    // Probability of dying between day d and d+1
    //
    if (mode == 1)
      prob = prob_nbinom((unsigned int)floor((*s)->data.popdev),p,r);
    else if (mode == 0)
      prob = prob_gamma((*s)->data.popdev,k,theta);
    else if (mode == 2)
      prob = prob_gamma_matrix((*s)->data.popdev,mean);
    (*s)->data.popsize *= 1.0 - prob;
    //
    if ((*s)->data.popsize < EPS) { // remove this batch
      (*s) = (*s)->next;
      free(prev->next);
      prev->next = (*s);
    } else { // move on to the next batch
      (*s)->data.popdev += 1.0;
      if ((*s)->data.popdev < thr) // Juvenile adults
        *valthr += (*s)->data.popsize;
      else // Mature adults
        *val += (*s)->data.popsize;        
      prev = (*s);
      (*s) = (*s)->next;
    }
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
