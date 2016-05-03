#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include "gamma.h"
#include "incubator.h"

extern double gamma_matrix[gamma_SIZE];
extern double gamma_matrix_done;

volatile clock_t istart = 0, idiff;
double itime2here(void) {
  idiff = clock() - istart;
  istart = clock();
  return idiff;
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
