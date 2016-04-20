#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_cdf.h>
#include "incubator.h"

#define EPS 1e-14

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

void incubator_develop(incubator *s, double time, double lambda, double sd, double *val) {
  /*
    Calculate daily survival according to 
    1. expected life span (lambda) and
    2. current age ((*s)->data.popdev)
   */
  double theta = sd * sd / lambda;
  double k = lambda / theta;
  double p4;
  //
  *val = 0; // total population size
  incubator dummy = (incubator)malloc(sizeof(struct incubator_st));
  dummy->next = *s;
  incubator prev = dummy;
  while (*s) {
    // Probability of dying between day d-1 and d
    // Q(x) = \int_{x}^{+\infty} dx' p(x')
    p4 = (gsl_cdf_gamma_P((*s)->data.popdev + time, k, theta)
          - gsl_cdf_gamma_P((*s)->data.popdev, k, theta))
      / gsl_cdf_gamma_Q((*s)->data.popdev, k, theta);
    //
    (*s)->data.popsize *= 1.0-p4;
    if ((*s)->data.popsize < EPS) { // remove this batch
      (*s) = (*s)->next;
      free(prev->next);
      prev->next = (*s);
    } else { // move on to the next batch
      *val += (*s)->data.popsize;
      (*s)->data.popdev += time;
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
