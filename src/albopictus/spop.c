#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include "gamma.h"
#include "spop.h"
#include "ran_gen.h"

#define MAX_DAYS 1000

extern gsl_rng *RAND_GSL;

// This is for optimisation
void set_prdev(spop s) {
  int n;
  for (n=0; n<MAX_DAYS; n++) {
    // Probability of developing from day tmpn->development until tmpn->development+1
    s->prdev[n] = gamma_mean_prob(s->development_time,s->development_time_sd,n);
    // printf("%d %g\n",n,s->prdev[n]);
  }
}
  
spop spop_init(double surv, double dev_time_mean, double dev_time_sd) {
  spop pop = (spop)malloc(sizeof(struct population_st));
  people_data tmp = (people_data)malloc(sizeof(struct people_st));
  tmp->development = 0;
  tmp->batchsize = 0;
  tmp->next = 0;
  pop->people = tmp;
  pop->development_time = dev_time_mean;
  pop->development_time_sd = dev_time_sd;
  pop->survival = surv;
  pop->popsize = 0;
  pop->prdev = (double *)malloc(MAX_DAYS*sizeof(double));
  set_prdev(pop);
  return pop;
}

void spop_setdevtime(spop s, double dev_time_mean, double dev_time_sd) {
  s->development_time = dev_time_mean;
  s->development_time_sd = dev_time_sd;
  set_prdev(s);
}

void spop_setsurv(spop s, double surv) {
  s->survival = surv;
}

void spop_add(spop s, int dev, int size){
  people_data tmp;
  for (tmp = s->people;
       tmp && tmp->next;
       tmp = tmp->next) {
    if (tmp->next->development == dev) {
      tmp->next->batchsize += size;
      s->popsize += size;
      return;
    }
  }
  tmp = (people_data)malloc(sizeof(struct people_st));
  tmp->development = dev;
  tmp->batchsize = size;
  tmp->next  = s->people->next;
  s->people->next = tmp;
  s->popsize += tmp->batchsize;
}

void spop_remove_next(spop s, people_data p){
  if (!p) return;
  people_data tmp = p->next;
  p->next = tmp->next;
  s->popsize -= tmp->batchsize;
  free(tmp);
}

void spop_empty(spop s) {
  while (s->popsize && s->people->next) {
    spop_remove_next(s,s->people);
  }
  s->popsize = 0;
}

void spop_destroy(spop *s) {
  spop_empty(*s);
  free((*s)->people);
  free((*s)->prdev);
  free((*s));
}

void spop_print(spop s) {
  people_data p = s->people->next;
  int count;
  printf("/------------------>\n");
  printf("Population size: %d\n",s->popsize);
  printf("Development time: %g\n",s->development_time);
  printf("Development time (std): %g\n",s->development_time_sd);
  printf("Daily survival: %g\n",s->survival);
  for (count=1; p; count++) {
    printf("%d\t%d\t%d\n",count,p->development,p->batchsize);
    p = p->next;
  }
  printf("\\------------------>\n");
}

int spop_survive(spop s) {
  int ret = 0;
  char remove = 0;
  int k;
  people_data tmp, tmpn;
  for (tmp = s->people;
       tmp && tmp->next;
       ) {
    tmpn = tmp->next;
    // Survive
    k = tmpn->batchsize - gsl_ran_binomial(RAND_GSL,s->survival,tmpn->batchsize);
    // printf("%d will die from %g,%d,%d\n",k,s->survival,tmpn->development,tmpn->batchsize);
    tmpn->batchsize -= k;
    s->popsize -= k;
    if (tmpn->batchsize==0) {
      remove = 1;
    } else if (tmpn->batchsize<0) {
      printf("ERROR: Error in population!\n");
      spop_print(s);
      exit(1);
    } else {
      // Develop
      k = gsl_ran_binomial(RAND_GSL,
                           s->prdev[tmpn->development],
                           tmpn->batchsize);
      // printf("%d will develop with %g probability in %d,%d\n",k,s->prdev[tmpn->development],tmpn->development,tmpn->batchsize);
      tmpn->development += 1;
      tmpn->batchsize -= k;
      s->popsize -= k;
      ret += k;
      if (tmpn->batchsize==0) {
        remove = 1;
      } else if (tmpn->batchsize<0) {
        printf("ERROR: Error in population!\n");
        spop_print(s);
        exit(1);
      }
      if (tmpn->development >= MAX_DAYS) {
        printf("ERROR: Development time is too high!\n");
        spop_print(s);
        exit(1);
      }
    }
    if (remove) {
      tmp->next = tmpn->next;
      free(tmpn);
      if (s->popsize == 0) break;
      remove = 0;
    } else {
      tmp = tmp->next;
    }
  }
  return ret;
}

void spop_swap(spop pfrom, spop pto) {
  if (pfrom->popsize <= 0) return;
  int i = gsl_rng_uniform_int(RAND_GSL,pfrom->popsize)+1;
  people_data tmp;
  for (tmp = pfrom->people;
       tmp && tmp->next;
       tmp = tmp->next) {
    if ((i -= tmp->next->batchsize) <= 0) break;
  };
  spop_add(pto,tmp->next->development,1);
  tmp->next->batchsize -= 1;
  pfrom->popsize -= 1;
  if (tmp->next->batchsize <= 0)
    spop_remove_next(pfrom,tmp);
}
