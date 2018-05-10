#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include "ran_gen.h"
#include "gamma.h"
#include "sdpop.h"

// Memory monitoring
/*
#define calloc(number,size) printmem(1,(void *)calloc((number),(size)),(number)*(size),__FILE__,__LINE__)
#define malloc(size) printmem(1,(void *)malloc((size)),(size),__FILE__,__LINE__)
#define free(pointer) {printmem(0,0,sizeof(*(pointer)),__FILE__,__LINE__); free((pointer));}
void *printmem(int, void *, size_t, char *, int);

void *printmem(int type, void *pointer, size_t size, char *filen, int linen) {
  static size_t total = 0;
  total += (type?1:-1)*size;
  printf("%d: %s: %d: %d: %d\n",type,filen,linen,(int)(size),(int)(total));
  return pointer;
}
*/
// ---

extern gsl_rng *RAND_GSL;

sdpop sdpop_init(unsigned char stochastic, unsigned char gamma_mode) {
  sdpop pop = (sdpop)malloc(sizeof(struct population_st));
  individual_data tmp = (individual_data)malloc(sizeof(struct individual_st));
  tmp->next = 0;
  tmp->age = 0;
  if (stochastic) {
    tmp->number.i = 0;
    pop->size.i = 0;
    pop->dead.i = 0;
    pop->developed.i = 0;
  } else {
    tmp->number.d = 0.0;
    pop->size.d = 0.0;
    pop->dead.d = 0.0;
    pop->developed.d = 0.0;
  }
  pop->individual = tmp;
  pop->gamma_mode = gamma_mode;
  pop->stochastic = stochastic;
  return pop;
}

void sdpop_sdadd(sdpop s, unsigned int age, sdnum number){
  individual_data tmp;
  for (tmp = s->individual;
       tmp && tmp->next;
       tmp = tmp->next) {
    if (tmp->next->age == age) {
      if (s->stochastic) {
        tmp->next->number.i += number.i;
        s->size.i += number.i;
      } else {
        tmp->next->number.d += number.d;
        s->size.d += number.d;
      }
      return;
    }
  }
  tmp = (individual_data)malloc(sizeof(struct individual_st));
  tmp->age = age;
  if (s->stochastic) {
    tmp->number.i = number.i;
    s->size.i += tmp->number.i;
  } else {
    tmp->number.d = number.d;
    s->size.d += tmp->number.d;
  }
  tmp->next  = s->individual->next;
  s->individual->next = tmp;
}

void sdpop_print(sdpop s) {
  individual_data p = s->individual->next;
  unsigned int count;
  printf("/------------------>\nGamma mode: %d\nPropulation type: ",s->gamma_mode);
  if (s->stochastic) {
    printf("Stochastic\nPopulation size: %d\nDead: %d\nDeveloped: %d\nAge\tNumber\n",s->size.i,s->dead.i,s->developed.i);
    for (count=1; p; count++, p = p->next)
      printf("%d\t%d\n",p->age,p->number.i);
  } else {
    printf("Deterministic\nPopulation size: %g\nDead: %g\nDeveloped: %g\nAge\tNumber\n",s->size.d,s->dead.d,s->developed.d);
    for (count=1; p; count++, p = p->next)
      printf("%d\t%g\n",p->age,p->number.d);
  }
  printf("\\------------------>\n");
}

void sdpop_remove_next(sdpop s, individual_data p){
  if (!p) return;
  individual_data tmp = p->next;
  p->next = tmp->next;
  if (s->stochastic)
    s->size.i -= tmp->number.i;
  else
    s->size.d -= tmp->number.d;
  free(tmp);
}

void sdpop_empty(sdpop s) {
  if (s->stochastic) {
    while (s->size.i > 0 && s->individual->next)
      sdpop_remove_next(s,s->individual);
    s->size.i = 0;
    s->dead.i = 0;
    s->developed.i = 0;
  } else {
    while (s->size.d > DPOP_EPS && s->individual->next)
      sdpop_remove_next(s,s->individual);
    s->size.d = 0;
    s->dead.d = 0;
    s->developed.d = 0;
  }
}

void sdpop_destroy(sdpop *s) {
  sdpop_empty(*s);
  if ((*s)->individual->next)
    free((*s)->individual->next);
  free((*s)->individual);
  free((*s));
}

char calc_sdpop(sdpop s, individual_data tmpn, double prob, sdnum *k) {
  if (s->stochastic) {
    (*k).i = gsl_ran_binomial(RAND_GSL,
                              prob,
                              tmpn->number.i);
    tmpn->number.i -= (*k).i;
    s->size.i -= (*k).i;
    s->dead.i += (*k).i;
    if (tmpn->number.i == 0)
      return 1;
  } else {
    (*k).d = tmpn->number.d * prob;
    tmpn->number.d -= (*k).d;
    s->size.d -= (*k).d;
    s->dead.d += (*k).d;
    if (tmpn->number.d <= DPOP_EPS) {
      if (tmpn->number.d < 0)
        return -1;
      return 1;
    }
  }
  return 0;
}

void sdpop_kill(sdpop  s,
                double prob) { // probability of death
  sdnum k;
  char remove = -1;
  individual_data tmp, tmpn;
  if (s->stochastic)
    s->dead.i = 0;
  else
    s->dead.d = 0.0;
  for (tmp = s->individual;
       tmp && tmp->next;
       ) {
    tmpn = tmp->next;
    //
    remove = calc_sdpop(s, tmpn, prob, &k);
    if (s->stochastic)
      s->dead.i = k.i;
    else
      s->dead.d = k.d;
    //
    switch (remove) {
    default:
    case -1:
      printf("ERROR: Error in population!\n");
      sdpop_print(s);
      s->dead.d = NAN;
      break;
    case 0:
      tmp = tmp->next;
      break;
    case 1:
      tmp->next = tmpn->next;
      free(tmpn);
      if ((s->stochastic && s->size.i == 0) || (!(s->stochastic) && s->size.d <= DPOP_EPS)) break;
      remove = 0;
      break;
    }
  }
}

double calc_prob_gamma_raw(unsigned int age, double prob, double mean, double sd) {
  return gamma_dist_prob(mean,sd,age);
}

double calc_prob_gamma_hash(unsigned int age, double prob, double mean, double sd) {
  if (!gamma_dist_hash(mean,sd,age,&prob)) {
    printf("ERROR: Gamma distribution failed!\n");
    return NAN;
  }
  return prob;
}

double calc_prob_fixed(unsigned int age, double prob, double mean, double sd) {
  return age >= mean - 1.0;
}

double calc_prob_raw(unsigned int age, double prob, double mean, double sd) {
  return prob;
}

prob_func assign_prob_func (sdpop s, double prob, double mean, double sd) {
  if (mean == 0 && sd == 0 && prob >= 0 && prob <= 1) {
    // printf("Assigning calc_prob_raw\n");
    return calc_prob_raw;
  }
  //
  if (mean > 0 && sd > 0) {
    switch (s->gamma_mode) {
    case MODE_GAMMA_RAW:
      // printf("Assigning calc_prob_gamma_raw\n");
      return calc_prob_gamma_raw;
      break;
    case MODE_GAMMA_HASH:
      // printf("Assigning calc_prob_gamma_hash\n");
      return calc_prob_gamma_hash;
      break;
    default:
      printf("ERROR: Wrong distribution option selected: %d\n",s->gamma_mode);
      return 0;
      break;
    }
  }
  //
  if (mean > 0 && sd == 0) {
    // printf("Assigning calc_prob_fixed\n");
    return calc_prob_fixed;
  }
  //
  printf("ERROR: Wrong options for probability function!\n");
  return 0;
}

void sdpop_iterate(sdpop  s,
                   double dev_prob,   // fixed development probability
                   double dev_mean,   // mean development time
                   double dev_sd,     // sd development time
                   double death_prob, // fixed death probability
                   double death_mean, // mean time of death
                   double death_sd) { // sd time of death
  if (s->stochastic) {
    s->dead.i = 0;
    s->developed.i = 0;
  } else {
    s->dead.d = 0.0;
    s->developed.d = 0.0;
  }
  //
  prob_func calc_prob_death = assign_prob_func(s, death_prob, death_mean, death_sd);
  prob_func calc_prob_dev = assign_prob_func(s, dev_prob, dev_mean, dev_sd);
  //
  char remove = -1;
  sdnum k;
  double prob;
  individual_data tmp, tmpn;
  for (tmp = s->individual;
       tmp && tmp->next;
       ) {
    tmpn = tmp->next;
    // Death
    prob = calc_prob_death(tmpn->age, death_prob, death_mean, death_sd);
    remove = calc_sdpop(s, tmpn, prob, &k);
    if (s->stochastic)
      s->dead.i = k.i;
    else
      s->dead.d = k.d;
    //
    if (remove == -1) {
      printf("ERROR: Error in population!\n");
      sdpop_print(s);
      s->dead.d = NAN;
      return;
    }
    // Develop
    prob = calc_prob_dev(tmpn->age, dev_prob, dev_mean, dev_sd);
    remove = calc_sdpop(s, tmpn, prob, &k);
    //
    tmpn->age += 1;
    if (s->stochastic)
      s->developed.i = k.i;
    else
      s->developed.d = k.d;
    //
    if (tmpn->age >= MAX_DAYS) {
      printf("ERROR: Development time is too high!\n");
      sdpop_print(s);
      s->developed.d = NAN;
      return;
    }
    switch (remove) {
    default:
    case -1:
      printf("ERROR: Error in population!\n");
      sdpop_print(s);
      s->developed.d = NAN;
      break;
    case 0:
      tmp = tmp->next;
      break;
    case 1:
      tmp->next = tmpn->next;
      free(tmpn);
      if ((s->stochastic && s->size.i == 0) || (!(s->stochastic) && s->size.d <= DPOP_EPS)) break;
      remove = 0;
      break;
    }
  }
}
