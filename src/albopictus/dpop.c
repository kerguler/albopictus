#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include "gamma.h"
#include "dpop.h"


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

dpop dpop_init(void) {
  dpop pop = (dpop)malloc(sizeof(struct population_det_st));
  people_det_data tmp = (people_det_data)malloc(sizeof(struct people_det_st));
  tmp->development = 0;
  tmp->batchsize = 0;
  tmp->next = 0;
  pop->people = tmp;
  pop->popsize = 0;
  return pop;
}

void dpop_add(dpop s, int dev, double size){
  people_det_data tmp;
  for (tmp = s->people;
       tmp && tmp->next;
       tmp = tmp->next) {
    if (tmp->next->development == dev) {
      tmp->next->batchsize += size;
      s->popsize += size;
      return;
    }
  }
  tmp = (people_det_data)malloc(sizeof(struct people_det_st));
  tmp->development = dev;
  tmp->batchsize = size;
  tmp->next  = s->people->next;
  s->people->next = tmp;
  s->popsize += tmp->batchsize;
}

void dpop_remove_next(dpop s, people_det_data p){
  if (!p) return;
  people_det_data tmp = p->next;
  p->next = tmp->next;
  s->popsize -= tmp->batchsize;
  free(tmp);
}

void dpop_empty(dpop s) {
  while (s->popsize > DPOP_EPS && s->people->next) {
    dpop_remove_next(s,s->people);
  }
  s->popsize = 0;
}

void spop_destroy(dpop *s) {
  dpop_empty(*s);
  if ((*s)->people->next)
    free((*s)->people->next);
  free((*s)->people);
  free((*s));
}

void dpop_print(dpop s) {
  people_det_data p = s->people->next;
  int count;
  printf("/------------------>\n");
  printf("Population size: %g\n",s->popsize);
  for (count=1; p; count++) {
    printf("%d\t%d\t%g\n",count,p->development,p->batchsize);
    p = p->next;
  }
  printf("\\------------------>\n");
}

double dpop_kill(dpop   s,
                 double prob) {  // fraction of death
  double ret = 0;
  char remove = 0;
  doubble k;
  people_det_data tmp, tmpn;
  for (tmp = s->people;
       tmp && tmp->next;
       ) {
    tmpn = tmp->next;
    // Survive
    k = tmpn->batchsize * prob;
    //
    tmpn->batchsize -= k;
    s->popsize -= k;
    ret += k;
    if (tmpn->batchsize <= DPOP_EPS) {
      remove = 1;
      if (tmpn->batchsize < 0) {
        printf("ERROR: Error in population!\n");
        dpop_print(s);
        return nan;
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

double dpop_survive(dpop   s,
                    double d_mean, // development time
                    double d_sd, // (if <= 0, fixed development, -d_mean)
                    double p_mean, // survival time
                    double p_sd, // (if <=0, fixed survival, -p_mean)
                    char   mode) {
  double ret = 0;
  char remove = 0;
  double k;
  double prob;
  people_det_data tmp, tmpn;
  for (tmp = s->people;
       tmp && tmp->next;
       ) {
    tmpn = tmp->next;
    // Survive
    if (p_sd > 0) { // survival time provided
      switch (mode) { // fraction of "death" happenning
      case MODE_GAMMA_RAW:
        prob = gamma_dist_prob(p_mean,p_sd,tmpn->development);
        break;
      case MODE_GAMMA_HASH:
        if (!gamma_dist_hash(p_mean,p_sd,tmpn->development,&prob)) {
          printf("ERROR: Gamma distribution failed!\n");
          return nan;
        }
        break;
      default:
        printf("ERROR: Wrong distribution option selected: %d\n",mode);
        break;
      }
    } else { // fixed daily probability of death
      prob = -p_mean;
    }
    k = tmpn->batchsize * prob;
    //    
    tmpn->batchsize -= k;
    s->popsize -= k;
    if (tmpn->batchsize <= DPOP_EPS) {
      remove = 1;
    } else if (tmpn->batchsize < 0) {
      printf("ERROR: Error in population!\n");
      dpop_print(s);
      return nan;
    } else {
      // Develop
      if (d_sd > 0) { // development time provided
        switch (mode) { // fraction of population developing
        case MODE_GAMMA_RAW:
          prob = gamma_dist_prob(d_mean,d_sd,tmpn->development);
          break;
        case MODE_GAMMA_HASH:
          if (!gamma_dist_hash(d_mean,d_sd,tmpn->development,&prob)) {
            printf("ERROR: Gamma distribution failed!\n");
            return nan;
          }
          break;
        default:
          printf("ERROR: Wrong distribution option selected: %d\n",mode);
          break;
        }
      } else { // fixed daily development probability
        prob = -d_mean;
      }
      k = tmpn->batchsize * prob;
      //
      tmpn->development += 1;
      tmpn->batchsize -= k;
      s->popsize -= k;
      ret += k;
      if (tmpn->batchsize <= DPOP_EPS) {
        remove = 1;
        if (tmpn->batchsize<0) {
          printf("ERROR: Error in population!\n");
          dpop_print(s);
          return nan;
        }
      }
      if (tmpn->development >= MAX_DAYS) {
        printf("ERROR: Development time is too high!\n");
        dpop_print(s);
        return nan;
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
