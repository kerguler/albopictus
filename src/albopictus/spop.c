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

spop spop_init(void) {
  spop pop = (spop)malloc(sizeof(struct population_st));
  people_data tmp = (people_data)malloc(sizeof(struct people_st));
  tmp->development = 0;
  tmp->batchsize = 0;
  tmp->next = 0;
  pop->people = tmp;
  pop->popsize = 0;
  return pop;
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
  if ((*s)->people->next)
    free((*s)->people->next);
  free((*s)->people);
  free((*s));
}

void spop_print(spop s) {
  people_data p = s->people->next;
  int count;
  printf("/------------------>\n");
  printf("Population size: %d\n",s->popsize);
  for (count=1; p; count++) {
    printf("%d\t%d\t%d\n",count,p->development,p->batchsize);
    p = p->next;
  }
  printf("\\------------------>\n");
}

int spop_survive(spop   s,
                 double d_mean, // development time
                 double d_sd, // (if <= 0, fixed development, -d_mean)
                 double p_mean, // survival time
                 double p_sd, // (if <=0, fixed survival, -p_mean)
                 char   mode) {
  int ret = 0;
  char remove = 0;
  int k;
  double prob;
  people_data tmp, tmpn;
  for (tmp = s->people;
       tmp && tmp->next;
       ) {
    tmpn = tmp->next;
    // Survive
    if (p_sd > 0) { // survival time provided
      switch (mode) { // probability of "death" happenning
      case MODE_GAMMA_RAW:
        prob = gamma_dist_prob(p_mean,p_sd,tmpn->development);
        break;
      case MODE_GAMMA_HASH:
        if (!gamma_dist_hash(p_mean,p_sd,tmpn->development,&prob)) {
          printf("ERROR: Gamma distribution failed!\n");
          exit(1);
        }
        break;
      default:
        printf("ERROR: Wrong distribution option selected: %d\n",mode);
        break;
      }
    } else { // fixed daily probability of death
      prob = -p_mean;
    }
    k = gsl_ran_binomial(RAND_GSL,
                         prob,
                         tmpn->batchsize);
    //
    // printf("%d will die from %g,%d,%d\n",k,prob,tmpn->development,tmpn->batchsize);
    //    
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
      if (d_sd > 0) { // development time provided
        switch (mode) { // probability of "development" happenning
        case MODE_GAMMA_RAW:
          prob = gamma_dist_prob(d_mean,d_sd,tmpn->development);
          break;
        case MODE_GAMMA_HASH:
          if (!gamma_dist_hash(d_mean,d_sd,tmpn->development,&prob)) {
            printf("ERROR: Gamma distribution failed!\n");
            exit(1);
          }
          break;
        default:
          printf("ERROR: Wrong distribution option selected: %d\n",mode);
          break;
        }
      } else { // fixed daily development probability
        prob = -d_mean;
      }
      k = gsl_ran_binomial(RAND_GSL,
                           prob,
                           tmpn->batchsize);
      // printf("%d will develop with %g probability in %d,%d\n",k,prob,tmpn->development,tmpn->batchsize);
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
