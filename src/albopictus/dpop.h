#ifndef DPOP_H
#define DPOP_H

#define DPOP_EPS 1e-14

typedef struct people_det_st
{
  int development;
  double batchsize;
  struct people_st *next;
} *people_det_data;

typedef struct population_det_st
{
  people_det_data people;
  double popsize;
} *dpop;

dpop dpop_init(void);
void dpop_add(dpop, int, double);
void dpop_remove_next(dpop, people_det_data);
void dpop_empty(dpop);
void dpop_destroy(dpop *);
void dpop_print(dpop);
double dpop_kill(dpop, double);
double dpop_survive(dpop, double, double, double, double, char);

#endif
