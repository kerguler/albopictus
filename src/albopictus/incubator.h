#define EPS 1e-14

typedef struct
{
  double popsize;
  double popdev;
} incubator_data;

typedef struct incubator_st
{
  incubator_data data;
  struct incubator_st *next;
} *incubator;

void incubator_empty(incubator *);
void incubator_add(incubator *, double, double);
void incubator_remove(incubator *, double *);
void incubator_update(incubator, void (*)(double *, double *, double *), double *);
void incubator_develop_survive(incubator *, double, double, double, double, double, double *, double *, double *, char);
double incubator_sum(incubator);
void incubator_print(incubator);
