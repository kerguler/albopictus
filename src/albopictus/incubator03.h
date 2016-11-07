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

void incubator_empty(incubator *s);
void incubator_add(incubator *s, double popsize, double popdev);
void incubator_remove(incubator *s, double *val);
void incubator_update(incubator s, void (*fun)(double *size, double *dev, double *par), double *par);
double incubator_sum(incubator s);
void incubator_print(incubator s);
