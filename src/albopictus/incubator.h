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

void test();

void prepare_gamma();
double prob_gamma(double n, double k, double theta);
double prob_nbinom(unsigned int k, double p, double n);
double prob_gamma_matrix(double n, double mean);
void incubator_empty(incubator *s);
void incubator_add(incubator *s, double popsize, double popdev);
void incubator_remove(incubator *s, double *val);
void incubator_update(incubator s, void (*fun)(double *size, double *dev, double *par), double *par);
void incubator_develop(incubator *s, double lambda, double sd, double *val, double thr, double *valthr, char mode);
double incubator_sum(incubator s);
void incubator_print(incubator s);
