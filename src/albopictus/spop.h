typedef struct people_st
{
  int development;
  int batchsize;
  struct people_st *next;
} *people_data;

typedef struct population_st
{
  people_data people;
  double development_time;
  double development_time_sd;
  double survival;
  int popsize;
  double *prdev;
} *spop;

double prob_gamma(double n, double k, double theta);
spop spop_init(double surv, double dev_time_mean, double dev_time_sd);
void spop_setdevtime(spop s, double dev_time_mean, double dev_time_sd);
void spop_setsurv(spop s, double surv);
void spop_add(spop s, int dev, int size);
void spop_remove_next(spop s, people_data p);
void spop_empty(spop s);
void spop_destroy(spop *s);
void spop_print(spop s);
int spop_survive(spop s);
void spop_swap(spop pfrom, spop pto);
