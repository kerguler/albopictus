#define EPS 1e-14

#define n_MAX 400.0
#define n_STEP 1.0
#define n_i_MAX 400
#define row_SIZE 400
#define mean_MAX 201.0
#define mean_STEP 0.01
#define mean_i_MAX 20100
#define gamma_SIZE 8040000

#define gamma_sd 0.375

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

void test_gamma_matrix(void);
void prepare_gamma(void);
double prob_gamma(double, double, double);
double prob_nbinom(unsigned int, double, double);
double prob_gamma_matrix(double, double);
void incubator_empty(incubator *);
void incubator_add(incubator *, double, double);
void incubator_remove(incubator *, double *);
void incubator_update(incubator, void (*)(double *, double *, double *), double *);
void incubator_develop_survive(incubator *, double, double, double, double, double, double *, double *, double *, char);
double incubator_sum(incubator);
void incubator_print(incubator);
