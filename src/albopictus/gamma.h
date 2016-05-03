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

void test_gamma_matrix(void);
void prepare_gamma(void);
double prob_gamma(double, double, double);
double prob_nbinom(unsigned int, double, double);
double prob_gamma_matrix(double, double);
