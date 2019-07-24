#ifndef GAMMA_H
#define GAMMA_H

#include <stdint.h>

#define EPS 1e-14

#define MODE_GAMMA_RAW    0
#define MODE_GAMMA_HASH   1
#define MODE_NBINOM_RAW   2
#define MODE_GAMMA_MATRIX 3
#define MODE_BINOM_RAW    4

#define n_MAX 400.0
#define n_STEP 1.0
#define n_i_MAX 400
#define row_SIZE 400
#define mean_MAX 201.0
#define mean_STEP 0.01
#define mean_i_MAX 20100
#define gamma_SIZE 8040000

#define MAX_DAYS 1000

#define gamma_matrix_sd 0.375

void set_gamma_mem(uint64_t);
void gamma_dist_destroy(void);
void gamma_dist_check(void);
double gamma_pdf(double, double, double);
char gamma_dist_hash(double, double, double, double *);
double gamma_dist_prob(double, double, double);
double gamma_dist_matrix(double, double);
void prepare_gamma_matrix(void);
double nbinom_prob(unsigned int, double, double);
double nbinom_dist_prob(double, double, unsigned int);

#endif
