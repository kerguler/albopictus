#define EPS 1e-14

#define gamma_sd 0.375

void gamma_mean_destroy(void);
char gamma_mean_hash(double, double, double, double *);
double gamma_mean_prob(double, double, double);
double nbinom_prob(unsigned int, double, double);
