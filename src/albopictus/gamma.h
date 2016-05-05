#define EPS 1e-14

void gamma_mean_destroy(void);
void gamma_mean_trim(void);
char gamma_mean_hash(double, double, double, double *);
double gamma_mean_prob(double, double, double);
double nbinom_prob(unsigned int, double, double);
