#define EPS 1e-14

#define MODE_GAMMA_RAW   0
#define MODE_GAMMA_HASH  1
#define MODE_GAMMA_ALBO  2
#define MODE_BINOM_RAW   3

void gamma_dist_destroy(void);
void gamma_dist_check(void);
char gamma_dist_hash(double, double, double, double *);
double gamma_dist_prob(double, double, double);
double nbinom_prob(unsigned int, double, double);
