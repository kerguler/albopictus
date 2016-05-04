#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "ran_gen.h"

//Random number variables
const gsl_rng_type *R_TYPE;
gsl_rng *RAND_GSL;

void rng_setup()
//Can also run like:
//GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./GDirect
{
  fprintf(stderr,"Setting up RNG\n");
  FILE *dev;

  unsigned int i;

  dev=fopen("/dev/urandom","r");

  fread(&i,sizeof(unsigned int),1,dev);

  fclose(dev);

  gsl_rng_env_setup();

  R_TYPE = gsl_rng_mt19937;
  RAND_GSL = gsl_rng_alloc(R_TYPE);
  gsl_rng_set(RAND_GSL,i);
}

void rng_setup_seed(unsigned int seed)
//Can also run like:
//GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./GDirect
{
  fprintf(stderr,"Setting up RNG with seed = %u\n",seed);

  gsl_rng_env_setup();

  R_TYPE = gsl_rng_mt19937;
  RAND_GSL = gsl_rng_alloc(R_TYPE);
  gsl_rng_set(RAND_GSL,seed);
}

void rng_destroy()
{
  gsl_rng_free(RAND_GSL);

  fprintf(stderr,"RNG destroyed\n");
}

double rng_exponential(const double mu)
{
  double u = gsl_rng_uniform(RAND_GSL);
  return -(1.0 / mu) * (u == 0.0 ? 0.0 : log (u));
}
