#include "Python.h"
#include <math.h>
#include <time.h>
#include "gamma.h"
#include "spop.h"

volatile clock_t start = 0, diff;
double time2here(void) {
  diff = clock() - start;
  start = clock();
  return diff;
}

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

#define NumParAea      48
#define NumMetAea      5

// --------------------------------------------
// Gamma distribution
// --------------------------------------------

char gamma_mode = MODE_GAMMA_HASH;
void set_gamma_mode(char mode) {
  gamma_mode = mode;
}

// --------------------------------------------
// modelDelayAalbopictus
// --------------------------------------------

#define alpha_init_n1     0
#define alpha_init_n2     1
#define alpha_init_n3     2
#define alpha_init_n4f    3
#define alpha_init_nBS    4

#define alpha_F4_1        5
#define alpha_F4_2        6
#define alpha_F4_3        7

#define alpha_n1_death_1  8
#define alpha_n1_death_2  9
#define alpha_n1_death_3  10

#define alpha_n2_death_1  11
#define alpha_n2_death_2  12
#define alpha_n2_death_3  13

#define alpha_n3_death_1  14
#define alpha_n3_death_2  15
#define alpha_n3_death_3  16

#define alpha_n1_mean_1   17
#define alpha_n1_mean_2   18
#define alpha_n1_mean_3   19
#define alpha_n1_std_1    20
#define alpha_n1_std_2    21
#define alpha_n1_std_3    22

#define alpha_n2_mean_1   23
#define alpha_n2_mean_2   24
#define alpha_n2_mean_3   25
#define alpha_n2_std_1    26
#define alpha_n2_std_2    27
#define alpha_n2_std_3    28

#define alpha_n3_mean_1   29
#define alpha_n3_mean_2   30
#define alpha_n3_mean_3   31
#define alpha_n3_std_1    32
#define alpha_n3_std_2    33
#define alpha_n3_std_3    34

#define alpha_n4_mean_1   35
#define alpha_n4_mean_2   36
#define alpha_n4_mean_3   37
#define alpha_n4_std_1    38
#define alpha_n4_std_2    39
#define alpha_n4_std_3    40

#define alpha_deltaT      41
#define alpha_avoidT      42

#define alpha_BS_dprec    43
#define alpha_BS_nevap    44
#define alpha_BS_const    45

#define alpha_pdens_1     46
#define alpha_pdens_2     47

#define f_egg(T,x1,x2,x3) (max(0.0, (x1) + (x2)*(T) + (x3)*pow((T),2.0)))
#define f_death(T,x1,x2,x3) (max(0.0, min(1.0, (x1) + (x2)*(T) + (x3)*pow((T),2.0))))
#define f_dev(T,x1,x2,x3) (max(0.0, (x1) + (x2)*(T) + (x3)*pow((T),2.0)))
#define f_exp(T,x1,x2,x3) (max(0.0, (x1)*exp((x2)+(x3)*(T))))
#define f_dsigi(T,x1,x2,x3) (max(0.0, min(1.0, 1.0-(x1)/((1.0+exp((x2)-(T)))*(1.0+exp((T)-(x3)))) )))

//double timebefore = 0;
//double timeof = 0;
//double timeafter = 0;

void calculate(double *air_temp,
               double *precipitation,
               double *nevaporation,
               double *param,
               spop   *conn1,
               spop   *conn2,
               spop   *conn3,
               spop   *conn4,
               double *n1,
               double *n2,
               double *n3,
               double *n4f,
               double *nBS,
               int    TIME) {
  // ---------------------
  // modelDelayAalbopictus
  // ---------------------

  double deltaT = param[alpha_deltaT];

  double Ta = air_temp[TIME];

  double Tw = Ta + deltaT;

  /*
   * Update the number of breeding sites
   */
  double eva = param[alpha_BS_nevap] * nevaporation[TIME];
  double revap = eva ? exp(-1.0 / eva) : 0.0;
  (*nBS) = param[alpha_BS_const] + param[alpha_BS_dprec] * precipitation[TIME] + revap * (*nBS);
  // (*K) = revap == 1.0 ? (*nBS) / (TIME+1.0) : (*nBS) * (revap-1.0) / (pow(revap, (TIME+1))-1.0);

  /*
   * Update survival and development rates/probabilities according to
   * the new environmental conditions and breeding sites
   */

  // Density-dependent decrease of survival with unscaled carrying capacity
  double fdens = ((*n1) + (*n2) + (*n3)) / (*nBS);
  double pdens = exp(-param[alpha_pdens_1]*pow(fdens,param[alpha_pdens_2]));

  // Fecundity
  double bigF4 = (0.25)*f_egg(Ta,param[alpha_F4_1],param[alpha_F4_2],param[alpha_F4_3]);

  // Egg mortality
  double p1_Tw = 1.0 - pdens*(pow(1.0-f_death(Tw,param[alpha_n1_death_1],param[alpha_n1_death_2],param[alpha_n1_death_3]),0.25));
  // Larva mortality
  double p2_Tw = 1.0 - pdens*(pow(1.0-f_dsigi(Tw,param[alpha_n2_death_1],param[alpha_n2_death_2],param[alpha_n2_death_3]),0.25));
  // Pupa mortality
  double p3_Tw = 1.0 - pdens*(pow(1.0-f_dsigi(Tw,param[alpha_n3_death_1],param[alpha_n3_death_2],param[alpha_n3_death_3]),0.25));

  // Egg development time
  double d1mean = (4.0)*f_dev(Tw,param[alpha_n1_mean_1],param[alpha_n1_mean_2],param[alpha_n1_mean_3]);
  double d1sd   = (4.0)*f_dev(Tw,param[alpha_n1_std_1], param[alpha_n1_std_2], param[alpha_n1_std_3]);
  // Larva development time
  double d2mean = (4.0)*f_dev(Tw,param[alpha_n2_mean_1],param[alpha_n2_mean_2],param[alpha_n2_mean_3]);
  double d2sd   = (4.0)*f_dev(Tw,param[alpha_n2_std_1], param[alpha_n2_std_2], param[alpha_n2_std_3]);
  // Pupa development time
  double d3mean = (4.0)*f_dev(Tw,param[alpha_n3_mean_1],param[alpha_n3_mean_2],param[alpha_n3_mean_3]);
  double d3sd   = (4.0)*f_dev(Tw,param[alpha_n3_std_1], param[alpha_n3_std_2], param[alpha_n3_std_3]);

  // Adult lifetime (from emergence)
  double Teff = Ta > param[alpha_avoidT] ? param[alpha_avoidT] : Ta;
  double p4mean = (4.0)*f_exp(Teff,param[alpha_n4_mean_1],param[alpha_n4_mean_2],param[alpha_n4_mean_3]);
  double p4sd   = (4.0)*f_exp(Teff,param[alpha_n4_std_1], param[alpha_n4_std_2], param[alpha_n4_std_3]);

  // printf("%g,%g %g %g,%g,%g %g,%g %g,%g %g,%g %g,%g\n",Ta,Tw,bigF4,p1_Tw,p2_Tw,p3_Tw,d1mean,d1sd,d2mean,d2sd,d3mean,d3sd,p4mean,p4sd);

  /*
   * Survival first, and then, development.
   * Update development times and survival proportions
   * of the immature stages and juvenile adults
   */
  //
  // Eggs
  spop_iterate((*conn1),
               0,
               d1mean, d1sd, // development (gamma-distributed)
               0,
               p1_Tw, // death (fixed rate)
               0, 0,
               0,
               0);
  //
  // Larvae
  spop_iterate((*conn2),
               0,
               d2mean, d2sd, // development (gamma-distributed)
               0,
               p2_Tw, // death (fixed rate)
               0, 0,
               0,
               0);
  //
  // Pupae
  spop_iterate((*conn3),
               0,
               d3mean, d3sd, // development (gamma-distributed)
               0,
               p3_Tw, // death (fixed rate)
               0, 0,
               0,
               0);
  //
  // Adult females
  spop_iterate((*conn4),
               0,
               0, 0, // development (no development)
               0,
               0,
               p4mean, p4sd, // death (gamma-distributed)
               0,
               0);
  //
  // Adults survived to produce eggs
  double n4_reproduce = (*conn4)->size.d;
  //
  // Perform state transformations
  spop_add((*conn4), 0, 0, 0, 0.5 * (*conn3)->developed.d);
  spop_add((*conn3), 0, 0, 0, (*conn2)->developed.d);
  spop_add((*conn2), 0, 0, 0, (*conn1)->developed.d);
  //
  // Lay eggs
  double neweggs = bigF4*n4_reproduce; // Total number of eggs that will be laid that day
  //
  spop_add((*conn1),0,0,0,neweggs); // Eggs
  //
  // Update all stages and state variables
  (*n1) = (*conn1)->size.d;
  (*n2) = (*conn2)->size.d;
  (*n3) = (*conn3)->size.d;
  (*n4f) = (*conn4)->size.d;
}

// --------------------------------------------
// sim_model
// --------------------------------------------

void numparModel(int *np, int *nm) {
  *np = NumParAea;
  *nm = NumMetAea;
}

void param_model(char **names, double *param) {
  char temp[NumMetAea+NumParAea][256] = {
    "coln1","coln2","coln3","coln4f","colnBS",
    "init_n1",
    "init_n2",
    "init_n3",
    "init_n4f",
    "init_nBS",
    "F4_1",
    "F4_2",
    "F4_3",
    "n1_death_1",
    "n1_death_2",
    "n1_death_3",
    "n2_death_1",
    "n2_death_2",
    "n2_death_3",
    "n3_death_1",
    "n3_death_2",
    "n3_death_3",
    "n1_mean_1",
    "n1_mean_2",
    "n1_mean_3",
    "n1_std_1",
    "n1_std_2",
    "n1_std_3",
    "n2_mean_1",
    "n2_mean_2",
    "n2_mean_3",
    "n2_std_1",
    "n2_std_2",
    "n2_std_3",
    "n3_mean_1",
    "n3_mean_2",
    "n3_mean_3",
    "n3_std_1",
    "n3_std_2",
    "n3_std_3",
    "n4_mean_1",
    "n4_mean_2",
    "n4_mean_3",
    "n4_std_1",
    "n4_std_2",
    "n4_std_3",
    "deltaT",
    "avoidT",
    "BS_dprec",
    "BS_nevap",
    "BS_const",
    "pdens_1",
    "pdens_2"
  };
  int i;
  for (i=0; i<(NumMetAea+NumParAea); i++)
    names[i] = strdup(temp[i]);
  // ---
  param[alpha_init_n1]    = 1000.0;
  param[alpha_init_n2]    = 1000.0;
  param[alpha_init_n3]    = 1000.0;
  param[alpha_init_n4f]   = 1000.0;
  param[alpha_init_nBS]   = 1000.0;

  param[alpha_F4_1]       = -25.49944116;
  param[alpha_F4_2]       = 2.64761782;
  param[alpha_F4_3]       = -0.05736217;

  param[alpha_n1_death_1] = 1.057782972;
  param[alpha_n1_death_2] = -0.117725165;
  param[alpha_n1_death_3] = 0.003552023;

  param[alpha_n2_death_1] = 0.9878398;
  param[alpha_n2_death_2] = 10.6100650;
  param[alpha_n2_death_3] = 33.2006971;

  param[alpha_n3_death_1] = 0.9681704;
  param[alpha_n3_death_2] = 11.7116140;
  param[alpha_n3_death_3] = 32.3451635;

  param[alpha_n1_mean_1]  = 27.11593866;
  param[alpha_n1_mean_2]  = -1.84221694;
  param[alpha_n1_mean_3]  = 0.03242374;
  param[alpha_n1_std_1]   = 2.012053822;
  param[alpha_n1_std_2]   = -0.139289973;
  param[alpha_n1_std_3]   = 0.002593632;

  param[alpha_n2_mean_1]  = 126.4857775;
  param[alpha_n2_mean_2]  = -7.6294751;
  param[alpha_n2_mean_3]  = 0.1290408;
  param[alpha_n2_std_1]   = 27.46596873;
  param[alpha_n2_std_2]   = -1.50113421;
  param[alpha_n2_std_3]   = 0.02370162;

  param[alpha_n3_mean_1]  = 19.41735844;
  param[alpha_n3_mean_2]  = -1.25195579;
  param[alpha_n3_mean_3]  = 0.02097719;
  param[alpha_n3_std_1]   = 0.0;
  param[alpha_n3_std_2]   = 0.0;
  param[alpha_n3_std_3]   = 0.0;

  param[alpha_n4_mean_1]  = 2.0548142;
  param[alpha_n4_mean_2]  = 5.9171821;
  param[alpha_n4_mean_3]  = -0.1165213;
  param[alpha_n4_std_1]   = 2.7238169;
  param[alpha_n4_std_2]   = 4.8608867;
  param[alpha_n4_std_3]   = -0.1049982;

  param[alpha_deltaT]     = 0.0;
  param[alpha_avoidT]     = 100.0;

  param[alpha_BS_dprec]   = 100.0;
  param[alpha_BS_nevap]   = 100.0;
  param[alpha_BS_const]   = 0.0;

  param[alpha_pdens_1]    = 0.25;
  param[alpha_pdens_2]    = 4.0;
}

void sim_model(double               *envar,
               double               *param,
               int                 *finalT,
               int                *control,
               double              *result,
               int                *success) {
  double *air_temp             = envar + 0*(*finalT);
  double *precipitation        = envar + 1*(*finalT);
  double *nevaporation         = envar + 2*(*finalT);

  double *controlpar = 0;
  if ((*control)) {
    controlpar = param + NumParAea;
  }
	       
  double *colT    = result + 0*(*finalT);
  double *coln1   = result + 1*(*finalT);
  double *coln2   = result + 2*(*finalT);
  double *coln3   = result + 3*(*finalT);
  double *coln4f  = result + 4*(*finalT);
  double *colnBS  = result + 5*(*finalT);

  int TIME = 0;
  (*success) = 2;
  // Set initial conditions
  double n1 = param[alpha_init_n1];
  double n2 = param[alpha_init_n2];
  double n3 = param[alpha_init_n3];
  double n4f = param[alpha_init_n4f];
  double nBS = param[alpha_init_nBS];
  //
  spop conn1 = spop_init(0,gamma_mode);
  spop conn2 = spop_init(0,gamma_mode);
  spop conn3 = spop_init(0,gamma_mode);
  spop conn4 = spop_init(0,gamma_mode);
  if (n1 > DPOP_EPS) spop_add(conn1,0,0,0,n1);
  if (n2 > DPOP_EPS) spop_add(conn2,0,0,0,n2);
  if (n3 > DPOP_EPS) spop_add(conn3,0,0,0,n3);
  if (n4f > DPOP_EPS) spop_add(conn4,0,0,0,n4f);
  // Record state
  colT[TIME] = TIME;
  coln1[TIME] = n1;
  coln2[TIME] = n2;
  coln3[TIME] = n3;
  coln4f[TIME] = n4f;
  colnBS[TIME] = nBS;
  //
  for (TIME=1; TIME<(*finalT); TIME++) {
    // Take a step
    calculate(air_temp,
              precipitation,
              nevaporation,
              param,
              &conn1,
              &conn2,
              &conn3,
              &conn4,
              &n1,
              &n2,
              &n3,
              &n4f,
              &nBS,
              TIME);
    //
    // Record state
    colT[TIME] = TIME;
    coln1[TIME] = n1;
    coln2[TIME] = n2;
    coln3[TIME] = n3;
    coln4f[TIME] = n4f;
    colnBS[TIME] = nBS;
    if (isnan(n1) ||
        isnan(n2) ||
        isnan(n3) ||
        isnan(n4f) ||
        isnan(nBS)) {
      //
      printf("ERROR_NAN: %d,%g,%g,%g,%g,%g\n",TIME,n1,n2,n3,n4f,nBS);
      printf("ERROR_PAR: ");
      int i;
      for (i=0; i<NumParAea; i++)
        printf("%s%g",i==0?"":",",param[i]);
      printf("\n");
      //
      (*success) = 0;
      goto endall;
    }
  }
  (*success) = 1;
  //
 endall:
  //
  spop_destroy(&conn1);
  spop_destroy(&conn2);
  spop_destroy(&conn3);
  spop_destroy(&conn4);
  //
  gamma_dist_check();
  //
  //printf("Times: %g %g %g\n",timebefore,timeof,timeafter);
}

// ---------------------------------------------------------------------------


