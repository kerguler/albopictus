#ifdef R
#define modelHeader
#include "R.h"
#include "Rmath.h"
#else
#include "Python.h"
#include <math.h>
#endif
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "gamma.h"
#include "spop.h"
#include "ran_gen.h"

volatile clock_t start = 0, diff;
double time2here(void) {
  diff = clock() - start;
  start = clock();
  return diff;
}

#define ZERO    0
#define MAX_DEV 1000
#define CHECK(x) ((x)<0)

extern gsl_rng *RAND_GSL;

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

#define NumParAea      45
#define NumMetAea      12

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

#define alpha_p1_1        0
#define alpha_p1_2        1
#define alpha_p1_3        2
#define alpha_p2_1        3
#define alpha_p2_2        4
#define alpha_p2_3        5
#define alpha_p3_1        6
#define alpha_p3_2        7
#define alpha_p3_3        8
#define alpha_d4_1        9
#define alpha_d4_2        10
#define alpha_d4_3        11
#define alpha_F4_1        12
#define alpha_F4_2        13
#define alpha_F4_3        14
#define alpha_d1_1        15
#define alpha_d1_2        16
#define alpha_d1_3        17
#define alpha_d2_1        18
#define alpha_d2_2        19
#define alpha_d2_3        20
#define alpha_d3_1        21
#define alpha_d3_2        22
#define alpha_d3_3        23

#define alpha_n23_surv    24

#define alpha_deltaT      25

#define alpha_BS_pdens    26
#define alpha_BS_dprec    27
#define alpha_BS_nevap    28
#define alpha_BS_Klulc    29

#define alpha_tbm_1       30
#define alpha_tbm_2       31
#define alpha_tbm_3       32

#define alpha_gtc_1       33
#define alpha_gtc_2       34
#define alpha_gtc_3       35

#define alpha_n23_1       36
#define alpha_n23_2       37
#define alpha_n23_3       38
#define alpha_n23_4       39
#define alpha_n23_5       40

#define alpha_capture     41

#define alpha_init_egg    42

#define alpha_flush_thr   43
#define alpha_flush_surv  44

#define flin(x, a1, a2) (max(0.0, min(1.0, (a1) + (a2)*(x))))
#define dsig(x, a1, a2, a3) (max(0.0, min(1.0, (a1)/((1.0+exp((a2)-(x)))*(1.0+exp((x)-(a3)))))))
#define dsig2(x, a1, a2, a3) (max(0.0, min(100.0, (a1)/((1.0+exp((a2)-(x)))*(1.0+exp((x)-(a3)))))))
#define poly(x, a1, a2, a3) (max(0.0, (a1) + (a2)*(x) + (a3)*pow((x),2)))
#define expd(x, k) (exp(-(k)*(x)))
#define fpow(x, a1, a2) (min(1.0, (a1)*pow((x),(a2))))
#define Tthr(T,Tmn,Tstd) (1.0/(1.0+exp(((Tmn)-(T))/(Tstd))))
#define pdiap(P,Pmn,Pstd,tt) (1.0-(tt)-(1.0-(tt))/((1.0+exp(((Pmn)-(P))/(Pstd)))))

//double timebefore = 0;
//double timeof = 0;
//double timeafter = 0;

double prob_capture = 0;

void set_param(double *param) {
  prob_capture = pow(10,param[alpha_capture]);
}

char calculate(double *mean_air_temp,
               double *daily_precipitation,
               double *popdens,
               double *lulc,
               double *daily_capture,
               double *param,
               spop   *conn1,
               spop   *conn2,
               spop   *conn3,
               spop   *conn4j,
               spop   *conn4,
               int    *n1,
               int    *n2,
               int    *n3,
               int    *n4fj,
               int    *n4f,
               double *nBS,
               double *K,
               int    *d4,
               int    *d4s,
               int    *F4,
               int    *egg,
               int    *cap,
               int    TIME) {
  // ---------------------
  // modelDelayAalbopictus
  // ---------------------

  double deltaT = param[alpha_deltaT];

  double Ta = mean_air_temp[TIME];
  
  double Tw = Ta + deltaT;

  /*
   * TO DO: Check how this can be implemented with a stochastic dynamics
   * Land use incorporated (the resolution is 0.01 - modis dataset)
   *
   * Update the number of breeding sites
   */
  (*nBS) = param[alpha_BS_pdens] * (*popdens) / (param[alpha_BS_Klulc] + (*lulc)) +
    param[alpha_BS_dprec] * daily_precipitation[TIME] +
    param[alpha_BS_nevap] * (*nBS);
  (*K) = param[alpha_BS_nevap] == 1.0 ? (*nBS) / (TIME+1.0) : (*nBS) * (param[alpha_BS_nevap]-1.0) / (pow(param[alpha_BS_nevap], (TIME+1))-1.0);

  /*
   * Update survival and development rates/probabilities according to
   * the new environmental conditions and breeding sites
   */
  // Density of the immature stages
  // Assume uniform distribution across all breeding sites
  double n23dens = (double)((*n2) + (*n3)) / (double)(*K);
  double densd = expd(n23dens,param[alpha_n23_surv]);
  double densdev = fpow(n23dens,param[alpha_n23_1],param[alpha_n23_2]*poly(Tw,param[alpha_n23_3],param[alpha_n23_4],param[alpha_n23_5]));
  // Fecundity (per adult female)
  // WARNING: Existing parameters are fixed for daily egg laying
  //          The following is aimed at egg laying per gonotrophic cycle
  //          Experimental data need to be updated
  double bigF4 = poly(Ta,param[alpha_F4_1],param[alpha_F4_2],param[alpha_F4_3]);
  //
  // Excess rain flushes the aquatic life stages
  double flush_surv = daily_precipitation[TIME] > param[alpha_flush_thr] ? param[alpha_flush_surv] : 1.0;
  //
  // Egg mortality
  double p1_Tw = 1.0 - dsig(Tw,param[alpha_p1_1],param[alpha_p1_2],param[alpha_p1_3])*flush_surv;
  // Larval mortality
  double p2_Tw = 1.0 - dsig(Tw,param[alpha_p2_1],param[alpha_p2_2],param[alpha_p2_3])*densd*flush_surv;
  // Pupal mortality
  double p3_Tw = 1.0 - dsig(Tw,param[alpha_p3_1],param[alpha_p3_2],param[alpha_p3_3])*densd*flush_surv;
  
  // Egg development time
  double d1 = max(1.0, poly(Tw,param[alpha_d1_1],param[alpha_d1_2],param[alpha_d1_3]));
  // Larval development time
  double d2 = max(1.0, poly(Tw,param[alpha_d2_1],param[alpha_d2_2],param[alpha_d2_3])*densdev);
  // Pupal development time
  double d3 = max(1.0, poly(Tw,param[alpha_d3_1],param[alpha_d3_2],param[alpha_d3_3])*densdev);
  // Time to first blood meal
  double alpha_blood = poly(Ta,param[alpha_tbm_1],param[alpha_tbm_2],param[alpha_tbm_3]);
  // Gonotrophic cycle length
  double alpha_ovipos = poly(Ta,param[alpha_gtc_1],param[alpha_gtc_2],param[alpha_gtc_3]);
  // Adult lifetime (from emergence)
  double dd4 = dsig2(Ta,param[alpha_d4_1],param[alpha_d4_2],param[alpha_d4_3]);
  double dd4s = gamma_matrix_sd*dd4;

  /*
   * Survival first, and then, development.
   * Update development times and survival proportions
   * of the immature stages and juvenile adults
   */
  //
  char test = 0;
  // Normal eggs
  test = spop_iterate((*conn1),
                      0,
                      d1, 0, // development (fixed-length)
                      0,
                      p1_Tw, // death (fixed daily rate)
                      0, 0,
                      0,
                      0);
  if (test) return 1;
  //
  // Larvae
  test = spop_iterate((*conn2),
                      0,
                      d2, 0, // development (fixed-length)
                      0,
                      p2_Tw, // death (fixed daily rate)
                      0, 0,
                      0,
                      0);
  if (test) return 1;
  //
  // Pupae
  test = spop_iterate((*conn3),
                      0,
                      d3, 0, // development (fixed-length)
                      0,
                      p3_Tw, // death (fixed daily rate)
                      0, 0,
                      0,
                      0);
  if (test) return 1;
  //
  // Adult naive females
  test = spop_iterate((*conn4j),
                      0,
                      alpha_blood, 0, // development (fixed-length)
                      0,
                      0,
                      dd4, dd4s, // death (gamma-distributed)
                      0,
                      0);
  if (test) return 1;
  //
  // Adult mature females
  test = spop_iterate((*conn4),
                      0,
                      alpha_ovipos, 0, // gonotrophic cycle (fixed length)
                      0,
                      0,
                      dd4, dd4s, // death (gamma-distributed)
                      0,
                      0);
  if (test) return 1;
  //
  // Adult females surviving and completing ovicyle
  int n4_reproduce = (*conn4)->developed.i;
  // Reintroduce ovipositioning adult females to the population
  spop_popadd((*conn4), (*conn4)->devtable);
  //
  // Perform state transformations
  spop_popadd((*conn4), (*conn4j)->devtable);
  spop_add((*conn4j), 0, 0, 0, gsl_ran_binomial(RAND_GSL,0.5,(*conn3)->developed.i));
  spop_add((*conn3), 0, 0, 0, (*conn2)->developed.i);
  spop_add((*conn2), 0, 0, 0, (*conn1)->developed.i);
  //
  // Lay eggs
  int neweggs = gsl_ran_poisson(RAND_GSL,bigF4 * n4_reproduce); // Total number of eggs that will be laid that day
  int captured = 0; // Eggs captured with ovitraps
  int freeeggs = neweggs; // Eggs included in the life cycle
  if (daily_capture[TIME] > 0.5) {
    captured = gsl_ran_binomial(RAND_GSL,prob_capture,neweggs);
    freeeggs = neweggs - captured;
  }
  //
  spop_add((*conn1),0,0,0,freeeggs); // Normal eggs
  //
  // Update all stages and state variables
  (*n1) = (*conn1)->size.i;
  (*n2) = (*conn2)->size.i;
  (*n3) = (*conn3)->size.i;
  (*n4fj) = (*conn4j)->size.i;
  (*n4f) = (*conn4)->size.i;
  // Update indicators
  (*d4) = dd4;
  (*d4s) = dd4s;
  (*F4) = bigF4;
  (*cap) = captured;
  (*egg) = neweggs;
  //
  return 0;
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
    "coln1","coln2","coln3","coln4fj","coln4f","colnBS","colK","cold4","cold4s","colF4","colegg","colcap",
    "p1.1","p1.2","p1.3",
    "p2.1","p2.2","p2.3",
    "p3.1","p3.2","p3.3",
    "d4.1","d4.2","d4.3",
    "F4.1","F4.2","F4.3",
    "d1.1","d1.2","d1.3",
    "d2.1","d2.2","d2.3",
    "d3.1","d3.2","d3.3",
    "n23.surv",
    "deltaT",
    "BS.pdens","BS.dprec","BS.nevap","BS.Klulc",
    "tbm.1","tbm.2","tbm.3",
    "gtc.1","gtc.2","gtc.3",
    "n23.1","n23.2","n23.3","n23.4","n23.5",
    "alpha_capture",
    "alpha_init_egg",
    "alpha_flush_thr","alpha_flush_surv"
  };
  int i;
  for (i=0; i<(NumMetAea+NumParAea); i++)
    names[i] = strdup(temp[i]);
  // ---
  param[alpha_p1_1] = 0.9120081298856487;
  param[alpha_p1_2] = -5.362550627291735;
  param[alpha_p1_3] = 40.392649926349065;

  param[alpha_p2_1] = 0.9630337102917345;
  param[alpha_p2_2] = 3.1535505906512573;
  param[alpha_p2_3] = 37.39321020554058;

  param[alpha_p3_1] = 0.9989963317439844;
  param[alpha_p3_2] = 12.682001028351872;
  param[alpha_p3_3] = 35.88839257040825;

  param[alpha_d4_1] = 40.546569688329726;
  param[alpha_d4_2] = -3.0;
  param[alpha_d4_3] = 34.5;

  param[alpha_F4_1] = -29.983659263234784;
  param[alpha_F4_2] = 2.547065674987852;
  param[alpha_F4_3] = -0.03703376793024484;

  param[alpha_d1_1] = 13.468941558731592;
  param[alpha_d1_2] = -0.8652474262795845;
  param[alpha_d1_3] = 0.01891979289744318;

  param[alpha_d2_1] = 44.59021955694031;
  param[alpha_d2_2] = -2.613778489396158;
  param[alpha_d2_3] = 0.04627568522818735;

  param[alpha_d3_1] = 8.92627601829611;
  param[alpha_d3_2] = -0.36161169585829006;
  param[alpha_d3_3] = 0.0039784277144151396;

  param[alpha_n23_surv] = 3.873133206986801e-05;

  param[alpha_deltaT] = 0.0;

  param[alpha_BS_pdens] = 0.00001;
  param[alpha_BS_dprec] = 0.001;
  param[alpha_BS_nevap] = 0.9;
  param[alpha_BS_Klulc] = 0.5;

  param[alpha_tbm_1] = 29.99806440902287;
  param[alpha_tbm_2] = -1.8712756925767178;
  param[alpha_tbm_3] = 0.03514278986605401;

  param[alpha_gtc_1] = 29.99806440902287;
  param[alpha_gtc_2] = -1.8712756925767178;
  param[alpha_gtc_3] = 0.03514278986605401;

  param[alpha_n23_1] = 0.6080929959830579;
  param[alpha_n23_2] = 0.1059944669627986;
  param[alpha_n23_3] = -1.0925997735665425;
  param[alpha_n23_4] = 0.25464919919205214;
  param[alpha_n23_5] = -0.00504358886713356;

  param[alpha_capture] = -1.0;

  param[alpha_init_egg] = 100.0;

  param[alpha_flush_thr] = 10.0;
  param[alpha_flush_surv] = 0.25;
}

void sim_model(double               *envar,
               double               *param,
               int                 *finalT,
               int                *control,
               double              *result,
               int                *success) {
  double *mean_air_temp        = envar + 0*(*finalT);
  double *daily_precipitation  = envar + 1*(*finalT);
  double *popdens              = envar + 2*(*finalT);
  double *lulc                 = envar + 2*(*finalT) + 1;
  double *daily_capture        = envar + 2*(*finalT) + 2;
  double *controlpar = 0;
  if ((*control)) {
    controlpar = param + NumParAea;
  }
  //
  set_param(param);
  //
  double *colT      = result + 0*(*finalT);
  double *coln1     = result + 1*(*finalT);
  double *coln2     = result + 2*(*finalT);
  double *coln3     = result + 3*(*finalT);
  double *coln4fj   = result + 4*(*finalT);
  double *coln4f    = result + 5*(*finalT);
  double *colnBS    = result + 6*(*finalT);
  double *colK      = result + 7*(*finalT);
  double *cold4     = result + 8*(*finalT);
  double *cold4s    = result + 9*(*finalT);
  double *colF4     = result + 10*(*finalT);
  double *colegg    = result + 11*(*finalT);
  double *colcap    = result + 12*(*finalT);

  int TIME = 0;
  (*success) = 2;
  // Set initial conditions
  int n1 = max(0,round(param[alpha_init_egg]));
  int n2 = 0;
  int n3 = 0;
  int n4fj = 0;
  int n4f = 0;
  double nBS = 0.0;
  double K = 0.0;
  int d4 = 0;
  int d4s = 0;
  int F4 = 0;
  int egg = 0;
  int cap = 0;
  //
  spop conn1 = spop_init(1,gamma_mode);
  spop conn2 = spop_init(1,gamma_mode);
  spop conn3 = spop_init(1,gamma_mode);
  spop conn4 = spop_init(1,gamma_mode);
  spop conn4j = spop_init(1,gamma_mode);
  if (n1) spop_add(conn1,0,0,0,n1);
  if (n2) spop_add(conn2,0,0,0,n2);
  if (n3) spop_add(conn3,0,0,0,n3);
  if (n4fj) spop_add(conn4j,0,0,0,n4fj);
  if (n4f) spop_add(conn4,0,0,0,n4f);
  // Record state
  colT[TIME] = (double)TIME;
  coln1[TIME] = (double)n1;
  coln2[TIME] = (double)n2;
  coln3[TIME] = (double)n3;
  coln4fj[TIME] = (double)n4fj;
  coln4f[TIME] = (double)n4f;
  colnBS[TIME] = (double)nBS;
  colK[TIME] = (double)K;
  cold4[TIME] = (double)d4;
  cold4s[TIME] = (double)d4s;
  colF4[TIME] = (double)F4;
  colegg[TIME] = (double)egg;
  colcap[TIME] = (double)cap;
  //
  char test = 0;
  for (TIME=1; TIME<(*finalT); TIME++) {
    // Take a step
    test = calculate(mean_air_temp,
                     daily_precipitation,
                     popdens,
                     lulc,
                     daily_capture,
                     param,
                     &conn1,
                     &conn2,
                     &conn3,
                     &conn4j,
                     &conn4,
                     &n1,
                     &n2,
                     &n3,
                     &n4fj,
                     &n4f,
                     &nBS,
                     &K,
                     &d4,
                     &d4s,
                     &F4,
                     &egg,
                     &cap,
                     TIME);
    if (test) {
      (*success) = 0;
      goto endall;
    }
    //
    if ((*control)) { // Apply control measures
      double tmp;
      //
      // Breeding site reduction (date0 - date1 - daily fraction)
      if (TIME>=controlpar[0] && TIME<controlpar[1]) {
        nBS = gsl_ran_binomial(RAND_GSL,1.0 - controlpar[2],nBS);
      }
      // Egg reduction
      if (TIME>=controlpar[3] && TIME<controlpar[4]) {
        tmp = 1.0 - controlpar[5];
        n1 = gsl_ran_binomial(RAND_GSL,tmp,n1);
        spop_iterate(conn1,
                     0,
                     0, 0,
                     0,
                     tmp,
                     0, 0,
                     0,
                     1);
      }
      // Larva reduction
      if (TIME>=controlpar[6] && TIME<controlpar[7]) {
        tmp = 1.0 - controlpar[8];
        n2 = gsl_ran_binomial(RAND_GSL,tmp,n2);
        spop_iterate(conn2,
                     0,
                     0, 0,
                     0,
                     tmp,
                     0, 0,
                     0,
                     1);
      }
      // Pupa reduction
      if (TIME>=controlpar[9] && TIME<controlpar[10]) {
        tmp = 1.0 - controlpar[11];
        n3 = gsl_ran_binomial(RAND_GSL,tmp,n3);
        spop_iterate(conn3,
                     0,
                     0, 0,
                     0,
                     tmp,
                     0, 0,
                     0,
                     1);
      }
      // Adult reduction
      if (TIME>=controlpar[12] && TIME<controlpar[13]) {
        tmp = 1.0 - controlpar[14];
        n4fj = gsl_ran_binomial(RAND_GSL,tmp,n4fj);
        n4f = gsl_ran_binomial(RAND_GSL,tmp,n4f);
        spop_iterate(conn4j,
                     0,
                     0, 0,
                     0,
                     tmp,
                     0, 0,
                     0,
                     1);
        spop_iterate(conn4,
                     0,
                     0, 0,
                     0,
                     tmp,
                     0, 0,
                     0,
                     1);
      }
    }
    // Record state
    colT[TIME] = (double)TIME;
    coln1[TIME] = (double)n1;
    coln2[TIME] = (double)n2;
    coln3[TIME] = (double)n3;
    coln4fj[TIME] = (double)n4fj;
    coln4f[TIME] = (double)n4f;
    colnBS[TIME] = (double)nBS;
    colK[TIME] = (double)K;
    cold4[TIME] = (double)d4;
    cold4s[TIME] = (double)d4s;
    colF4[TIME] = (double)F4;
    colegg[TIME] = (double)egg;
    colcap[TIME] = (double)cap;
    if (CHECK(n1) ||
        CHECK(n2) ||
        CHECK(n3) ||
        CHECK(n4fj) ||
        CHECK(n4f) ||
        isnan(nBS) ||
        isnan(K) ||
        CHECK(d4) ||
        CHECK(d4s) ||
        CHECK(F4) ||
        CHECK(egg) ||
        CHECK(cap)) {
      (*success) = 0;
      goto endall;
    }
  }
  (*success) = 1;
  //
 endall:
  //
  (*success) = ((*success)==0) ? 0 : 1;
  spop_destroy(&conn1);
  spop_destroy(&conn2);
  spop_destroy(&conn3);
  spop_destroy(&conn4j);
  spop_destroy(&conn4);
  //
  gamma_dist_check();
  //
  //printf("Times: %g %g %g\n",timebefore,timeof,timeafter);
}

// ---------------------------------------------------------------------------


