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

#define NumParAea      43
#define NumMetAea      18

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
#define alpha_F4_4        15
#define alpha_d1_1        16
#define alpha_d1_2        17
#define alpha_d1_3        18
#define alpha_d2_1        19
#define alpha_d2_2        20
#define alpha_d2_3        21
#define alpha_d3_1        22
#define alpha_d3_2        23
#define alpha_d3_3        24
#define alpha_n23_surv    25
#define alpha_deltaT      26

#define alpha_BS_pdens    27
#define alpha_BS_dprec    28
#define alpha_BS_nevap    29
#define alpha_K           30

#define alpha_dp_egg      31
#define alpha_dp_thr      32
#define alpha_ta_thr      33
#define alpha_tq_thr      34

#define alpha_tbm_1       35
#define alpha_tbm_2       36
#define alpha_tbm_3       37

#define alpha_p0_1        38
#define alpha_p0_2        39

#define alpha_capture     40

#define alpha_flush_thr   41
#define alpha_flush_surv  42

#define alpha_ta_thr_std  (0.5)
#define alpha_dp_thr_std  (10.0/24.0/60.0)
#define alpha_tq_thr_std  (0.5)

#define flin(x, a1, a2) (max(0.0, min(1.0, (a1) + (a2)*(x))))
#define dsig(x, a1, a2, a3) (max(0.0, min(1.0, (a1)/((1.0+exp((a2)-(x)))*(1.0+exp((x)-(a3)))))))
#define dsig2(x, a1, a2, a3) (max(0.0, min(100.0, (a1)/((1.0+exp((a2)-(x)))*(1.0+exp((x)-(a3)))))))
#define poly(x, a1, a2, a3) (max(0.0, (a1) + (a2)*(x) + (a3)*pow((x),2)))
#define Tthr(T,Tmn,Tstd) (1.0/(1.0+exp(((Tmn)-(T))/(Tstd))))
#define pdiap(P,Pmn,Pstd,tt) (1.0-(tt)-(1.0-(tt))/((1.0+exp(((Pmn)-(P))/(Pstd)))))
#define expo(x, n) (exp(-pow((x),(n))))
#define ratk(x, T1, T2, K, c) (pow((c)*((x)-(T1))*(1.0-exp((K)*((x)-(T2)))),2))

double fun_p0Ta(double Ta, double *param) {return flin(Ta,param[alpha_p0_1],param[alpha_p0_2]);}
double fun_p1Tw(double Tw, double *param) {return dsig(Tw,param[alpha_p1_1],param[alpha_p1_2],param[alpha_p1_3]);}
double fun_p2Tw_raw(double Tw, double *param) {return dsig(Tw,param[alpha_p2_1],param[alpha_p2_2],param[alpha_p2_3]);}
double fun_p3Tw_raw(double Tw, double *param) {return dsig(Tw,param[alpha_p3_1],param[alpha_p3_2],param[alpha_p3_3]);}
double fun_d1Tw(double Tw, double *param) {return max(1.0, poly(Tw,param[alpha_d1_1],param[alpha_d1_2],param[alpha_d1_3]));}
double fun_d2Tw(double Tw, double *param) {return max(1.0, poly(Tw,param[alpha_d2_1],param[alpha_d2_2],param[alpha_d2_3]));}
double fun_d3Tw(double Tw, double *param) {return max(1.0, poly(Tw,param[alpha_d3_1],param[alpha_d3_2],param[alpha_d3_3]));}
double fun_F4(double Ta, double *param) {return Ta<=param[alpha_F4_1] || Ta>=param[alpha_F4_2] ? 0.0 : ratk(Ta,param[alpha_F4_1],param[alpha_F4_2],param[alpha_F4_3],param[alpha_F4_4]);}
double fun_dens(double n23dens, double *param) {return expo(n23dens,param[alpha_n23_surv]);}
double fun_tbmTa(double Ta, double *param) {return poly(Ta,param[alpha_tbm_1],param[alpha_tbm_2],param[alpha_tbm_3]);}
double fun_d4Ta_mean(double Ta, double *param) {return dsig2(Ta,param[alpha_d4_1],param[alpha_d4_2],param[alpha_d4_3]);}
double fun_d4Ta_std(double dd4) {return gamma_matrix_sd*dd4;}

//double timebefore = 0;
//double timeof = 0;
//double timeafter = 0;

double prob_capture = 0;

void set_param(double *param) {
  prob_capture = pow(10,param[alpha_capture]);
}

char calculate(double *photoperiod,
               double *mean_air_temp,
               double *daily_precipitation,
               double *popdens,
               double *daily_capture,
               double *param,
               spop   *conn10,
               spop   *conn1,
               spop   *conn2,
               spop   *conn3,
               spop   *conn4j,
               spop   *conn4,
               int    *n0,
               int    *n10,
               int    *n1,
               int    *nh,
               int    *n2,
               int    *n3,
               int    *n4fj,
               int    *n4f,
               double *nBS,
               int    *d4,
               int    *d4s,
               int    *F4,
               int    *egg,
               int    *cap,
               double *fhatch,
               double *fdiap,
               double *fquie,
               int    TIME) {
  // ---------------------
  // modelDelayAalbopictus
  // ---------------------

  double deltaT = param[alpha_deltaT];

  double Ta = mean_air_temp[TIME];
  
  double Tw = Ta + deltaT;

  /*
   * TO DO: Check how this can be implemented with a stochastic dynamics
   *
   * Update the number of breeding sites
   */
  (*nBS) = param[alpha_BS_pdens] * (popdens[TIME]) +
    param[alpha_BS_dprec] * daily_precipitation[TIME] +
    param[alpha_BS_nevap] * (*nBS);
  // (*K) = param[alpha_BS_nevap] == 1.0 ? (*nBS) / (TIME+1.0) : (*nBS) * (param[alpha_BS_nevap]-1.0) / (pow(param[alpha_BS_nevap], (TIME+1))-1.0);

  /*
   * Update survival and development rates/probabilities according to
   * the new environmental conditions and breeding sites
   */
  // Density of the immature stages
  // Assume uniform distribution across all breeding sites
  // 2019/09/23 - ErgulerK: K -> nBS
  // 2019/09/24 - ErgulerK: Instead of the average of previous precipitation events, I propose to use the aggregate, i.e. the weighted sum
  double n23dens = (double)((*n2) + (*n3)) / (param[alpha_K] * (double)(*nBS));
  //
  // The effect of density on the survival of immature stages
  double densd = fun_dens(n23dens,param);
  //
  // Fecundity (per adult female per day)
  // Since oviposition frequency is different for Aedes albopictus,
  // it is best to model fecundity as a daily event.
  double bigF4 = fun_F4(Ta,param);
  //
  // Excess rain flushes the aquatic life stages
  double flush_surv = daily_precipitation[TIME] > param[alpha_flush_thr] ? param[alpha_flush_surv] : 1.0;
  //
  // Egg survival (diapausing eggs)
  double p0_Ta = fun_p0Ta(Ta,param);
  // Egg mortality (non-diapausing eggs)
  double p1_Tw = 1.0 - fun_p1Tw(Tw,param);
  // Larval mortality
  double p2_Tw = 1.0 - fun_p2Tw_raw(Tw,param)*densd*flush_surv;
  // Pupal mortality
  double p3_Tw = 1.0 - fun_p3Tw_raw(Tw,param)*densd*flush_surv;
  
  // Egg development time
  double d1 = fun_d1Tw(Tw,param);
  // Larval development time
  double d2 = fun_d2Tw(Tw,param);
  // Pupal development time
  double d3 = fun_d3Tw(Tw,param);

  // Time to first blood meal
  double alpha_ovipos = fun_tbmTa(Ta,param);

  // Adult lifetime (from emergence)
  double dd4 = fun_d4Ta_mean(Ta,param);
  double dd4s = fun_d4Ta_std(dd4);

  // Diapause and quiescence
  // quie -> larvae
  double fracHatch = Tthr(Ta,param[alpha_ta_thr],alpha_ta_thr_std);
  // neweggs -> n10
  double fracDiap = pdiap(photoperiod[TIME],param[alpha_dp_thr],alpha_dp_thr_std,fracHatch);
  // n0 -> quie
  double fracQuie = 1.0 - Tthr(Ta,param[alpha_tq_thr],alpha_tq_thr_std);
  fracQuie *= 1.0 - fracDiap; // Do not allow quiescence during diapause conditions
  
  /*
   * Survival first, and then, development.
   * Update development times and survival proportions
   * of the immature stages and juvenile adults
   */
  //
  char test = 0;
  // Normal and tagged eggs
  test = spop_iterate((*conn1),
                      0,
                      d1, 0, // development (fixed-length)
                      0,
                      p1_Tw, // death (fixed daily rate)
                      0, 0,
                      0,
                      0);
  if (test) return 1;
  test = spop_iterate((*conn10),
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
                      alpha_ovipos, 0, // development (fixed-length)
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
                      0, 0, // daily fecundity
                      0,
                      0,
                      dd4, dd4s, // death (gamma-distributed)
                      0,
                      0);
  if (test) return 1;
  //
  int quie = gsl_ran_binomial(RAND_GSL,p0_Ta,(*nh)); // Number of surviving quiescent eggs (survival is similar to diapausing eggs)
  int quieh = gsl_ran_binomial(RAND_GSL,fracHatch,quie); // Number of surviving quiescent eggs to hatch
  quie -= quieh;
  //
  // Diapausing egg survival
  // Tagged eggs always become diapausing eggs
  int dp_eggs = gsl_ran_binomial(RAND_GSL,p0_Ta,(*n0));
  dp_eggs += (*conn10)->developed.i;
  //
  // WARNING: Here, the eggs are produced and allowed to become quiescent immediately
  //          However, they need to be exposed to cold temperatures before being allowed to hatch
  int nquie = gsl_ran_binomial(RAND_GSL,fracQuie,dp_eggs);
  quie += nquie;
  dp_eggs -= nquie;
  //
  // Adult females surviving and ovipositioning
  int n4_reproduce = (*conn4)->size.i;
  //
  // Perform state transformations
  spop_popadd((*conn4), (*conn4j)->devtable);
  spop_add((*conn4j), 0, 0, 0, gsl_ran_binomial(RAND_GSL,0.5,(*conn3)->developed.i));
  spop_add((*conn3), 0, 0, 0, (*conn2)->developed.i);
  spop_add((*conn2), 0, 0, 0, (*conn1)->developed.i + quieh);
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
  int tagged = gsl_ran_binomial(RAND_GSL,fracDiap,freeeggs);
  spop_add((*conn10),0,0,0,tagged); // Tagged eggs
  spop_add((*conn1),0,0,0,(freeeggs-tagged)); // Normal eggs
  //
  // Update all stages and state variables
  (*n0) = dp_eggs;
  (*n10) = (*conn10)->size.i;
  (*n1) = (*conn1)->size.i;
  (*nh) = quie;
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
  (*fhatch) = fracHatch;
  (*fdiap) = fracDiap;
  (*fquie) = fracQuie;
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
    "coln0","coln10","coln1","colnh","coln2","coln3","coln4fj","coln4f","colnBS","cold4","cold4s","colF4","colegg","colcap","colhatch","coldiap","colquie","alpha_flush_thr","alpha_flush_surv",
    "p1.1","p1.2","p1.3","p2.1","p2.2","p2.3","p3.1","p3.2","p3.3","d4.1","d4.2","d4.3","F4.1","F4.2","F4.3","F4.4","d1.1","d1.2","d1.3","d2.1","d2.2","d2.3","d3.1","d3.2","d3.3","n23.surv","deltaT","BS.pdens","BS.dprec","BS.nevap","PP.init","PP.thr","PP.ta.thr","PP.tq.thr","tbm.1","tbm.2","tbm.3","p0.1","p0.2","alpha_capture"
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
  param[alpha_F4_1] = 15.4466728;
  param[alpha_F4_2] = 44.8260538;
  param[alpha_F4_3] = -0.0137863654;
  param[alpha_F4_4] = 1.09767332;
  param[alpha_d1_1] = 13.468941558731592;
  param[alpha_d1_2] = -0.8652474262795845;
  param[alpha_d1_3] = 0.01891979289744318;
  param[alpha_d2_1] = 44.59021955694031;
  param[alpha_d2_2] = -2.613778489396158;
  param[alpha_d2_3] = 0.04627568522818735;
  param[alpha_d3_1] = 8.92627602;
  param[alpha_d3_2] = -0.377139162;
  param[alpha_d3_3] = 0.00464125318;
  param[alpha_n23_surv] = 10.0;
  param[alpha_deltaT] = 0.0;

  param[alpha_BS_pdens] = 0.00001;
  param[alpha_BS_dprec] = 0.001;
  param[alpha_BS_nevap] = 0.9;
  param[alpha_K]        = 1.0;

  param[alpha_dp_egg] = 1.0;
  param[alpha_dp_thr] = 0.4;
  param[alpha_ta_thr] = 20.0;
  param[alpha_tq_thr] = 10.0;

  param[alpha_tbm_1] = 52.54632144;
  param[alpha_tbm_2] = -3.81802776;
  param[alpha_tbm_3] = 0.07403983;

  param[alpha_p0_1] = 0.7975103960182536;
  param[alpha_p0_2] = 0.07409437745556843;

  param[alpha_capture] = -1.0;

  param[alpha_flush_thr] = 10.0;
  param[alpha_flush_surv] = 0.25;
}

void sim_model(double               *envar,
               double               *param,
               int                 *finalT,
               int                *control,
               double              *result,
               int                *success) {
  double *photoperiod          = envar + 0*(*finalT);
  double *mean_air_temp        = envar + 1*(*finalT);
  double *daily_precipitation  = envar + 2*(*finalT);
  double *popdens              = envar + 3*(*finalT);
  double *daily_capture        = envar + 4*(*finalT);
  double *controlpar = 0;
  if ((*control)) {
    controlpar = param + NumParAea;
  }
  //
  set_param(param);
  //
  double *colT      = result + 0*(*finalT);
  double *coln0     = result + 1*(*finalT);
  double *coln10    = result + 2*(*finalT);
  double *coln1     = result + 3*(*finalT);
  double *colnh     = result + 4*(*finalT);
  double *coln2     = result + 5*(*finalT);
  double *coln3     = result + 6*(*finalT);
  double *coln4fj   = result + 7*(*finalT);
  double *coln4f    = result + 8*(*finalT);
  double *colnBS    = result + 9*(*finalT);
  double *cold4     = result + 10*(*finalT);
  double *cold4s    = result + 11*(*finalT);
  double *colF4     = result + 12*(*finalT);
  double *colegg    = result + 13*(*finalT);
  double *colcap    = result + 14*(*finalT);
  double *colhatch  = result + 15*(*finalT);
  double *coldiap   = result + 16*(*finalT);
  double *colquie   = result + 17*(*finalT);

  int TIME = 0;
  (*success) = 2;
  // Set initial conditions
  int n0 = max(0,round(param[alpha_dp_egg]));
  int n10 = 0;
  int n1 = 0;
  int nh = 0;
  int n2 = 0;
  int n3 = 0;
  int n4fj = 0;
  int n4f = 0;
  double nBS = 0.0;
  int d4 = 0;
  int d4s = 0;
  int F4 = 0;
  int egg = 0;
  int cap = 0;
  double fhatch = 0.0;
  double fdiap = 0.0;
  double fquie = 0.0;
  //
  spop conn10 = spop_init(1,gamma_mode);
  spop conn1 = spop_init(1,gamma_mode);
  spop conn2 = spop_init(1,gamma_mode);
  spop conn3 = spop_init(1,gamma_mode);
  spop conn4 = spop_init(1,gamma_mode);
  spop conn4j = spop_init(1,gamma_mode);
  if (n10) spop_add(conn10,0,0,0,n10);
  if (n1) spop_add(conn1,0,0,0,n1);
  if (n2) spop_add(conn2,0,0,0,n2);
  if (n3) spop_add(conn3,0,0,0,n3);
  if (n4fj) spop_add(conn4j,0,0,0,n4fj);
  if (n4f) spop_add(conn4,0,0,0,n4f);
  // Record state
  colT[TIME] = (double)TIME;
  coln0[TIME] = (double)n0;
  coln10[TIME] = (double)n10;
  coln1[TIME] = (double)n1;
  colnh[TIME] = (double)nh;
  coln2[TIME] = (double)n2;
  coln3[TIME] = (double)n3;
  coln4fj[TIME] = (double)n4fj;
  coln4f[TIME] = (double)n4f;
  colnBS[TIME] = (double)nBS;
  cold4[TIME] = (double)d4;
  cold4s[TIME] = (double)d4s;
  colF4[TIME] = (double)F4;
  colegg[TIME] = (double)egg;
  colcap[TIME] = (double)cap;
  colhatch[TIME] = (double)fhatch;
  coldiap[TIME] = (double)fdiap;
  colquie[TIME] = (double)fquie;
  //
  char test = 0;
  for (TIME=1; TIME<(*finalT); TIME++) {
    // Take a step
    test = calculate(photoperiod,
                     mean_air_temp,
                     daily_precipitation,
                     popdens,
                     daily_capture,
                     param,
                     &conn10,
                     &conn1,
                     &conn2,
                     &conn3,
                     &conn4j,
                     &conn4,
                     &n0,
                     &n10,
                     &n1,
                     &nh,
                     &n2,
                     &n3,
                     &n4fj,
                     &n4f,
                     &nBS,
                     &d4,
                     &d4s,
                     &F4,
                     &egg,
                     &cap,
                     &fhatch,
                     &fdiap,
                     &fquie,
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
        n0 = gsl_ran_binomial(RAND_GSL,tmp,n0);
        n10 = gsl_ran_binomial(RAND_GSL,tmp,n10);
        n1 = gsl_ran_binomial(RAND_GSL,tmp,n1);
        nh = gsl_ran_binomial(RAND_GSL,tmp,nh);
        spop_iterate(conn1,
                     0,
                     0, 0,
                     0,
                     tmp,
                     0, 0,
                     0,
                     1);
        spop_iterate(conn10,
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
    coln0[TIME] = (double)n0;
    coln10[TIME] = (double)n10;
    coln1[TIME] = (double)n1;
    colnh[TIME] = (double)nh;
    coln2[TIME] = (double)n2;
    coln3[TIME] = (double)n3;
    coln4fj[TIME] = (double)n4fj;
    coln4f[TIME] = (double)n4f;
    colnBS[TIME] = (double)nBS;
    cold4[TIME] = (double)d4;
    cold4s[TIME] = (double)d4s;
    colF4[TIME] = (double)F4;
    colegg[TIME] = (double)egg;
    colcap[TIME] = (double)cap;
    colhatch[TIME] = (double)fhatch;
    coldiap[TIME] = (double)fdiap;
    colquie[TIME] = (double)fquie;
    if (CHECK(n0) || CHECK(n10) || CHECK(n1) || CHECK(nh) || CHECK(n2) || CHECK(n3) || CHECK(n4fj) || CHECK(n4f) || isnan(nBS) || CHECK(d4) || CHECK(d4s) || CHECK(F4) || CHECK(egg) || CHECK(cap) || isnan(fhatch) || isnan(fdiap) || isnan(fquie)) {
      (*success) = 0;
      goto endall;
    }
  }
  (*success) = 1;
  //
 endall:
  //
  (*success) = ((*success)==0) ? 0 : 1;
  spop_destroy(&conn10);
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


