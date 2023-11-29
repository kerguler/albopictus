#include "Python.h"
#include <math.h>
#include <time.h>
#include "gamma.h"
#include "incubator.h"

volatile clock_t start = 0, diff;
double time2here(void) {
  diff = clock() - start;
  start = clock();
  return diff;
}

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

#define NumParAea      48
#define NumMetAea      12

// --------------------------------------------
// Gamma distribution
// --------------------------------------------

char gamma_mode = MODE_GAMMA_HASH;
void set_gamma_mode(char mode) {
  gamma_mode = mode;
}

// --------------------------------------------
// Incubators
// --------------------------------------------

incubator conn1;
incubator conn10;
incubator conn2;
incubator conn3;
incubator conn4;
//
void update(double *size, double *dev, double *par) {
  *size *= par[0];
  *dev += par[1];
}
//
void empty_incubators(void) {
  //
  incubator_empty(&conn10);
  incubator_empty(&conn1);
  incubator_empty(&conn2);
  incubator_empty(&conn3);
  incubator_empty(&conn4);
  //
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

#define alpha_0           25
#define alpha_1           26
#define alpha_2           27
#define alpha_3           28
#define alpha_deltaT      29

#define alpha_BS_pdens    30
#define alpha_BS_dprec    31
#define alpha_BS_nevap    32

#define alpha_dp_egg      33
#define alpha_dp_thr      34
#define alpha_ta_thr      35
#define alpha_rate_strong 36
#define alpha_rate_normal 37

#define alpha_tbm_1       38
#define alpha_tbm_2       39
#define alpha_tbm_3       40

#define alpha_p0_1        41
#define alpha_p0_2        42

#define alpha_n23_1       43
#define alpha_n23_2       44
#define alpha_n23_3       45
#define alpha_n23_4       46
#define alpha_n23_5       47

#define flin(x, a1, a2) (max(0.0, min(1.0, (a1) + (a2)*(x))))
#define dsig(x, a1, a2, a3) (max(0.0, min(1.0, (a1)/((1.0+exp((a2)-(x)))*(1.0+exp((x)-(a3)))))))
#define dsig2(x, a1, a2, a3) (max(0.0, min(100.0, (a1)/((1.0+exp((a2)-(x)))*(1.0+exp((x)-(a3)))))))
#define poly(x, a1, a2, a3) (max(0.0, (a1) + (a2)*(x) + (a3)*pow((x),2)))
#define expd(x, k) (exp(-(k)*(x)))
#define fpow(x, a1, a2) (min(1.0, (a1)*pow((x),(a2))))

//double timebefore = 0;
//double timeof = 0;
//double timeafter = 0;

void calculate(double *photoperiod,
               double *mean_air_temp,
               double *daily_precipitation,
               double *popdens,
               double *param,
               double *n0,
               double *n10,
               double *n1,
               double *n2,
               double *n3,
               double *n4fj,
               double *n4f,
               double *nBS,
               double *K,
               double *d4,
               double *d4s,
               double *F4,
               double *egg,
               double *percent_strong,
               double *diap,
               int    TIME) {
  double par[2];

  if (TIME==-1) {
    (*n0) = param[alpha_dp_egg];
    (*n10) = 0.0;
    (*n1) = param[alpha_0];
    (*n2) = param[alpha_1];
    (*n3) = param[alpha_2];
    (*n4fj) = 0.0;
    (*n4f) = param[alpha_3];
    (*nBS) = 0.0;
    (*K) = 0.0;
    (*d4) = 0.0;
    (*d4s) = 0.0;
    (*F4) = 0.0;
    (*egg) = (*n0) + (*n10) + (*n1);
    (*diap) = 0;
    (*percent_strong) = 0.0;

    if ((*n10) > 0) incubator_add(&conn10,(*n10),0.0);
    if ((*n1) > 0) incubator_add(&conn1,(*n1),0.0);
    if ((*n2) > 0) incubator_add(&conn2,(*n2),0.0);
    if ((*n3) > 0) incubator_add(&conn3,(*n3),0.0);
    if ((*n4f) > 0) incubator_add(&conn4,(*n4f),0.0);

    return;
  }

  // ---------------------
  // modelDelayAalbopictus
  // ---------------------

  double deltaT = param[alpha_deltaT];

  double Ta = mean_air_temp[TIME];
  
  double Tw = Ta+deltaT;

  // Number of breeding sites
  (*nBS) = param[alpha_BS_pdens]*(popdens[TIME]) + param[alpha_BS_dprec]*daily_precipitation[TIME] + param[alpha_BS_nevap]*(*nBS);
  (*K) = param[alpha_BS_nevap]==1.0 ? (*nBS)/(TIME+1.0) : (*nBS)*(param[alpha_BS_nevap]-1.0)/(pow(param[alpha_BS_nevap],(TIME+1))-1.0);

  // Density of the immature stages
  // Assume uniform distribution across all breeding sites
  double n23dens = ((incubator_sum(conn2)+incubator_sum(conn3))/(*K));

  // Fecundity
  double bigF4 = poly(Ta,param[alpha_F4_1],param[alpha_F4_2],param[alpha_F4_3]);
  (*F4) = bigF4;

  double densd = expd(n23dens,param[alpha_n23_surv]);
  // Egg survival (diapausing eggs)
  double p0_Ta = flin(Ta,param[alpha_p0_1],param[alpha_p0_2]);
  // Egg survival (non-diapausing eggs)
  double p1_Tw = dsig(Tw,param[alpha_p1_1],param[alpha_p1_2],param[alpha_p1_3]);
  // Larval survival
  double p2_Tw = dsig(Tw,param[alpha_p2_1],param[alpha_p2_2],param[alpha_p2_3])*densd;
  // Pupal survival
  double p3_Tw = dsig(Tw,param[alpha_p3_1],param[alpha_p3_2],param[alpha_p3_3])*densd;

  double densdev = fpow(n23dens,param[alpha_n23_1],param[alpha_n23_2]*poly(Tw,param[alpha_n23_3],param[alpha_n23_4],param[alpha_n23_5]));
  // Egg development time
  double d1 = poly(Tw,param[alpha_d1_1],param[alpha_d1_2],param[alpha_d1_3]);
  // Larval development time
  double d2 = poly(Tw,param[alpha_d2_1],param[alpha_d2_2],param[alpha_d2_3])*densdev;
  // Pupal development time
  double d3 = poly(Tw,param[alpha_d3_1],param[alpha_d3_2],param[alpha_d3_3])*densdev;
  // Time to first blood meal
  double alpha_blood = poly(Ta,param[alpha_tbm_1],param[alpha_tbm_2],param[alpha_tbm_3]);
  // Adult lifetime (from emergence)
  double dd4 = dsig2(Ta,param[alpha_d4_1],param[alpha_d4_2],param[alpha_d4_3]);
  double dd4s = gamma_matrix_sd*dd4;
  (*d4) = dd4;
  (*d4s) = dd4s;

  // Update development times and survival proportions
  // of the immature stages and juvenile adults
  //
  // Diapausing eggs
  (*n0) = p0_Ta*(*n0);
  //
  // Normal and tagged eggs
  par[0] = p1_Tw; par[1] = 1.0/max(1.0,d1);
  incubator_update(conn1,update,par);
  incubator_update(conn10,update,par);
  //
  // Larvae
  par[0] = p2_Tw; par[1] = 1.0/max(1.0,d2);
  incubator_update(conn2,update,par);
  //
  // Pupae
  par[0] = p3_Tw; par[1] = 1.0/max(1.0,d3);
  incubator_update(conn3,update,par);
  //
  // Adult females
  incubator_develop_survive(&conn4,-1,0,dd4,dd4s,alpha_blood,n4fj,n4f,0,gamma_mode);
  //
  // Check if it is winter or summer
  char short_days = photoperiod[TIME] < param[alpha_dp_thr];
  char cold_days = Ta < param[alpha_ta_thr];
  //
  // Lay eggs
  (*egg) = bigF4*(*n4f); // Total number of eggs that will be laid that day
  double hatch = 0;
  //
  if (short_days && cold_days) { // DIAPAUSE
    (*diap) = 1.0;
    // The fraction of tagged eggs increases linearly
    (*percent_strong) = min(1.0,(*percent_strong)+param[alpha_rate_strong]);
    //
    incubator_add(&conn10,bigF4*(*n4f)*(*percent_strong),0.0); // Tagged eggs
    incubator_add(&conn1,bigF4*(*n4f)*(1.0-(*percent_strong)),0.0); // Normal eggs
    //
  } else { // NO DIAPAUSE
    (*diap) = 0.0;
    if (!short_days && !cold_days) { // EXIT FROM DIAPAUSE
      (*diap) = 0.5;
      // Prepare diapausing eggs for hatching
      double vn0_n10 = (*n0)*param[alpha_rate_normal];
      (*n0) -= vn0_n10; // Diapausing eggs
      hatch += vn0_n10; // Eggs to hatch
    }
    // Lay normal eggs
    incubator_add(&conn1,bigF4*(*n4f),0.0); // Normal eggs
    //
    (*percent_strong) = 0.0;
  }
  // Develop
  double nn1;
  double nn2;
  double nn3;
  incubator_remove(&conn3,&nn3);
  incubator_remove(&conn2,&nn2); incubator_add(&conn3,nn2,0.0);
  incubator_remove(&conn1,&nn1); hatch += nn1;
  // Harvest developed tagged eggs
  double vn10_hn0;
  incubator_remove(&conn10,&vn10_hn0);
  // Tagged eggs always become diapausing eggs
  (*n0) = (*n0) + vn10_hn0;
  // All eggs, which are ready to hatch, become larvae
  incubator_add(&conn2,hatch,0.0);
  //
  // Update immature stage counts
  (*n3) = incubator_sum(conn3);
  (*n2) = incubator_sum(conn2);
  (*n1) = incubator_sum(conn1) + incubator_sum(conn10);
  //
  // Add newly developed females to adults
  nn3 *= 0.5;
  incubator_add(&conn4,nn3,0.0);
  (*n4f) += nn3;
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
    "coln0","coln1","coln2","coln3","coln4fj","coln4f","colK","cold4","cold4s","colF4","colegg","coldiap",
    "p1.1","p1.2","p1.3","p2.1","p2.2","p2.3","p3.1","p3.2","p3.3","d4.1","d4.2","d4.3","F4.1","F4.2","F4.3","d1.1","d1.2","d1.3","d2.1","d2.2","d2.3","d3.1","d3.2","d3.3","n23.surv","new.init1","new.init2","new.init3","new.init4","new.deltaT","BS.pdens","BS.dprec","BS.nevap","PP.init","PP.thr","PP.ta.thr","PP.strong","PP.normal","tbm.1","tbm.2","tbm.3","p0.1","p0.2","n23.1","n23.2","n23.3","n23.4","n23.5"
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

  param[alpha_0] = 0;
  param[alpha_1] = 0;
  param[alpha_2] = 0;
  param[alpha_3] = 0;
  param[alpha_deltaT] = 0.0;

  param[alpha_BS_pdens] = 0.00001;
  param[alpha_BS_dprec] = 0.001;
  param[alpha_BS_nevap] = 0.9;

  param[alpha_dp_egg] = 1.0;
  param[alpha_dp_thr] = 0.5;
  param[alpha_ta_thr] = 21.0;
  param[alpha_rate_strong] = 1.0;
  param[alpha_rate_normal] = 1.0;

  param[alpha_tbm_1] = 29.99806440902287;
  param[alpha_tbm_2] = -1.8712756925767178;
  param[alpha_tbm_3] = 0.03514278986605401;

  param[alpha_p0_1] = 0.7975103960182536;
  param[alpha_p0_2] = 0.07409437745556843;

  param[alpha_n23_1] = 0.6080929959830579;
  param[alpha_n23_2] = 0.1059944669627986;
  param[alpha_n23_3] = -1.0925997735665425;
  param[alpha_n23_4] = 0.25464919919205214;
  param[alpha_n23_5] = -0.00504358886713356;
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
  double *controlpar = 0;
  if ((*control)) {
    controlpar = param + NumParAea;
  }
	       
  double *colT    = result + 0*(*finalT);
  double *coln0   = result + 1*(*finalT);
  double *coln1   = result + 2*(*finalT);
  double *coln2   = result + 3*(*finalT);
  double *coln3   = result + 4*(*finalT);
  double *coln4fj = result + 5*(*finalT);
  double *coln4f  = result + 6*(*finalT);
  double *colK    = result + 7*(*finalT);
  double *cold4   = result + 8*(*finalT);
  double *cold4s  = result + 9*(*finalT);
  double *colF4   = result + 10*(*finalT);
  double *colegg  = result + 11*(*finalT);
  double *coldiap = result + 12*(*finalT);

  double n0 = 0;
  double n1 = 0;
  double n2 = 0;
  double n3 = 0;
  double n4fj = 0;
  double n4f = 0;
  double nBS = 0;
  double K = 0;
  double d4 = 0;
  double d4s = 0;
  double F4 = 0;
  double egg = 0;
  double diap = 0;

  double percent_strong = 0;
  double n10 = 0;

  int TIME;
  for (TIME=-1; TIME<(*finalT)-1; ) {
    //
    calculate(photoperiod,mean_air_temp,daily_precipitation,popdens,param,&n0,&n10,&n1,&n2,&n3,&n4fj,&n4f,&nBS,&K,&d4,&d4s,&F4,&egg,&percent_strong,&diap,TIME);
    TIME++;
    //
    if ((*control)) { // Apply control measures
      // Breeding site reduction (date0 - date1 - daily fraction)
      if (TIME>=controlpar[0] && TIME<controlpar[1]) {
        nBS *= controlpar[2];
      }
      // Egg reduction
      if (TIME>=controlpar[3] && TIME<controlpar[4]) {
        n0 *= controlpar[5];
        n10 *= controlpar[5];
        n1 *= controlpar[5];
        double par[2];
        par[0] = controlpar[5]; par[1] = 0;
        incubator_update(conn1,update,par);
        incubator_update(conn10,update,par);
      }
      // Larva reduction
      if (TIME>=controlpar[6] && TIME<controlpar[7]) {
        n2 *= controlpar[8];
        double par[2];
        par[0] = controlpar[8]; par[1] = 0;
        incubator_update(conn2,update,par);
      }
      // Pupa reduction
      if (TIME>=controlpar[9] && TIME<controlpar[10]) {
        n3 *= controlpar[11];
        double par[2];
        par[0] = controlpar[11]; par[1] = 0;
        incubator_update(conn3,update,par);
      }
      // Adult reduction
      if (TIME>=controlpar[12] && TIME<controlpar[13]) {
        n4fj *= controlpar[14];
        n4f *= controlpar[14];
        double par[2];
        par[0] = controlpar[14]; par[1] = 0;
        incubator_update(conn4,update,par);
      }
    }
    colT[TIME] = TIME;
    coln0[TIME] = n0;
    coln1[TIME] = n1;
    coln2[TIME] = n2;
    coln3[TIME] = n3;
    coln4fj[TIME] = n4fj;
    coln4f[TIME] = n4f;
    colK[TIME] = K;
    cold4[TIME] = d4;
    cold4s[TIME] = d4s;
    colF4[TIME] = F4;
    colegg[TIME] = egg;
    coldiap[TIME] = diap;
    if (isnan(n0) || isnan(n10) || isnan(n1) || isnan(n2) || isnan(n3) || isnan(n4fj) || isnan(n4f) || isnan(nBS) || isnan(K) || isnan(d4) || isnan(d4s) || isnan(F4) || isnan(egg) || isnan(percent_strong) || isnan(diap)) {
      //
      printf("ERROR_NAN: %d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",TIME,n0,n10,n1,n2,n3,n4fj,n4f,nBS,K,d4,d4s,F4,egg,percent_strong,diap);
      printf("ERROR_PAR: ");
      int i;
      for (i=0; i<NumParAea; i++)
        printf("%s%g",i==0?"":",",param[i]);
      printf("\n");
      // incubator_print(conn4);
      //
      (*success) = 0;
      goto endall;
    }
  }
  (*success) = 1;
  //
 endall:
  //
  empty_incubators();
  gamma_dist_check();
  //
  //printf("Times: %g %g %g\n",timebefore,timeof,timeafter);
}

// ---------------------------------------------------------------------------


