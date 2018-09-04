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

#define NumParAea      43
#define NumMetAea      11

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

#define alpha_dp_egg      29
#define alpha_dp_thr      30
#define alpha_ta_thr      31
#define alpha_tq_thr      32

#define alpha_tbm_1       33
#define alpha_tbm_2       34
#define alpha_tbm_3       35

#define alpha_p0_1        36
#define alpha_p0_2        37

#define alpha_n23_1       38
#define alpha_n23_2       39
#define alpha_n23_3       40
#define alpha_n23_4       41
#define alpha_n23_5       42

#define alpha_ta_thr_std  (0.5)
#define alpha_dp_thr_std  (10.0/24.0/60.0)
#define alpha_tq_thr_std  (0.5)

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

void calculate(double *photoperiod,
               double *mean_air_temp,
               double *daily_precipitation,
               double *popdens,
               double *param,
               spop   *conn10,
               spop   *conn1,
               spop   *conn2,
               spop   *conn3,
               spop   *conn4j,
               spop   *conn4,
               double *n0,
               double *n10,
               double *n1,
               double *nh,
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
               int    TIME) {
  // ---------------------
  // modelDelayAalbopictus
  // ---------------------

  double deltaT = param[alpha_deltaT];

  double Ta = mean_air_temp[TIME];
  
  double Tw = Ta + deltaT;

  /*
   * Update the number of breeding sites
   */
  (*nBS) = param[alpha_BS_pdens] * (*popdens) +
    param[alpha_BS_dprec] * daily_precipitation[TIME] +
    param[alpha_BS_nevap] * (*nBS);
  (*K) = param[alpha_BS_nevap] == 1.0 ? (*nBS) / (TIME+1.0) : (*nBS) * (param[alpha_BS_nevap]-1.0) / (pow(param[alpha_BS_nevap], (TIME+1))-1.0);

  /*
   * Update survival and development rates/probabilities according to
   * the new environmental conditions and breeding sites
   */
  // Density of the immature stages
  // Assume uniform distribution across all breeding sites
  double n23dens = ((*n2) + (*n3)) / (*K);
  double densd = expd(n23dens,param[alpha_n23_surv]);
  double densdev = fpow(n23dens,param[alpha_n23_1],param[alpha_n23_2]*poly(Tw,param[alpha_n23_3],param[alpha_n23_4],param[alpha_n23_5]));
  // Fecundity
  double bigF4 = poly(Ta,param[alpha_F4_1],param[alpha_F4_2],param[alpha_F4_3]);
  
  // Egg survival (diapausing eggs)
  double p0_Ta = flin(Ta,param[alpha_p0_1],param[alpha_p0_2]);
  // Egg mortality (non-diapausing eggs)
  double p1_Tw = 1.0 - dsig(Tw,param[alpha_p1_1],param[alpha_p1_2],param[alpha_p1_3]);
  // Larval mortality
  double p2_Tw = 1.0 - dsig(Tw,param[alpha_p2_1],param[alpha_p2_2],param[alpha_p2_3])*densd;
  // Pupal mortality
  double p3_Tw = 1.0 - dsig(Tw,param[alpha_p3_1],param[alpha_p3_2],param[alpha_p3_3])*densd;
  
  // Egg development time
  double d1 = max(1.0, poly(Tw,param[alpha_d1_1],param[alpha_d1_2],param[alpha_d1_3]));
  // Larval development time
  double d2 = max(1.0, poly(Tw,param[alpha_d2_1],param[alpha_d2_2],param[alpha_d2_3])*densdev);
  // Pupal development time
  double d3 = max(1.0, poly(Tw,param[alpha_d3_1],param[alpha_d3_2],param[alpha_d3_3])*densdev);
  // Time to first blood meal
  double alpha_blood = poly(Ta,param[alpha_tbm_1],param[alpha_tbm_2],param[alpha_tbm_3]);
  // Adult lifetime (from emergence)
  double dd4 = dsig2(Ta,param[alpha_d4_1],param[alpha_d4_2],param[alpha_d4_3]);
  double dd4s = gamma_matrix_sd*dd4;

  // Diapause and quiescence
  // quie -> larvae
  double fracHatch = Tthr(Ta,param[alpha_ta_thr],alpha_ta_thr_std);
  // neweggs -> n10
  double fracDiap = pdiap(photoperiod[TIME],param[alpha_dp_thr],alpha_dp_thr_std,fracHatch);
  // n0 -> quie
  double fracQuie = 1.0 - Tthr(Ta,param[alpha_tq_thr],alpha_tq_thr_std);
  
  /*
   * Survival first, and then, development.
   * Update development times and survival proportions
   * of the immature stages and juvenile adults
   */
  //
  // Normal and tagged eggs
  spop_iterate((*conn1),
               0,
               d1, 0, // development (fixed-length)
               0,
               p1_Tw, // death (fixed daily rate)
               0, 0,
               0,
               0);
  spop_iterate((*conn10),
               0,
               d1, 0, // development (fixed-length)
               0,
               p1_Tw, // death (fixed daily rate)
               0, 0,
               0,
               0);
  double quie = (1.0 - p1_Tw) * (*nh); // Number of surviving quiescent eggs
  //
  // Larvae
  spop_iterate((*conn2),
               0,
               d2, 0, // development (fixed-length)
               0,
               p2_Tw, // death (fixed daily rate)
               0, 0,
               0,
               0);
  //
  // Pupae
  spop_iterate((*conn3),
               0,
               d3, 0, // development (fixed-length)
               0,
               p3_Tw, // death (fixed daily rate)
               0, 0,
               0,
               0);
  //
  // Adult naive females
  spop_iterate((*conn4j),
               0,
               alpha_blood, 0, // development (fixed-length)
               0,
               0,
               dd4, dd4s, // death (gamma-distributed)
               0,
               0);
  //
  // Adult mature females
  spop_iterate((*conn4),
               0,
               0, 0, // development (no development)
               0,
               0,
               dd4, dd4s, // death (gamma-distributed)
               0,
               0);
  //
  // Diapausing egg survival
  // Tagged eggs always become diapausing eggs
  double dp_eggs = p0_Ta*(*n0) + (*conn10)->developed.d;
  // Adults survived to produce eggs
  double n4_reproduce = (*conn4)->size.d;
  //
  // Perform state transformations
  spop_popadd((*conn4), (*conn4j)->devtable);
  spop_add((*conn4j), 0, 0, 0, 0.5 * (*conn3)->developed.d);
  spop_add((*conn3), 0, 0, 0, (*conn2)->developed.d);
  spop_add((*conn2), 0, 0, 0, (*conn1)->developed.d + (fracHatch * quie));
  quie = (1.0 - fracHatch) * quie;
  //  
  // Lay eggs
  double neweggs = bigF4*n4_reproduce; // Total number of eggs that will be laid that day
  //
  spop_add((*conn10),0,0,0,fracDiap * neweggs); // Tagged eggs
  spop_add((*conn1),0,0,0,(1.0-fracDiap) * neweggs); // Normal eggs
  //
  quie += fracQuie * dp_eggs;
  dp_eggs = (1.0 - fracQuie) * dp_eggs;
  //
  // Update all stages and state variables
  (*n0) = dp_eggs;
  (*n10) = (*conn10)->size.d;
  (*n1) = (*conn1)->size.d;
  (*nh) = quie;
  (*n2) = (*conn2)->size.d;
  (*n3) = (*conn3)->size.d;
  (*n4fj) = (*conn4j)->size.d;
  (*n4f) = (*conn4)->size.d;
  // Update indicators
  (*d4) = dd4;
  (*d4s) = dd4s;
  (*F4) = bigF4;
  (*egg) = neweggs;
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
    "coln0","coln1","coln2","coln3","coln4fj","coln4f","colK","cold4","cold4s","colF4","colegg",
    "p1.1","p1.2","p1.3","p2.1","p2.2","p2.3","p3.1","p3.2","p3.3","d4.1","d4.2","d4.3","F4.1","F4.2","F4.3","d1.1","d1.2","d1.3","d2.1","d2.2","d2.3","d3.1","d3.2","d3.3","n23.surv","deltaT","BS.pdens","BS.dprec","BS.nevap","PP.init","PP.thr","PP.ta.thr","PP.tq.thr","tbm.1","tbm.2","tbm.3","p0.1","p0.2","n23.1","n23.2","n23.3","n23.4","n23.5"
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

  param[alpha_dp_egg] = 1.0;
  param[alpha_dp_thr] = 0.5;
  param[alpha_ta_thr] = 21.0;
  param[alpha_tq_thr] = 4.0;

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

  int TIME = 0;
  (*success) = 2;
  // Set initial conditions
  double n0 = param[alpha_dp_egg];
  double n10 = 0.0;
  double n1 = 0.0;
  double nh = 0.0;
  double n2 = 0.0;
  double n3 = 0.0;
  double n4fj = 0.0;
  double n4f = 0.0;
  double nBS = 0.0;
  double K = 0.0;
  double d4 = 0.0;
  double d4s = 0.0;
  double F4 = 0.0;
  double egg = 0.0;
  //
  spop conn10 = spop_init(0,gamma_mode);
  spop conn1 = spop_init(0,gamma_mode);
  spop conn2 = spop_init(0,gamma_mode);
  spop conn3 = spop_init(0,gamma_mode);
  spop conn4 = spop_init(0,gamma_mode);
  spop conn4j = spop_init(0,gamma_mode);
  if (n10 > DPOP_EPS) spop_add(conn10,0,0,0,n10);
  if (n1 > DPOP_EPS) spop_add(conn1,0,0,0,n1);
  if (n2 > DPOP_EPS) spop_add(conn2,0,0,0,n2);
  if (n3 > DPOP_EPS) spop_add(conn3,0,0,0,n3);
  if (n4fj > DPOP_EPS) spop_add(conn4j,0,0,0,n4fj);
  if (n4f > DPOP_EPS) spop_add(conn4,0,0,0,n4f);
  // Record state
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
  //
  for (TIME=1; TIME<(*finalT); TIME++) {
    // Take a step
    calculate(photoperiod,
              mean_air_temp,
              daily_precipitation,
              popdens,
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
              &K,
              &d4,
              &d4s,
              &F4,
              &egg,
              TIME);
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
        nh *= controlpar[5];
        spop_iterate(conn1,
                     0,
                     0, 0,
                     0,
                     1.0 - controlpar[5],
                     0, 0,
                     0,
                     1);
        spop_iterate(conn10,
                     0,
                     0, 0,
                     0,
                     1.0 - controlpar[5],
                     0, 0,
                     0,
                     1);
      }
      // Larva reduction
      if (TIME>=controlpar[6] && TIME<controlpar[7]) {
        n2 *= controlpar[8];
        spop_iterate(conn2,
                     0,
                     0, 0,
                     0,
                     1.0 - controlpar[8],
                     0, 0,
                     0,
                     1);
      }
      // Pupa reduction
      if (TIME>=controlpar[9] && TIME<controlpar[10]) {
        n3 *= controlpar[11];
        spop_iterate(conn3,
                     0,
                     0, 0,
                     0,
                     1.0 - controlpar[11],
                     0, 0,
                     0,
                     1);
      }
      // Adult reduction
      if (TIME>=controlpar[12] && TIME<controlpar[13]) {
        n4fj *= controlpar[14];
        n4f *= controlpar[14];
        spop_iterate(conn4j,
                     0,
                     0, 0,
                     0,
                     1.0 - controlpar[14],
                     0, 0,
                     0,
                     1);
        spop_iterate(conn4,
                     0,
                     0, 0,
                     0,
                     1.0 - controlpar[14],
                     0, 0,
                     0,
                     1);
      }
    }
    // Record state
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
    if (isnan(n0) || isnan(n10) || isnan(n1) || isnan(nh) || isnan(n2) || isnan(n3) || isnan(n4fj) || isnan(n4f) || isnan(nBS) || isnan(K) || isnan(d4) || isnan(d4s) || isnan(F4) || isnan(egg)) {
      //
      printf("ERROR_NAN: %d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",TIME,n0,n10,n1,nh,n2,n3,n4fj,n4f,nBS,K,d4,d4s,F4,egg);
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


