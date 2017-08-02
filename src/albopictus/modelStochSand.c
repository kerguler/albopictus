#ifdef R#define modelHeader#include "R.h"#include "Rmath.h"#else#include "math.h"#endif#include <time.h>#include <stdio.h>#include <stdlib.h>#include <string.h>#include <gsl/gsl_rng.h>#include <gsl/gsl_randist.h>#include "gamma.h"#include "spop.h"#include "ran_gen.h"#define max(a,b) ((a)>(b)?(a):(b))#define min(a,b) ((a)<(b)?(a):(b))#define ZERO    0#define MAX_DEV 1000extern gsl_rng *RAND_GSL;// --------------------------------------------// Gamma distribution// --------------------------------------------char gamma_mode = MODE_GAMMA_HASH;void set_gamma_mode(char mode) {  gamma_mode = mode;}// --------------------------------------------// Implementation of the dynamics// --------------------------------------------double *param;#define NumPar      35#define NumMet      7#define init_E      0#define init_L      1#define init_P      2#define init_A      3#define param_dE_1  4#define param_dE_2  5#define param_dE_3  6#define param_dE_4  7#define param_dL_1  8#define param_dL_2  9#define param_dL_3  10#define param_dL_4  11#define param_dP_1  12#define param_dP_2  13#define param_dP_3  14#define param_dP_4  15#define param_dA_1  16#define param_pE_1  17#define param_pE_2  18#define param_pE_3  19#define param_pL_1  20#define param_pL_2  21#define param_pL_3  22#define param_pP_1  23#define param_pP_2  24#define param_pP_3  25#define param_pA_1  26#define param_pA_2  27#define param_pA_3  28#define param_FA_1  29#define param_ovipos 30#define param_female 31#define param_PP_ta 32#define param_PP_hu 33#define param_eta   34void numparModel(int *np, int *nm) {  *np = NumPar;  *nm = NumMet;}void param_model(char **names, double *par) {  char temp[NumMet+NumPar][256] = {    "colnvE","colnvL","colnvP","colnvAf","colnvAfcap","colnvAm","colnvAmcap",    "init_E","init_L","init_P","init_A","param_dE_1","param_dE_2","param_dE_3","param_dE_4","param_dL_1","param_dL_2","param_dL_3","param_dL_4","param_dP_1","param_dP_2","param_dP_3","param_dP_4","param_dA_1","param_pE_1","param_pE_2","param_pE_3","param_pL_1","param_pL_2","param_pL_3","param_pP_1","param_pP_2","param_pP_3","param_pA_1","param_pA_2","param_pA_3","param_FA_1","param_ovipos","param_female","param_PP_ta","param_PP_hu","param_eta"  };  int i;  for (i=0; i<(NumMet+NumPar); i++)    names[i] = strdup(temp[i]);  // ---  //   par[init_E]      = 100; // initial number of egg  par[init_L]      = 100; // initial number of larva  par[init_P]      = 100; // initial number of pupa  par[init_A]      = 100; // initial number of adult    par[param_dE_1]  = 5.0; // min  par[param_dE_2]  = 35.0; // max  par[param_dE_3]  = 10.0; // Tmin  par[param_dE_4]  = 30.0; // Tmax    par[param_dL_1]  = 18.0;  par[param_dL_2]  = 30.0;  par[param_dL_3]  = 10.0;  par[param_dL_4]  = 30.0;  par[param_dP_1]  = 5.0;  par[param_dP_2]  = 32.0;  par[param_dP_3]  = 10.0;  par[param_dP_4]  = 30.0;    par[param_dA_1]  = 6.0;  par[param_pE_1]  = 50.0; // max  par[param_pE_2]  = 10.0; // Tmin  par[param_pE_3]  = 30.0; // Tmax    par[param_pL_1]  = 50.0;  par[param_pL_2]  = 10.0;  par[param_pL_3]  = 30.0;    par[param_pP_1]  = 50.0;  par[param_pP_2]  = 10.0;  par[param_pP_3]  = 30.0;    par[param_pA_1]  = 5.0;  par[param_pA_2]  = 10.0;  par[param_pA_3]  = 30.0;    par[param_FA_1]  = 80.0; // eggs per female in 4-8 days  par[param_ovipos] = 0.5; // surviving oviposition    par[param_female] = 0.5;    par[param_PP_ta] = 19.81;  par[param_PP_hu] = 0.7;    par[param_eta]   = 1.0;}// --------------------------------------------------------------------#define dsigL(x, a2) (max(0.0, min(1.0, 1.0/(1.0+exp(100.0*((a2)-(x)))))))#define dsig(x, a1, a2, a3) (max(0.0, min(0.999, (a1)/((1.0+exp((a2)-(x)))*(1.0+exp((x)-(a3)))))))#define dsigi(x, a1l, a1h, a2, a3) (max((a1l), min((a1h), ((a1h)-((a1h)-(a1l))/((1.0+exp((a2)-(x)))*(1.0+exp((x)-(a3))))))))int fun_capture(int    *nvA,                double *par) {  return gsl_ran_binomial(RAND_GSL,par[param_eta],*nvA);}int fun_female(int    *nvP,            double *par) {  return gsl_ran_binomial(RAND_GSL,par[param_female],*nvP);}void fun_devsur(double temp,                double hum,                double *par,                double *vec) {  // Diapause flag  char diapause = temp < par[param_PP_ta];  // Calculate development and survival times  /* d1 */  vec[0] = dsigi(temp,par[param_dE_1],par[param_dE_2],par[param_dE_3],par[param_dE_4]);  /* d2 */  vec[1] = diapause ? MAX_DEV : dsigi(temp,par[param_dL_1],par[param_dL_2],par[param_dL_3],par[param_dL_4]);  /* d3 */  vec[2] = dsigi(temp,par[param_dP_1],par[param_dP_2],par[param_dP_3],par[param_dP_4]);  /* d4 */  vec[3] = par[param_dA_1];  //  /* p1 */  vec[4] = 1.0-dsig(temp,1.0-(1.0/(1.0+par[param_pE_1])),par[param_pE_2],par[param_pE_3]);  /* p2 */  vec[5] = 1.0-dsig(temp,1.0-(1.0/(1.0+par[param_pL_1])),par[param_pL_2],par[param_pL_3]);  /* p3 */  vec[6] = 1.0-dsig(temp,1.0-(1.0/(1.0+par[param_pP_1])),par[param_pP_2],par[param_pP_3]);  /* p4 */  vec[7] = 1.0-dsig(temp,1.0-(1.0/(1.0+par[param_pA_1])),par[param_pA_2],par[param_pA_3]);  //  /* FA */  vec[8] = par[param_FA_1];  //  /* p1hu */  vec[9] = 1.0-dsigL(hum,par[param_PP_hu]);}// --------------------------------------------------------------------void calculateSand(double mean_air_temp,                   double daily_humidity,                   double daily_capture,                   spop   *vecE,                   spop   *vecL,                   spop   *vecP,                   spop   *vecAf,                   int    *nvE,                   int    *nvL,                   int    *nvP,                   int    *nvAf,                   int    *nvAfcap,                   int    *nvAm,                   int    *nvAmcap) {  // Calculate development and survival times  double vec[10];  fun_devsur(mean_air_temp,             daily_humidity,             param,             vec);  double d1=vec[0],    d2=vec[1],    d3=vec[2],    d4=vec[3],    p1=vec[4],    p2=vec[5],    p3=vec[6],    p4=vec[7],    FA=vec[8],    p1hu=vec[9];  // Check the results of the development processes and survival  int nvecE = spop_survive((*vecE),                           d1, sqrt(d1), // development                           -(p1*p1hu), 0, // survival (temperature and humidity)                           gamma_mode);  int nvecL = spop_survive((*vecL),                           d2, sqrt(d2), // development                           -p2, 0, // survival                           gamma_mode);  int nvecP = spop_survive((*vecP),                           d3, sqrt(d3), // development                           -p3, 0, // survival                           gamma_mode);  int nvecAf = spop_survive((*vecAf),                           d4, sqrt(d4), // development (oviposition)                           -p4, 0, // survival                           gamma_mode);  int tmp;  double tmpd = FA*nvecAf;  int nFA = (tmpd>1000.0) ? max(0.0,round(tmpd + gsl_ran_gaussian(RAND_GSL,sqrt(tmpd)))) : gsl_ran_poisson(RAND_GSL,tmpd);  // Newly emerged adult females  int newf = fun_female(&nvecP,param);  //  spop_add((*vecE),0,nFA);  spop_add((*vecL),0,nvecE);  spop_add((*vecP),0,nvecL);  //  tmp = gsl_ran_binomial(RAND_GSL,param[param_ovipos],nvecAf); // females surviving oviposition  spop_add((*vecAf),0,newf+tmp); // newly emerged ones and  // subsequent gonotrophic cycle of the oviposited ones  //  char capture = daily_capture > 0.5;  //  if (capture) {    tmp = spop_kill((*vecAf),                    param[param_eta]); // probability of capture    (*nvAfcap) = tmp;  } else {    (*nvAfcap) = 0;  }  //  tmp = gsl_ran_binomial(RAND_GSL,p4,(*nvAm));  (*nvAm) += nvecP - tmp - newf;  //  if (capture) {    tmp = fun_capture(nvAm,param);    (*nvAmcap) = tmp;    (*nvAm) -= tmp;  } else {    (*nvAmcap) = 0;  }  //  (*nvE) = (*vecE)->popsize;  (*nvL) = (*vecL)->popsize;  (*nvP) = (*vecP)->popsize;  (*nvAf) = (*vecAf)->popsize;  //  // printf("%d %d %d %d\n",(*nvE),(*nvL),(*nvP),(*nvAf));}// --------------------------------------------// sim_model// --------------------------------------------void sim_model(double   *envar,               double     *par,               int     *finalT,               int    *control,               double  *result,               int    *success) {  param = par;  double *controlpar = 0;  if ((*control)) {    controlpar = param + NumPar;  }    double *mean_air_temp = envar + 0*(*finalT);;  double *daily_humidity = envar + 1*(*finalT);;  double *daily_capture = envar + 2*(*finalT);;  double *colT   = result + 0*(*finalT);;  double *colnvE = result + 1*(*finalT);;  double *colnvL = result + 2*(*finalT);;  double *colnvP = result + 3*(*finalT);;  double *colnvAf = result + 4*(*finalT);;  double *colnvAfcap = result + 5*(*finalT);;  double *colnvAm = result + 6*(*finalT);;  double *colnvAmcap = result + 7*(*finalT);;  int TIME = 0;  (*success) = 2;  // Set initial conditions  int nvE = round(param[init_E]);  int nvL = round(param[init_L]);  int nvP = round(param[init_P]);  int nvA = round(param[init_A]);  int nvAf = fun_female(&nvA,param);  int nvAm = nvA - nvAf;  //  int nvAfcap;  int nvAmcap;  if (daily_capture[0] > 0.5) {    nvAfcap = fun_capture(&nvAf,param);    nvAmcap = fun_capture(&nvAm,param);    nvAf -= nvAfcap;    nvAm -= nvAmcap;  }  //  spop pop_vecE = spop_init();  spop pop_vecL = spop_init();  spop pop_vecP = spop_init();  spop pop_vecAf = spop_init();  if (nvE>0) spop_add(pop_vecE,0,nvE);  if (nvL>0) spop_add(pop_vecL,0,nvL);  if (nvP>0) spop_add(pop_vecP,0,nvP);  if (nvAf>0) spop_add(pop_vecAf,0,nvAf);  // Record state  colT[TIME] = (double)TIME;  colnvE[TIME] = (double)nvE;  colnvL[TIME] = (double)nvL;  colnvP[TIME] = (double)nvP;  colnvAf[TIME] = (double)nvAf;  colnvAfcap[TIME] = (double)nvAfcap;  colnvAm[TIME] = (double)nvAm;  colnvAmcap[TIME] = (double)nvAmcap;  //  for (TIME=1; TIME<(*finalT); TIME++) {    // Take a step    calculateSand(mean_air_temp[TIME-1],                  daily_humidity[TIME-1],                  daily_capture[TIME], // Capture using the previous day's results (i.e. traps stay overnight)                  &pop_vecE,                  &pop_vecL,                  &pop_vecP,                  &pop_vecAf,                  &nvE,                  &nvL,                  &nvP,                  &nvAf,                  &nvAfcap,                  &nvAm,                  &nvAmcap);    // Record state    colT[TIME] = (double)TIME;    colnvE[TIME] = (double)nvE;    colnvL[TIME] = (double)nvL;    colnvP[TIME] = (double)nvP;    colnvAf[TIME] = (double)nvAf;    colnvAfcap[TIME] = (double)nvAfcap;    colnvAm[TIME] = (double)nvAm;    colnvAmcap[TIME] = (double)nvAmcap;    if (isnan(nvE) || isnan(nvL) || isnan(nvP) || isnan(nvAf) || isnan(nvAfcap) || isnan(nvAm) || isnan(nvAmcap)) {      (*success) = 0;      goto endall;    }  }  // endall:  //  (*success) = ((*success)==0) ? 0 : 1;  spop_destroy(&pop_vecE);  spop_destroy(&pop_vecL);  spop_destroy(&pop_vecP);  spop_destroy(&pop_vecAf);  //  gamma_dist_check();}// ---------------------------------------------------------------------------