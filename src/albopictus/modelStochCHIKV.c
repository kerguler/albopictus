#ifdef R
#define modelHeader
#include "R.h"
#include "Rmath.h"
#else
#include "math.h"
#endif
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "gamma.h"
#include "spop01.h"
#include "ran_gen.h"

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

#define ZERO    0
#define CHECK(x) ((x)<0)

#define NumParChk      31
#define NumMetChk      13

extern gsl_rng *RAND_GSL;

// --------------------------------------------
// Simulation mode
// --------------------------------------------

char quarantine = 0; // quarantine symptomatic infections
char sim_mode = 1; // 0: check for secondary infection, 1: full simulation
void set_sim_mode(char mode) {
  sim_mode = mode;
}

// --------------------------------------------
// Gamma distribution
// --------------------------------------------

char gamma_mode = MODE_GAMMA_HASH;
void set_gamma_mode(char mode) {
  gamma_mode = mode;
}

// --------------------------------------------
// Implementation of the dynamics
// --------------------------------------------

#define init_S    0
#define init_E    1
#define init_I    2
#define init_Ic   3
#define init_C    4
#define init_R    5
#define init_vAe  6
#define init_vAi  7
#define scale     8
#define gamma     9
#define gamma_sd  10
#define omega     11
#define omega_sd  12
#define susceptH  13
#define prob_chronic 14
#define prob_C_to_Ic 15
#define prob_C_to_R 16
#define susceptV  17
#define vec_omega 18
#define vec_omega_sd   19
#define ovitraps  20
#define introduce_time 21
#define introduce_location 22
#define introduce_infectious 23
#define introduce_ivector 24
#define fixed_bites 25
#define report    26
#define qmove     27
#define connect   28
#define Khum      29
#define Kvec      30

volatile clock_t start = 0, diff;
void time2here(int mark) {
  diff = clock() - start;
  start = clock();
  printf("timeit\t%d\t%u\n",mark,(unsigned int)diff);
}

double *param;
int *total_human;
double *tH;
double param_prob_C_to_R;
double param_prob_C_to_Ic;

void calc_tH(double *tprobs, int *numreg) {
  int region, regionB;
  for (region=0; region<(*numreg); region++) {
    tH[region] = 0;
    for (regionB=0; regionB<(*numreg); regionB++) {
      tH[region] += tprobs[regionB*(*numreg)+region]*total_human[regionB];
    }
  }
}

void calculateAeCHIKVinit(double daily_vector,
                          double      popdens,
                          spop          *humE,
                          spop          *humI,
                          spop         *humIc,
                          spop         *vecAe,
                          spop         *vecAi,
                          int             *nS,
                          int             *nE,
                          int             *nI,
                          int            *nIc,
                          int             *nC,
                          int             *nR,
                          int           *nvAs,
                          int           *nvAe,
                          int           *nvAi,
                          int          *nInew) {
  (*humE) = spop_init();
  (*humI) = spop_init();
  (*humIc) = spop_init();
  (*vecAe) = spop_init();
  (*vecAi) = spop_init();

  (*nS) = round(popdens*param[init_S]*param[scale]);
  (*nE) = round(popdens*param[init_E]*param[scale]);
  (*nI) = round(popdens*param[init_I]*param[scale]);
  (*nIc) = round(popdens*param[init_Ic]*param[scale]);
  (*nC) = round(popdens*param[init_C]*param[scale]);
  (*nR) = round(popdens*param[init_R]*param[scale]);
  (*nvAe) = round(param[init_vAe]*param[ovitraps]*param[scale]);
  (*nvAi) = round(param[init_vAi]*param[ovitraps]*param[scale]);
  (*nvAs) = max(0, round(daily_vector*param[ovitraps]*param[scale])-(*nvAe)-(*nvAi));
  (*nInew) = 0;
  
  if ((*nE)>0) spop_add((*humE),0,(*nE));
  if ((*nI)>0) spop_add((*humI),0,(*nI));
  if ((*nIc)>0) spop_add((*humIc),0,(*nIc));
  if ((*nvAe)>0) spop_add((*vecAe),0,(*nvAe));
  if ((*nvAi)>0) spop_add((*vecAi),0,(*nvAi));

  // Re-scaling is done for the inference
  param_prob_C_to_R = param[prob_C_to_R]/365.0;
  param_prob_C_to_Ic = param[prob_C_to_Ic]/365.0;
  // ---
}

void calculateAeCHIKVtransmission(int    numreg,
                                  int    regionA,
                                  double *tprobs,
                                  double *fecundity,
                                  int    *nS,
                                  int    *nE,
                                  int    *nI,
                                  int    *nIc,
                                  int    *nC,
                                  int    *nR,
                                  int    *nvAs,
                                  int    *nvAe,
                                  int    *nvAi,
                                  spop   *humE,
                                  spop   *vecAe,
                                  double *npV,
                                  double *npH,
                                  double *nbite) {
  double total_vector = nvAs[regionA]+nvAe[regionA]+nvAi[regionA];
  double daily_bites;
  char discrete = 0;
  double tmp;
  int regionB, regionT;
  //
  if (quarantine) // quarantine (variable human population size)
    calc_tH(tprobs,&numreg);
  //
  // time2here(50);
  double tHI = 0, tHS = 0;
  double *tHSc = (double *)malloc(numreg*sizeof(double));
  for (regionB=0; regionB<numreg; regionB++) {
    tmp = tprobs[regionB*numreg+regionA];
    tHI += tmp*(nI[regionB]+nIc[regionB]); // sum
    tHS += tmp*nS[regionB]; // sum
    tHSc[regionB] = tHS; // cumulative
  }
  //
  // time2here(51);
  if (param[fixed_bites]) {
    // Biting rate independent of humans
    tmp = param[Kvec]*total_vector;
  } else {
    // Fecundity-dependent bite-rate
    tmp = 2.0*fecundity[regionA]*param[Kvec]*param[Khum]*total_vector*tH[regionA]/(fecundity[regionA]*param[Kvec]*total_vector + param[Khum]*tH[regionA]); // Harmonic mean of two rates
  }
  if (tmp>1000.0) {
    //printf("tmp=%g,sqrt(tmp)=%g,",tmp,sqrt(tmp));
    daily_bites = max(0,tmp + gsl_ran_gaussian(RAND_GSL,sqrt(tmp)));
    //printf("daily_bites=%g\n",daily_bites);
  } else {
    discrete = 1;
    daily_bites = gsl_ran_poisson(RAND_GSL,tmp);
  }
  nbite[regionA] = daily_bites;
  // daily_bites includes bites to all regions
  //
  double probH = 0,
    probV = 0,
    probHV,
    prob;
  double counter_nbites = 0,
    counter_rng;
  //
  // time2here(52);
  while (counter_nbites < daily_bites) {
    probH = param[susceptH]*nvAi[regionA]*tHS/(total_vector*tH[regionA]);
    probV = param[susceptV]*nvAs[regionA]*tHI/(total_vector*tH[regionA]);
    if (counter_nbites==0) {
      npV[regionA] = probV;
      npH[regionA] = probH;
      /*
      printf("%d, %g,%g, %g, %g,%g, %g,%d,%g, %g,%d,%g\n",
             regionA,
             tH[regionA],total_vector,
             daily_bites,
             probH,probV,
             param[susceptH],nvAi[regionA],tHS,
             param[susceptV],nvAs[regionA],tHI);
      */
    }
    probHV = probH + probV;
    if (probHV <= 0) {
      break;
    }
    // Number of bites until next S->E (probH) or Vs->Vi (probV)
    if (discrete) {
      counter_rng = gsl_ran_geometric(RAND_GSL,probHV);
    } else {
      counter_rng = rng_exponential(probHV);
    }
    if (counter_rng <= 0) {
      break;
    }
    counter_nbites += counter_rng; // nbites represent the successful next bite
    if (counter_nbites > daily_bites) {
      break;
    }
    // probH or probV?
    prob = gsl_rng_uniform(RAND_GSL);
    if (prob < probH/probHV) { // vector->human transmission occured (probH: S->E)
      // Which region?
      tmp = tHS*gsl_rng_uniform(RAND_GSL);
      for (regionB=0; regionB<numreg; regionB++)
        if (tmp<=tHSc[regionB]) break;
      spop_add(humE[regionB],0,1);
      nE[regionB]++;
      nS[regionB]--;
      //
      tmp = tprobs[regionB*numreg+regionA];
      tHS -= tmp;
      for (regionT=regionB; regionT<numreg; regionT++) {
        tHSc[regionT] -= tmp;
      }
      //
      // printf("%d->%d, %d S -> %d E\n",regionA,regionB, nS[regionB],nE[regionB]);
      //
    } else { // human->vector transmission occured (probV: Vs->Vi)
      spop_add(vecAe[regionA],0,1);
      nvAe[regionA]++;
      nvAs[regionA]--;
    }
  }
  //
  free(tHSc);
  // time2here(53);
}

void calculateAeCHIKVexperiment(spop   *vecAi,
                                spop   *humI,
                                int    *nvAi,
                                int    *nI,
                                int    *nInew,
                                int    *nS) {
  int infect = round(param[introduce_infectious]);
  if (param[introduce_ivector]) {
    // Introduce an infectious vector
    spop_add((*vecAi),0,infect);
    (*nvAi) += infect;
  } else if ((*nS)>=infect) {
    // Introduce an infectious human
    spop_add((*humI),0,infect);
    (*nI) += infect;
    (*nInew) += infect;
    (*nS) -= infect;
  }
}

void calculateAeCHIKV(double daily_vector,
                      double vector_survival_times,
                      double vector_survival_times_sd,
                      double daily_vector_fecundity,
                      spop   *humE,
                      spop   *humI,
                      spop   *humIc,
                      spop   *vecAe,
                      spop   *vecAi,
                      int    *nS,
                      int    *nE,
                      int    *nI,
                      int    *nIc,
                      int    *nC,
                      int    *nR,
                      int    *nvAs,
                      int    *nvAe,
                      int    *nvAi,
                      double *npV,
                      double *npH,
                      double *nbite,
                      int    *nInew,
                      int    *totalH) {
  // Check the results of the development processes (and survival)
  // printf("B: %d, %d, %d, %d, %d\n",(*humE)->popsize,(*humI)->popsize,(*humIc)->popsize,(*vecAe)->popsize,(*vecAi)->popsize);
  int nhumE = spop_survive((*humE),
                           param[omega], param[omega_sd], // development
                           0.0, 0.0, // survival (immortal humans)
                           gamma_mode);
  int nhumI = spop_survive((*humI),
                           param[gamma], param[gamma_sd], // development
                           0.0, 0.0, // survival (immortal humans)
                           gamma_mode);
  int nhumIc = spop_survive((*humIc),
                            param[gamma], param[gamma_sd], // development
                            0.0, 0.0, // survival (immortal humans)
                            gamma_mode);
  int nvecAe = spop_survive((*vecAe),
                            param[vec_omega], param[vec_omega_sd], // development
                            vector_survival_times, vector_survival_times_sd, // survival
                            gamma_mode);
  int nvecAi = spop_survive((*vecAi),
                            0.0, 0.0, // development (persistent infectious stage)
                            vector_survival_times, vector_survival_times_sd, // survival
                            gamma_mode);
  if (nvecAi != 0) {
    printf("ERROR: Development observed for infectious vector!\n");
    exit(1);
  }
  // printf("A: %d, %d, %d, %d, %d\n",nhumE,nhumI,nhumIc,nvecAe,nvecAi);

  // Decide the fate of acute infections
  int nhumIC = gsl_ran_binomial(RAND_GSL,param[prob_chronic],nhumI); // end up in the chronic stage
  int nhumIR = nhumI - nhumIC; // recover completely

  // Update chronic stage recovery and relapse
  // This is temperature-independent relapse
  double prob_R_Ic = min(1.0, param_prob_C_to_R + param_prob_C_to_Ic);
  int C_to_R_Ic = gsl_ran_binomial(RAND_GSL,prob_R_Ic,(*nC)); // total exit from the chronic stage
  int nhumCR = gsl_ran_binomial(RAND_GSL,param_prob_C_to_R/prob_R_Ic,C_to_R_Ic); // recovery
  int nhumCIc = C_to_R_Ic - nhumCR; // relapsed infections

  // Perform state transformations
  if (!quarantine) // not quarantine
    spop_add((*humI),0,nhumE);
  spop_add((*humIc),0,nhumCIc);
  spop_add((*vecAi),0,nvecAe); // newly infected vector population
  // Number of newly infectious symptomatic cases which are reported
  (*nInew) = gsl_ran_binomial(RAND_GSL,param[report],nhumE);
  if (quarantine) { // quarantine
    spop_add((*humI),0,nhumE-(*nInew));
    (*totalH) -= (*nInew);
  }
  (*nE) = (*humE)->popsize; // exposed human population
  (*nI) = (*humI)->popsize; // infectious human population
  (*nIc) = (*humIc)->popsize; // human population with infectious relapse
  (*nC) += nhumIc + nhumIC - nhumCIc - nhumCR; // chronic human cases
  (*nR) += nhumIR + nhumCR; // recovered human population
  (*nvAi) = (*vecAi)->popsize; // infected vector population
  (*nvAe) = (*vecAe)->popsize; // exposed vector population
  // daily_vector from modelDelayAalbopictus is for the end of the day
  (*nvAs) = max(0, round(daily_vector*param[ovitraps]*param[scale])-(*nvAe)-(*nvAi));
  //
  // printf("B: %d %d %d\n",*nvAs,*nvAe,*nvAi);
}

// --------------------------------------------
// sim_model
// --------------------------------------------

void numparModel(int *np, int *nm) {
  *np = NumParChk;
  *nm = NumMetChk;
}

void param_model(char **names, double *par) {
  char temp[NumMetChk+NumParChk][256] = {
    "colnS","colnE","colnI","colnIc","colnC","colnR","colnvAs","colnvAe","colnvAi","colnpV","colnpH","colnbite","colnInew",
    "init_S","init_E","init_I","init_Ic","init_C","init_R","init_vAe","init_vAi","scale","gamma","gamma_sd","omega","omega_sd","susceptH","prob_chronic","prob_C_to_Ic","prob_C_to_R","susceptV","vec_omega","vec_omega_sd","ovitraps","introduce_time","introduce_location","introduce_infectious","introduce_ivector","fixed_bites","report","qmove","connect","Khum","Kvec"
  };
  int i;
  for (i=0; i<(NumMetChk+NumParChk); i++)
    names[i] = strdup(temp[i]);
  // ---
  // 
  par[init_S] = 1.0; // susceptible (% popdens)
  par[init_E] = 0; // latent (% popdens)
  par[init_I] = 0; // infectious (a/symptomatic) (% popdens)
  par[init_Ic] = 0; // infectious relapse (a/symptomatic) (% popdens)
  par[init_C] = 0; // chronic stage (% popdens)
  par[init_R] = 0; // recovered (% popdens)
  par[init_vAe] = 0; // adult vector (latent)
  par[init_vAi] = 0; // adult vector (infectious)
  par[scale] = 1.0; // simulation span (affects double -> int conversion)
  par[gamma] = 4.5; // infectious period in humans (2-7 days)
  par[gamma_sd] = 1.0; //
  par[omega] = 3.5; // latent period in humans (2-4 days)
  par[omega_sd] = 0.5; //
  par[susceptH] = 0.65; //0.8; // human susceptibility to infections (50-80%)
  par[prob_chronic] = 0.33; // fraction of infections entering into chronic stage
  par[prob_C_to_Ic] = 1.0; // daily probability of relapse (/365)
  par[prob_C_to_R] = 1.0; // daily probability of recovery from chronic stage (/365)
  //                           0 <= prob_C_to_R <= 1
  par[susceptV] = 0.85; //1.0; // vector susceptibility to infections (70-100%)
  par[vec_omega] = 3.0; // latent period in vectors (2-3 days)
  par[vec_omega_sd] = 0.1; //
  //
  par[ovitraps] = 0.32; // ovitrap / km^2 (~200 ovitraps per grid point, 25*25 km^2)
  // Also, each ovitrap samples 10% of the population in the vicinity (?)
  //
  par[introduce_time] = 0; // when to introduce an infection?
  par[introduce_location] = 0; // where to introduce the infection?
  par[introduce_infectious] = 0; // how many infected human/vector to be introduced?
  par[introduce_ivector] = 0; // introduce infectious vector? or human?
  //
  par[fixed_bites] = 0; // fixed or fecundity-dependent bite rates?
  //
  par[report] = 1.0; // fraction of cases reported (scale on the output)
  //
  par[qmove] = 1.0;
  par[connect] = 1.0;
  par[Khum] = 50.0;
  par[Kvec] = 0.25; // num. bites / mosquito / day (?)
}

// --------------------------------------------------------------------

void prepare_tprobs(int *numreg, double *ttprobs, double *tprobs) {
  int rA, rB, i=0, j=0;
  double sum, cen;
  for (rA=0; rA<(*numreg); rA++) {
    sum = 0;
    j = i;
    for (rB=0; rB<(*numreg); rB++) {
      if (rA==rB)
        cen = tprobs[i] = param[qmove]*ttprobs[i];
      else
        sum += tprobs[i] = pow(ttprobs[i],param[connect]);
      i++;
    }
    cen = 1.0 - cen;
    for (rB=0; rB<(*numreg); rB++) {
      if (rA!=rB)
        tprobs[j] *= cen / sum;
      // printf("%d %d %d %g %g\n",rA,rB,j,ttprobs[j],tprobs[j]);
      j++;
    }
  }
}

void sim_spread(double   *envar,
                double *ttprobs,
                double     *par,
                int     *finalT,
                int    *control,
                int     *numreg,
                double  *result,
                int    *success) {
  int r, tmp;
  param = par;
  double *controlpar = 0;
  if ((*control)) {
    controlpar = param + NumParChk;
  }
  //
  quarantine = 0;
  //
  double *tprobs = (double *)malloc((*numreg)*(*numreg)*sizeof(double));
  if (tprobs == NULL) {
    for (r=0; r<(*numreg); r++) success[r] = 0;
    return;
  }
  prepare_tprobs(numreg,ttprobs,tprobs);
  //
  total_human = (int *)malloc((*numreg)*sizeof(int));
  
  spop *pop_humE = (spop *)malloc((*numreg)*sizeof(spop));
  spop *pop_humI = (spop *)malloc((*numreg)*sizeof(spop));
  spop *pop_humIc = (spop *)malloc((*numreg)*sizeof(spop));
  spop *pop_vecAe = (spop *)malloc((*numreg)*sizeof(spop));
  spop *pop_vecAi = (spop *)malloc((*numreg)*sizeof(spop));

  char *experiment = (char *)calloc((*numreg),sizeof(char));

  int *nS = (int *)calloc((*numreg),sizeof(int));
  int *nE = (int *)calloc((*numreg),sizeof(int));
  int *nI = (int *)calloc((*numreg),sizeof(int));
  int *nIc = (int *)calloc((*numreg),sizeof(int));
  int *nC = (int *)calloc((*numreg),sizeof(int));
  int *nR = (int *)calloc((*numreg),sizeof(int));
  int *nvAs = (int *)calloc((*numreg),sizeof(int));
  int *nvAe = (int *)calloc((*numreg),sizeof(int));
  int *nvAi = (int *)calloc((*numreg),sizeof(int));
  double *npV = (double *)calloc((*numreg),sizeof(double));
  double *npH = (double *)calloc((*numreg),sizeof(double));
  double *nbite = (double *)calloc((*numreg),sizeof(double));
  int *nInew = (int *)calloc((*numreg),sizeof(int));

  double *environ;
  double *output;
  double *daily_vector;
  double *vector_survival_times;
  double *vector_survival_times_sd;
  double *daily_vector_fecundity;
  double *colT;
  double *colnS;
  double *colnE;
  double *colnI;
  double *colnIc;
  double *colnC;
  double *colnR;
  double *colnvAs;
  double *colnvAe;
  double *colnvAi;
  double *colnpV;
  double *colnpH;
  double *colnbite;
  double *colnInew;

  double *fecundity = (double *)malloc((*numreg)*sizeof(double));
  tH = (double *)malloc((*numreg)*sizeof(double));
  
  int TIME = 0;
  int region;
  for (region=0; region<(*numreg); region++) {
    if (region != (int)floor(param[introduce_location])) experiment[region] = 1;
    success[region] = 2;
  }
  // Set initial conditions
  /*
    WARNING!
    ========
    
    Information to decode environmental variables is hard coded here!
  */
  // time2here(0);
  for (region=0; region<(*numreg); region++) {
    calculateAeCHIKVinit(*(envar + region*(4*(*finalT)+1) + 0*(*finalT)), // +1 is for popdens,
                         *(envar + region*(4*(*finalT)+1) + 4*(*finalT)), // +1 is for popdens
                         &pop_humE[region],
                         &pop_humI[region],
                         &pop_humIc[region],
                         &pop_vecAe[region],
                         &pop_vecAi[region],
                         &nS[region],
                         &nE[region],
                         &nI[region],
                         &nIc[region],
                         &nC[region],
                         &nR[region],
                         &nvAs[region],
                         &nvAe[region],
                         &nvAi[region],
                         &nInew[region]);
    total_human[region] = nS[region]+nE[region]+nI[region]+nIc[region]+nC[region]+nR[region];
    //
    output   = result + region*14*(*finalT);
    colT     = output + 0*(*finalT);
    colnS    = output + 1*(*finalT);
    colnE    = output + 2*(*finalT);
    colnI    = output + 3*(*finalT);
    colnIc   = output + 4*(*finalT);
    colnC    = output + 5*(*finalT);
    colnR    = output + 6*(*finalT);
    colnvAs  = output + 7*(*finalT);
    colnvAe  = output + 8*(*finalT);
    colnvAi  = output + 9*(*finalT);
    colnpV   = output + 10*(*finalT);
    colnpH   = output + 11*(*finalT);
    colnbite = output + 12*(*finalT);
    colnInew = output + 13*(*finalT);
    colT[TIME] = (double)TIME;
    colnS[TIME] = (double)nS[region];
    colnE[TIME] = (double)nE[region];
    colnI[TIME] = (double)nI[region];
    colnIc[TIME] = (double)nIc[region];
    colnC[TIME] = (double)nC[region];
    colnR[TIME] = (double)nR[region];
    colnvAs[TIME] = (double)nvAs[region];
    colnvAe[TIME] = (double)nvAe[region];
    colnvAi[TIME] = (double)nvAi[region];
    colnpV[TIME] = npV[region];
    colnpH[TIME] = npH[region];
    colnbite[TIME] = nbite[region];
    colnInew[TIME] = (double)nInew[region];
  }
  //
  // time2here(1);
  calc_tH(tprobs,numreg);
  //
  for (TIME=1; TIME<(*finalT); TIME++) {
    // printf("TIME: %d\n",TIME);
    // time2here(2);
    for (region=0; region<(*numreg); region++) {
      // Take the first step
      /*
        WARNING!
        ========

        Information to decode environmental variables is hard coded here!
       */
      environ                   = envar + region*(4*(*finalT)+1); // +1 is for popdens
      daily_vector              = environ + 0*(*finalT);
      vector_survival_times     = environ + 1*(*finalT);
      vector_survival_times_sd  = environ + 2*(*finalT);
      daily_vector_fecundity    = environ + 3*(*finalT);
      
      calculateAeCHIKV(daily_vector[TIME-1],
                       vector_survival_times[TIME-1],
                       vector_survival_times_sd[TIME-1],
                       daily_vector_fecundity[TIME-1],
                       &pop_humE[region],
                       &pop_humI[region],
                       &pop_humIc[region],
                       &pop_vecAe[region],
                       &pop_vecAi[region],
                       &nS[region],
                       &nE[region],
                       &nI[region],
                       &nIc[region],
                       &nC[region],
                       &nR[region],
                       &nvAs[region],
                       &nvAe[region],
                       &nvAi[region],
                       &npV[region],
                       &npH[region],
                       &nbite[region],
                       &nInew[region],
                       &total_human[region]);
    }
    // Apply control measures
    if ((*control)) {
      // Adult reduction
      if (TIME>=controlpar[0] && TIME<controlpar[1]) {
        //
        for (region=0; region<(*numreg); region++) { // Apply to all regions at once!
          // printf("chikv before: %d %d %g %g %g   %g %d %d %d %d\n",region,TIME,controlpar[0],controlpar[1],controlpar[2],(envar + region*(4*(*finalT)+1))[TIME-1],nvAs[region],nvAe[region],nvAi[region],pop_vecAi[region]->popsize);
          //
          tmp = spop_survive(pop_vecAe[region],
                             0.0, 0.0, // no development
                             controlpar[2]-1.0, -1.0, // daily probability of death
                             gamma_mode);
          tmp += spop_survive(pop_vecAi[region],
                              0.0, 0.0, // no development
                              controlpar[2]-1.0, -1.0, // daily probability of death
                              gamma_mode);
          if (tmp) {
            printf("ERROR: Vector control failed in chikv!\n");
            exit(1);
          }
          nvAs[region] *= controlpar[2]; // This will be updated later
          nvAe[region] = pop_vecAe[region]->popsize;
          nvAi[region] = pop_vecAi[region]->popsize;
          //
          // printf("chikv after: %d %d %g %g %g   %g %d %d %d %d\n",region,TIME,controlpar[0],controlpar[1],controlpar[2],(envar + region*(4*(*finalT)+1))[TIME-1],nvAs[region],nvAe[region],nvAi[region],pop_vecAi[region]->popsize);
        }
      }
      // Quarantining symptomatic individuals
      if (controlpar[3]==1 && TIME>=controlpar[4]) {
        quarantine = 1;
      }
    }
    // Check for transmission
    /*
      WARNING!
      ========
      
      Information to decode environmental variables is hard coded here!
    */
    // time2here(3);
    for (region=0; region<(*numreg); region++)
      fecundity[region] = (envar + region*(4*(*finalT)+1) + 3*(*finalT))[TIME-1]; // +1 is for popdens
    // time2here(4);
    for (region=0; region<(*numreg); region++) {
      calculateAeCHIKVtransmission((*numreg),
                                   region,
                                   tprobs,
                                   fecundity,
                                   nS,
                                   nE,
                                   nI,
                                   nIc,
                                   nC,
                                   nR,
                                   nvAs,
                                   nvAe,
                                   nvAi,
                                   pop_humE,
                                   pop_vecAe,
                                   npV,
                                   npH,
                                   nbite);
    }
    // time2here(5);
    for (region=0; region<(*numreg); region++) {
      // Check if an experiment is scheduled
      if (param[introduce_time] > 0 &&
          param[introduce_infectious] > 0 &&
          !experiment[region] &&
          TIME >= param[introduce_time]) {
        calculateAeCHIKVexperiment(&pop_vecAi[region],
                                   &pop_humI[region],
                                   &nvAi[region],
                                   &nI[region],
                                   &nInew[region],
                                   &nS[region]);
        experiment[region] = 1;
      }
    }
    // time2here(6);
    for (region=0; region<(*numreg); region++) {
      // Recrd the output
      /*
        WARNING!
        ========
        
        Information to decode environmental variables is hard coded here!
      */
      output   = result + region*14*(*finalT);
      colT     = output + 0*(*finalT);
      colnS    = output + 1*(*finalT);
      colnE    = output + 2*(*finalT);
      colnI    = output + 3*(*finalT);
      colnIc   = output + 4*(*finalT);
      colnC    = output + 5*(*finalT);
      colnR    = output + 6*(*finalT);
      colnvAs  = output + 7*(*finalT);
      colnvAe  = output + 8*(*finalT);
      colnvAi  = output + 9*(*finalT);
      colnpV   = output + 10*(*finalT);
      colnpH   = output + 11*(*finalT);
      colnbite = output + 12*(*finalT);
      colnInew = output + 13*(*finalT);
      colT[TIME] = (double)TIME;
      colnS[TIME] = (double)nS[region];
      colnE[TIME] = (double)nE[region];
      colnI[TIME] = (double)nI[region];
      colnIc[TIME] = (double)nIc[region];
      colnC[TIME] = (double)nC[region];
      colnR[TIME] = (double)nR[region];
      colnvAs[TIME] = (double)nvAs[region];
      colnvAe[TIME] = (double)nvAe[region];
      colnvAi[TIME] = (double)nvAi[region];
      colnpV[TIME] = npV[region];
      colnpH[TIME] = npH[region];
      colnbite[TIME] = nbite[region];
      colnInew[TIME] = (double)nInew[region];
      if (CHECK(nS[region]) ||
          CHECK(nE[region]) ||
          CHECK(nI[region]) ||
          CHECK(nIc[region]) ||
          CHECK(nC[region]) ||
          CHECK(nR[region]) ||
          CHECK(nvAs[region]) ||
          CHECK(nvAe[region]) ||
          CHECK(nvAi[region]) ||
          CHECK(npV[region]) ||
          CHECK(npH[region]) ||
          (CHECK(nbite[region]) || CHECK(nbite[region])) ||
          CHECK(nInew[region])) {
        success[region] = 0;
        goto endall;
      }
      if (sim_mode==0 &&
          nInew[region]>0 &&
          param[introduce_time] > 0 &&
          param[introduce_infectious] > 0 &&
          TIME > param[introduce_time]) {
        success[region] = 3;
        goto endall;
      }
      /*
        printf("%d %d %d %d\n",
        region,
        total_human[region],
        nS[region]+nE[region]+nI[region]+nIc[region]+nC[region]+nR[region],
        nInew[region]);
      */
    }
  }
  //
 endall:
  //
  quarantine = 0;
  //
  for (region=0; region<(*numreg); region++) {
    // success[region] = (success[region]==0) ? 0 : 1;
    spop_destroy(&pop_humE[region]);
    spop_destroy(&pop_humI[region]);
    spop_destroy(&pop_humIc[region]);
    spop_destroy(&pop_vecAe[region]);
    spop_destroy(&pop_vecAi[region]);
  }
  free(tprobs);
  free(total_human);
  free(tH);
  free(fecundity);
  free(pop_humE);
  free(pop_humI);
  free(pop_humIc);
  free(pop_vecAe);
  free(pop_vecAi);
  free(experiment);
  free(nS);
  free(nE);
  free(nI);
  free(nIc);
  free(nC);
  free(nR);
  free(nvAs);
  free(nvAe);
  free(nvAi);
  free(npV);
  free(npH);
  free(nbite);
  free(nInew);
  //
  // gamma_dist_check();
  gamma_dist_destroy();
}

// ---------------------------------------------------------------------------

