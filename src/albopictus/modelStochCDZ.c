#ifdef R
#define modelHeader
#include "R.h"
#include "Rmath.h"
#else
#include "math.h"
#endif
// ---
#include <time.h>
#include <stdio.h>
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

#define NumParDis      32
#define NumMetDis      12
#define NumEnvVar      6

extern gsl_rng *RAND_GSL;

double gamma_rng(double mean, double sd) {
  double theta = sd * sd / mean;
  double k = mean / theta;
  return gsl_ran_gamma(RAND_GSL,k,theta);
}

// --------------------------------------------
// Simulation mode
// --------------------------------------------

char quarantine = 0; // quarantine symptomatic infections
char sim_mode = 1; // 0: check for secondary infection, 1: full simulation, 2: R0
void set_sim_mode(char *new_mode) {
  sim_mode = *new_mode;
}

// --------------------------------------------
// Gamma distribution
// --------------------------------------------

char gamma_mode_disease = MODE_GAMMA_HASH;
void set_gamma_mode_disease(char *new_mode) {
  gamma_mode_disease = *new_mode;
}

// --------------------------------------------
// Implementation of the dynamics
// --------------------------------------------

#define init_S               0
#define init_E               1
#define init_I               2
#define init_R               3
#define init_vAe             4
#define init_vAi             5
#define scale                6
#define gamma                7
#define gamma_sd             8
#define omega                9
#define omega_sd             10
#define susceptH             11
#define susceptV             12
#define vec_omega_1          13
#define vec_omega_2          14
#define vec_omega_3          15
#define vec_omega_4          16
#define vec_omega_5          17
#define vec_omega_6          18
#define vec_omega_sd         19
#define ovitraps             20
#define introduce_time       21
#define introduce_location   22
#define introduce_infectious 23
#define introduce_ivector    24
#define fixed_bites          25
#define report               26
#define qmove                27
#define connect              28
#define Khum                 29
#define Kvec                 30
#define Passim               31

double *param;
int *total_human;
double *tH;

void calc_tH(double *tprobs, int *numreg) {
  int region, regionB;
  for (region=0; region<(*numreg); region++) {
    tH[region] = 0;
    for (regionB=0; regionB<(*numreg); regionB++) {
      tH[region] += tprobs[regionB*(*numreg)+region]*total_human[regionB];
    }
  }
}

void calculateAeDiseaseInit(double daily_vector,
                          double      popdens,
                          spop          *humE,
                          spop          *humI,
                          spop         *vecAe,
                          spop         *vecAi,
                          int             *nS,
                          int             *nE,
                          int             *nI,
                          int             *nR,
                          int           *nvAs,
                          int           *nvAe,
                          int           *nvAi,
                          int          *nInew,
                          int          *nIimp) {
  (*humE) = spop_init();
  (*humI) = spop_init();
  (*vecAe) = spop_init();
  (*vecAi) = spop_init();

  (*nS) = round(popdens*param[init_S]*param[scale]);
  (*nE) = round(popdens*param[init_E]*param[scale]);
  (*nI) = round(popdens*param[init_I]*param[scale]);
  (*nR) = round(popdens*param[init_R]*param[scale]);
  (*nvAe) = round(param[init_vAe]*param[ovitraps]*param[scale]);
  (*nvAi) = round(param[init_vAi]*param[ovitraps]*param[scale]);
  (*nvAs) = max(0, round(daily_vector*param[ovitraps]*param[scale])-(*nvAe)-(*nvAi));
  (*nInew) = 0;
  (*nIimp) = 0;
  
  if ((*nE)>0) spop_add((*humE),0,(*nE));
  if ((*nI)>0) spop_add((*humI),0,(*nI));
  if ((*nvAe)>0) spop_add((*vecAe),0,(*nvAe));
  if ((*nvAi)>0) spop_add((*vecAi),0,(*nvAi));
}

void calculateAeDiseaseTransmission(int    numreg,
                                  int    regionA,
                                  double *tprobs,
                                  double *fecundity,
                                  int    *nS,
                                  int    *nE,
                                  int    *nI,
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
  double tHI = 0, tHS = 0;
  double *tHSc = (double *)malloc(numreg*sizeof(double));
  for (regionB=0; regionB<numreg; regionB++) {
    tmp = tprobs[regionB*numreg+regionA];
    tHI += tmp*nI[regionB]; // sum
    tHS += tmp*nS[regionB]; // sum
    tHSc[regionB] = tHS; // cumulative
  }
  //
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
}

void calculateAeDiseaseExperiment(spop   *vecAi,
                                spop   *humI,
                                int    *nvAi,
                                int    *nI,
                                int    *nInew,
                                int    *nIimp,
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
    (*nIimp) += infect;
    (*nS) -= infect;
  }
}

void calculateAeRandomImportation(double daily_importation_prob,
                                  spop   *humI,
                                  int    *nI,
                                  int    *nInew,
                                  int    *nIimp,
                                  int    *nS) {
  int infect = gsl_ran_poisson(RAND_GSL,daily_importation_prob);
  if (infect > (*nS))
    infect = (*nS);
  if (infect) {
    // Introduce an infectious human
    spop_add((*humI),0,infect);
    (*nI) += infect;
    (*nInew) += infect;
    (*nIimp) += infect;
    (*nS) -= infect;
  }
}

void calculateAeImportation(int     daily_importation,
                            spop   *humI,
                            int    *nI,
                            int    *nInew,
                            int    *nIimp,
                            int    *nS) {
  if (daily_importation > (*nS))
    daily_importation = (*nS);
  if (daily_importation) {
    // Introduce an infectious human
    spop_add((*humI),0,daily_importation);
    (*nI) += daily_importation;
    (*nInew) += daily_importation;
    (*nIimp) += daily_importation;
    (*nS) -= daily_importation;
  }
}

char applyAdultControl(spop   *vecAe,
                       spop   *vecAi,
                       int    *nvAs,
                       int    *nvAe,
                       int    *nvAi,
                       double frac) {
  int tmp;
  tmp = spop_survive((*vecAe),
                     0.0, 0.0, // no development
                     frac-1.0, -1.0, // daily probability of death
                     gamma_mode_disease);
  tmp += spop_survive((*vecAi),
                      0.0, 0.0, // no development
                      frac-1.0, -1.0, // daily probability of death
                      gamma_mode_disease);
  if (tmp) {
    printf("ERROR: Vector control failed in disease!\n");
    exit(1);
  }
  (*nvAs) *= frac; // This will be updated later
  (*nvAe) = (*vecAe)->popsize;
  (*nvAi) = (*vecAi)->popsize;
  //
  return 0;
}

char calculateAeDisease(double daily_mean_temperature,
                      double daily_vector,
                      double vector_survival_times,
                      double vector_survival_times_sd,
                      double daily_vector_fecundity,
                      spop   *humE,
                      spop   *humI,
                      spop   *vecAe,
                      spop   *vecAi,
                      int    *nS,
                      int    *nE,
                      int    *nI,
                      int    *nR,
                      int    *nvAs,
                      int    *nvAe,
                      int    *nvAi,
                      double *npV,
                      double *npH,
                      double *nbite,
                      int    *nInew,
                      int    *nIimp,
                      int    *totalH) {
  // Check the results of the development processes (and survival)
  // printf("B: %d, %d, %d, %d, %d\n",(*humE)->popsize,(*humI)->popsize,(*humIc)->popsize,(*vecAe)->popsize,(*vecAi)->popsize);
  int nhumE = spop_survive((*humE),
                           param[omega], param[omega_sd], // development
                           0.0, 0.0, // survival (immortal humans)
                           gamma_mode_disease);
  // Loose the asymptomatic cases
  int nhumIR = gsl_ran_binomial(RAND_GSL,param[Passim],nhumE);; // recover asymptomatic cases
  nhumE -= nhumIR;
  //
  int nhumI = spop_survive((*humI),
                           param[gamma], param[gamma_sd], // development
                           0.0, 0.0, // survival (immortal humans)
                           gamma_mode_disease);
  // Temperature-driven EIP
  double vec_omega = 0.0;
  if (param[vec_omega_4]==0.0 && param[vec_omega_5]==0.0 && param[vec_omega_5]==0.0) { // CHIKV - Dengue
    vec_omega = param[vec_omega_1] + exp(param[vec_omega_2] - param[vec_omega_3] * daily_mean_temperature);
  } else { // ZIKA
    vec_omega = param[vec_omega_1]+(param[vec_omega_4]-param[vec_omega_2]*(daily_mean_temperature-param[vec_omega_3]))/(param[vec_omega_5]+param[vec_omega_6]*(daily_mean_temperature-param[vec_omega_3]));
  }
  if (vec_omega < 0.0) vec_omega = 1e13;
  //
  int nvecAe = spop_survive((*vecAe),
                            vec_omega, param[vec_omega_sd], // development
                            vector_survival_times, vector_survival_times_sd, // survival
                            gamma_mode_disease);
  //
  int nvecAi = spop_survive((*vecAi),
                            0.0, 0.0, // development (persistent infectious stage)
                            vector_survival_times, vector_survival_times_sd, // survival
                            gamma_mode_disease);
  if (nvecAi != 0) {
    printf("ERROR: Development observed for infectious vector!\n");
    exit(1);
  }
  // printf("A: %d, %d, %d, %d, %d\n",nhumE,nhumI,nhumIc,nvecAe,nvecAi);

  // Decide the fate of acute infections
  nhumIR += nhumI; // recover completely

  if (sim_mode == 2) { // Calculate R0
    nhumIR += nhumE;
    nhumE = 0;
  }

  // Perform state transformations
  if (!quarantine) // not quarantine
    spop_add((*humI),0,nhumE);
  spop_add((*vecAi),0,nvecAe); // newly infected vector population
  // No imported cases as a result of disease progression
  (*nIimp) = 0;
  // Number of newly infectious symptomatic cases which are reported (endemic)
  (*nInew) = gsl_ran_binomial(RAND_GSL,param[report],nhumE);
  if (quarantine) { // quarantine
    spop_add((*humI),0,nhumE-(*nInew));
    (*totalH) -= (*nInew);
  }
  (*nE) = (*humE)->popsize; // exposed human population
  (*nI) = (*humI)->popsize; // infectious human population
  (*nR) += nhumIR; // recovered human population
  (*nvAi) = (*vecAi)->popsize; // infected vector population
  (*nvAe) = (*vecAe)->popsize; // exposed vector population
  // daily_vector from modelDelayAalbopictus is for the end of the day
  (*nvAs) = max(0, round(daily_vector*param[ovitraps]*param[scale])-(*nvAe)-(*nvAi));
  //
  // printf("B: %d %d %d\n",*nvAs,*nvAe,*nvAi);
  return 0;
}

// --------------------------------------------
// sim_model
// --------------------------------------------

void numparModel(int *np, int *nm) {
  *np = NumParDis;
  *nm = NumMetDis;
}

void param_model(char **names, double *par) {
  char temp[NumMetDis+NumParDis][256] = {
    "colnS","colnE","colnI","colnR","colnvAs","colnvAe","colnvAi","colnpV","colnpH","colnbite","colnInew","colnIimp",
    "init_S","init_E","init_I","init_R","init_vAe","init_vAi","scale","gamma","gamma_sd","omega","omega_sd","susceptH","susceptV","vec_omega_1","vec_omega_2","vec_omega_3","vec_omega_4","vec_omega_5","vec_omega_6","vec_omega_sd","ovitraps","introduce_time","introduce_location","introduce_infectious","introduce_ivector","fixed_bites","report","qmove","connect","Khum","Kvec","Passim"
  };
  int i;
  for (i=0; i<(NumMetDis+NumParDis); i++)
    names[i] = strdup(temp[i]);
  // ---
  // 
  par[init_S] = 1.0; // susceptible (% popdens)
  par[init_E] = 0; // latent (% popdens)
  par[init_I] = 0; // infectious (a/symptomatic) (% popdens)
  par[init_R] = 0; // recovered (% popdens)
  par[init_vAe] = 0; // adult vector (latent)
  par[init_vAi] = 0; // adult vector (infectious)
  par[scale] = 1.0; // simulation span (affects double -> int conversion)
  par[gamma] = 4.5; // infectious period in humans (2-7 days)
  par[gamma_sd] = 1.0; //
  par[omega] = 3.5; // latent period in humans (2-4 days)
  par[omega_sd] = 0.5; //
  par[susceptH] = 0.65; //0.8; // human susceptibility to infections (50-80%)
  par[susceptV] = 0.85; //1.0; // vector susceptibility to infections (70-100%)
  //
  /* The EIP is modified to accommodate regulation by temperature (from AedesRisk)
   * Model incorporating extrinsic incubation period = 4 + exp(5.15 - 0.123*T)
   * https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0089783
   */
  par[vec_omega_1] = 4.0; 
  par[vec_omega_2] = 5.15;
  par[vec_omega_3] = 0.123; 
  par[vec_omega_4] = 0.0; 
  par[vec_omega_5] = 0.0;
  par[vec_omega_6] = 0.0; 
  par[vec_omega_sd] = 0.1; // Std. dev. of the mean
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
  //
  par[Passim] = 0.0; // no asymptomatic cases
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
  // rng_setup();
  // rng_setup_seed(123);
  
  int r;
  param = par;
  double *controlpar = 0;
  if ((*control)) {
    controlpar = param + NumParDis;
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
  spop *pop_vecAe = (spop *)malloc((*numreg)*sizeof(spop));
  spop *pop_vecAi = (spop *)malloc((*numreg)*sizeof(spop));

  char case_aut = 0;
  char case_imp = 0;
  double *daily_vector_control = (double *)calloc((*finalT),sizeof(double));

  char *experiment = (char *)calloc((*numreg),sizeof(char));

  int *nS = (int *)calloc((*numreg),sizeof(int));
  int *nE = (int *)calloc((*numreg),sizeof(int));
  int *nI = (int *)calloc((*numreg),sizeof(int));
  int *nR = (int *)calloc((*numreg),sizeof(int));
  int *nvAs = (int *)calloc((*numreg),sizeof(int));
  int *nvAe = (int *)calloc((*numreg),sizeof(int));
  int *nvAi = (int *)calloc((*numreg),sizeof(int));
  double *npV = (double *)calloc((*numreg),sizeof(double));
  double *npH = (double *)calloc((*numreg),sizeof(double));
  double *nbite = (double *)calloc((*numreg),sizeof(double));
  int *nInew = (int *)calloc((*numreg),sizeof(int));
  int *nIimp = (int *)calloc((*numreg),sizeof(int));

  double *environ;
  double *output;
  double *daily_mean_temperature;
  double *daily_vector;
  double *vector_survival_times;
  double *vector_survival_times_sd;
  double *daily_vector_fecundity;
  double *daily_importation;
  double *colT;
  double *colnS;
  double *colnE;
  double *colnI;
  double *colnR;
  double *colnvAs;
  double *colnvAe;
  double *colnvAi;
  double *colnpV;
  double *colnpH;
  double *colnbite;
  double *colnInew;
  double *colnIimp;

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
  for (region=0; region<(*numreg); region++) {
    calculateAeDiseaseInit(*(envar + region*(NumEnvVar*(*finalT)+1) + 1*(*finalT)), // daily_vector (+1 is for popdens)
                           *(envar + region*(NumEnvVar*(*finalT)+1) + NumEnvVar*(*finalT)), // popdens (+1 is for popdens)
                           &pop_humE[region],
                           &pop_humI[region],
                           &pop_vecAe[region],
                           &pop_vecAi[region],
                           &nS[region],
                           &nE[region],
                           &nI[region],
                           &nR[region],
                           &nvAs[region],
                           &nvAe[region],
                           &nvAi[region],
                           &nInew[region],
                           &nIimp[region]);
    total_human[region] = nS[region]+nE[region]+nI[region]+nR[region];
    //
    output   = result + region*(NumMetDis+1)*(*finalT);
    colT     = output + 0*(*finalT);
    colnS    = output + 1*(*finalT);
    colnE    = output + 2*(*finalT);
    colnI    = output + 3*(*finalT);
    colnR    = output + 4*(*finalT);
    colnvAs  = output + 5*(*finalT);
    colnvAe  = output + 6*(*finalT);
    colnvAi  = output + 7*(*finalT);
    colnpV   = output + 8*(*finalT);
    colnpH   = output + 9*(*finalT);
    colnbite = output + 10*(*finalT);
    colnInew = output + 11*(*finalT);
    colnIimp = output + 12*(*finalT);
    colT[TIME] = (double)TIME;
    colnS[TIME] = (double)nS[region];
    colnE[TIME] = (double)nE[region];
    colnI[TIME] = (double)nI[region];
    colnR[TIME] = (double)nR[region];
    colnvAs[TIME] = (double)nvAs[region];
    colnvAe[TIME] = (double)nvAe[region];
    colnvAi[TIME] = (double)nvAi[region];
    colnpV[TIME] = npV[region];
    colnpH[TIME] = npH[region];
    colnbite[TIME] = nbite[region];
    colnInew[TIME] = (double)nInew[region];
    colnIimp[TIME] = (double)nIimp[region];
  }
  //
  calc_tH(tprobs,numreg);
  //
  for (TIME=1; TIME<(*finalT); TIME++) {
    // printf("TIME: %d\n",TIME);
    for (region=0; region<(*numreg); region++) {
      // Take the first step
      /*
        WARNING!
        ========

        Information to decode environmental variables is hard coded here!
       */
      environ                   = envar + region*(NumEnvVar*(*finalT)+1); // +1 is for popdens
      daily_mean_temperature    = environ + 0*(*finalT);
      daily_vector              = environ + 1*(*finalT);
      vector_survival_times     = environ + 2*(*finalT);
      vector_survival_times_sd  = environ + 3*(*finalT);
      daily_vector_fecundity    = environ + 4*(*finalT);
      daily_importation         = environ + 5*(*finalT);
      
      if (calculateAeDisease(daily_mean_temperature[TIME-1],
                       daily_vector[TIME-1],
                       vector_survival_times[TIME-1],
                       vector_survival_times_sd[TIME-1],
                       daily_vector_fecundity[TIME-1],
                       &pop_humE[region],
                       &pop_humI[region],
                       &pop_vecAe[region],
                       &pop_vecAi[region],
                       &nS[region],
                       &nE[region],
                       &nI[region],
                       &nR[region],
                       &nvAs[region],
                       &nvAe[region],
                       &nvAi[region],
                       &npV[region],
                       &npH[region],
                       &nbite[region],
                       &nInew[region],
                       &nIimp[region],
                       &total_human[region])) {
        success[region] = 0;
        goto endall;
      }
    }
    // Apply control measures
    if ((*control)) {
      // Adult reduction
      if (TIME>=controlpar[0] && TIME<controlpar[1]) {
        //
        for (region=0; region<(*numreg); region++) { // Apply to all regions at once!
          // printf("disease before: %d   %d %d %g %g %g   %g %d %d %d %d\n",TIME,region,TIME,controlpar[0],controlpar[1],controlpar[2],(envar + region*(NumEnvVar*(*finalT)+1))[TIME-1],nvAs[region],nvAe[region],nvAi[region],pop_vecAi[region]->popsize);
          if (applyAdultControl(&pop_vecAe[region],
                                &pop_vecAi[region],
                                &nvAs[region],
                                &nvAe[region],
                                &nvAi[region],
                                controlpar[2])) {
            success[region] = 0;
            goto endall;
          }
          // printf("disease after: %d   %d %d %g %g %g   %g %d %d %d %d\n",TIME,region,TIME,controlpar[0],controlpar[1],controlpar[2],(envar + region*(NumEnvVar*(*finalT)+1))[TIME-1],nvAs[region],nvAe[region],nvAi[region],pop_vecAi[region]->popsize);
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
    for (region=0; region<(*numreg); region++)
      fecundity[region] = (envar + region*(NumEnvVar*(*finalT)+1) + 4*(*finalT))[TIME-1]; // daily_vector_fecundity (+1 is for popdens)
    for (region=0; region<(*numreg); region++) {
      calculateAeDiseaseTransmission((*numreg),
                                   region,
                                   tprobs,
                                   fecundity,
                                   nS,
                                   nE,
                                   nI,
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
    for (region=0; region<(*numreg); region++) {
      // Check if an experiment is scheduled
      if (param[introduce_time] > 0 &&
          param[introduce_infectious] > 0 &&
          !experiment[region] &&
          TIME >= param[introduce_time]) {
        calculateAeDiseaseExperiment(&pop_vecAi[region],
                                   &pop_humI[region],
                                   &nvAi[region],
                                   &nI[region],
                                   &nInew[region],
                                   &nIimp[region],
                                   &nS[region]);
        experiment[region] = 1;
      }
      // Check if random importation is configured for the region
      daily_importation = envar + region*(NumEnvVar*(*finalT)+1) + 5*(*finalT); // daily_importation (+1 is for popdens)
      if (daily_importation[TIME] > 0) {
        calculateAeImportation((int)round(daily_importation[TIME]),
                               &pop_humI[region],
                               &nI[region],
                               &nInew[region],
                               &nIimp[region],
                               &nS[region]);
      }
    }
    //
    // Apply reactive control
    if ((*control)) {
      // Reactive control
      if (controlpar[5]>0 && controlpar[6]>0 && controlpar[7]>0 && controlpar[8]>0 && controlpar[9]>0 && controlpar[10]>=0) {
        // Check if new control should be scheduled
        for (case_aut=0, case_imp=0, region=0; region<(*numreg); region++) {
          // Assume all imported cases are detected and they trigger control
          case_imp = (nIimp[region] > 0);
          // Reported/detected autochthonous infectious cases could trigger control
          case_aut = (nInew[region] > nIimp[region]);
          //
          if ((controlpar[5] >= 2.5 && controlpar[5] < 3.5 && (case_imp || case_aut)) ||   // 3: After both imported and autochthonous cases
              (controlpar[5] >= 1.5 && controlpar[5] < 2.5 && case_aut) ||                 // 2: Only after autochthonous cases
              (controlpar[5] >= 0.5 && controlpar[5] < 1.5 && case_imp)) {                 // 1: Only after imported cases
              //
              // printf("%g %d %d %d %d %d\n",controlpar[5],TIME,case_aut,case_imp,nInew[region],nIimp[region]);
              //
              // Schedule control (in all regions)
              int tfn, tfm;
              int lnt = (int)round(controlpar[9]);
              int mnt = (int)round(controlpar[10]);
              int rnum = (int)round(gamma_rng(controlpar[6],controlpar[7]));
              int tf = min((*finalT)-1, TIME + rnum);
              char reg = 1;
              //
              tfn = max(0, tf - mnt);
              tfm = min((*finalT)-1, tf+lnt+mnt);
              for (; 
                   tfn < tfm; 
                   tfn++) {
                    if (daily_vector_control[tfn] > 0.0) { // Too close, do not register!
                      reg = 0;
                      break;
                    }
              }
              //
              if (reg) {
                tfn = tf;
                tfm = min((*finalT)-1, tf + lnt);
                for (; 
                     tfn < tfm;
                     tfn++) {
                  daily_vector_control[tfn] = controlpar[8];
                }
              }
          }
          // If any control is scheduled for today, apply (in this region)
          if (daily_vector_control[TIME] > 0) {
            // printf("disease before: %d %d %g %g %g   %g %d %d %d %d\n",region,TIME,controlpar[6],controlpar[7],controlpar[8],(envar + region*(NumEnvVar*(*finalT)+1))[TIME-1],nvAs[region],nvAe[region],nvAi[region],pop_vecAi[region]->popsize);
            if (applyAdultControl(&pop_vecAe[region],
                                  &pop_vecAi[region],
                                  &nvAs[region],
                                  &nvAe[region],
                                  &nvAi[region],
                                  daily_vector_control[TIME])) {
              success[region] = 0;
              goto endall;
            }
            // printf("disease after: %d %d %g %g %g   %g %d %d %d %d\n",region,TIME,controlpar[6],controlpar[7],controlpar[8],(envar + region*(NumEnvVar*(*finalT)+1))[TIME-1],nvAs[region],nvAe[region],nvAi[region],pop_vecAi[region]->popsize);
          }
        }
      }
    }
    //
    for (region=0; region<(*numreg); region++) {
      // Recrd the output
      /*
        WARNING!
        ========
        
        Information to decode environmental variables is hard coded here!
      */
      output   = result + region*(NumMetDis+1)*(*finalT);
      colT     = output + 0*(*finalT);
      colnS    = output + 1*(*finalT);
      colnE    = output + 2*(*finalT);
      colnI    = output + 3*(*finalT);
      colnR    = output + 4*(*finalT);
      colnvAs  = output + 5*(*finalT);
      colnvAe  = output + 6*(*finalT);
      colnvAi  = output + 7*(*finalT);
      colnpV   = output + 8*(*finalT);
      colnpH   = output + 9*(*finalT);
      colnbite = output + 10*(*finalT);
      colnInew = output + 11*(*finalT);
      colnIimp = output + 12*(*finalT);
      colT[TIME] = (double)TIME;
      colnS[TIME] = (double)nS[region];
      colnE[TIME] = (double)nE[region];
      colnI[TIME] = (double)nI[region];
      colnR[TIME] = (double)nR[region];
      colnvAs[TIME] = (double)nvAs[region];
      colnvAe[TIME] = (double)nvAe[region];
      colnvAi[TIME] = (double)nvAi[region];
      colnpV[TIME] = npV[region];
      colnpH[TIME] = npH[region];
      colnbite[TIME] = nbite[region];
      colnInew[TIME] = (double)nInew[region];
      colnIimp[TIME] = (double)nIimp[region];
      if (CHECK(nS[region]) ||
          CHECK(nE[region]) ||
          CHECK(nI[region]) ||
          CHECK(nR[region]) ||
          CHECK(nvAs[region]) ||
          CHECK(nvAe[region]) ||
          CHECK(nvAi[region]) ||
          CHECK(npV[region]) ||
          CHECK(npH[region]) ||
          (CHECK(nbite[region]) || CHECK(nbite[region])) ||
          CHECK(nInew[region]) ||
          CHECK(nIimp[region])) {
        success[region] = 0;
        goto endall;
      }
      if (sim_mode==0 &&
          nInew[region]>0 &&
          nInew[region]>nIimp[region] &&
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
    spop_destroy(&pop_vecAe[region]);
    spop_destroy(&pop_vecAi[region]);
  }
  free(tprobs);
  free(total_human);
  free(tH);
  free(fecundity);
  free(pop_humE);
  free(pop_humI);
  free(pop_vecAe);
  free(pop_vecAi);
  free(experiment);
  free(daily_vector_control);
  free(nS);
  free(nE);
  free(nI);
  free(nR);
  free(nvAs);
  free(nvAe);
  free(nvAi);
  free(npV);
  free(npH);
  free(nbite);
  free(nInew);
  free(nIimp);
  //
  // gamma_dist_check();
  gamma_dist_destroy();
  // rng_destroy();
}

// ---------------------------------------------------------------------------

