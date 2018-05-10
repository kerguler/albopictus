#ifndef SDPOP_H
#define SDPOP_H

#define DPOP_EPS 1e-14

typedef double (*prob_func)(unsigned int, double, double, double);

typedef union {
  unsigned int i;
  double d;
} sdnum;

typedef struct individual_st
{
  unsigned int age;
  sdnum number;
  struct individual_st *next;
} *individual_data;

typedef struct population_st
{
  individual_data individual;
  sdnum size;
  sdnum dead;
  sdnum developed;
  unsigned char gamma_mode;
  unsigned char stochastic;
} *sdpop;

#define sdpop_add(s,age,number) {               \
    if ((number)>DPOP_EPS) {                    \
      sdnum tmp;                                \
      if ((s)->stochastic)                      \
        tmp.i = (int)(number);                  \
      else                                      \
        tmp.d = (double)(number);               \
      sdpop_sdadd((s),(age),tmp);               \
    }                                           \
  }

sdpop sdpop_init(unsigned char, unsigned char);
void sdpop_sdadd(sdpop, unsigned int, sdnum);
void sdpop_print(sdpop);
void sdpop_remove_next(sdpop, individual_data);
void sdpop_empty(sdpop);
void sdpop_destroy(sdpop*);
void sdpop_kill(sdpop, double);

prob_func assign_prob_func(sdpop, double, double, double);

void sdpop_iterate(sdpop, double, double, double, double, double, double);

#endif
