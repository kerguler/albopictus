#ifndef SPOP_H
#define SPOP_H

#include <float.h>

#define DPOP_EPS      DBL_MIN
#define DPOP_MAX_DAYS 1000000

typedef union {
  unsigned int i;
  double d;
} sdnum;

typedef struct individual_st {
  unsigned int age;
  unsigned int devcycle;
  unsigned int development;
  sdnum number;
} individual_data;

#define is_greater(s,p) (((s).age > (p).age) || ((s).age == (p).age && (s).devcycle > (p).devcycle) || ((s).age == (p).age && (s).devcycle == (p).devcycle && (s).development > (p).development))
#define is_smaller(s,p) (((s).age < (p).age) || ((s).age == (p).age && (s).devcycle < (p).devcycle) || ((s).age == (p).age && (s).devcycle == (p).devcycle && (s).development < (p).development))
#define is_equal(s,p) ((s).age == (p).age && (s).devcycle == (p).devcycle && (s).development == (p).development)
#define is_empty(sp,s) (((sp)->stochastic && (s).number.i==0) || (!((sp)->stochastic) && (s).number.d<=DPOP_EPS))

#define get_lchild(num) (((num)<<1)+1)
#define get_rchild(num) (((num)<<1)+2)
#define get_parent(num) ((num) ? ((num)-1)>>1 : 0)

typedef struct population_st {
  individual_data *individuals;
  unsigned int ncat;
  unsigned int cat;
  sdnum size;
  sdnum dead;
  sdnum developed;
  void *devtable;
  unsigned char gamma_mode;
  unsigned char stochastic;
} *spop;

spop spop_init(unsigned char, unsigned char);
void spop_empty(spop);
void spop_destroy(spop*);
void spop_print(spop);

void swap(spop, individual_data *, individual_data *);

#define spop_add(s,age,devcycle,development,number) {   \
    sdnum tmp;                                          \
    if ((s)->stochastic)                                \
      tmp.i = (int)(number);                            \
    else                                                \
      tmp.d = (double)(number);                            \
    spop_sdadd((s),(age),(devcycle),(development),tmp);    \
  }
void spop_sdadd(spop, unsigned int, unsigned int, unsigned int, sdnum);

void spop_popadd(spop, spop);

typedef double (*prob_func)(unsigned int, double, double, double);
typedef void (*iter_func)(const individual_data*, double*, double*, double*);
void spop_iterate(spop, double, double, double, iter_func, double, double, double, iter_func, unsigned char);

#endif
