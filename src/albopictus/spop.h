typedef struct people_st
{
  int development;
  int batchsize;
  struct people_st *next;
} *people_data;

typedef struct population_st
{
  people_data people;
  int popsize;
} *spop;

spop spop_init(void);
void spop_add(spop, int, int);
void spop_remove_next(spop, people_data);
void spop_empty(spop);
void spop_destroy(spop *);
void spop_print(spop);
int spop_survive(spop, double, double, double, double, char);
void spop_swap(spop, spop);
