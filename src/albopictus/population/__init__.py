"""

This is a Python implementation of the age-structured population dynamics models described in Erguler et al. 2016 and 2017. The spop class records the number and age of individuals and implements two processes to exit from the population: development and death. The two processes act upon the population sequentially; survival is imposed prior to development. If the population survives for one day, then, it is allowed to grow and complete its development. Survival and development are defined either with an age-independent daily probability, or an age-dependent gamma- or negative binomial-distributed probability.

spop

     Parameters
     ----------

     stochastic: a logical value indicating a deterministic or a stochastic population dynamics

     prob: a character string indicating the basis of age-dependent survival or development (gamma: gamma-distributed (default), nbinom: negative binomial-distributed)

Examples

     from albopictus import population

     # Generate a population with stochastic dynamics
     s = population.spop(stochastic=True)

     # Add 1000 20-day-old individuals ([[age, devcycle, development, number]])
     s.add([[20,0,20,1000]])

     # Iterate one day without death and assume development in 20 (+-5) days (gamma-distributed)
     s.iterate(dev_mean=20,dev_sd=5,death=0)
     print(s.developed)

     # Iterate another day assuming no development but age-dependent survival
     # Let each individual survive for 20 days (+-5) (gamma-distributed)
     s.iterate(death_mean=20,death_sd=5,dev=0)
     print(s.dead)
     # Note that the previous values of developed and dead will be overwritten by this command

     # Generate a deterministic population and observe the difference
     s = population.spop(stochastic=False)
     s.add([[20,0,20,1000]])

     s.iterate(dev_mean=20,dev_sd=5,death=0)
     print(s.developed)

     s.iterate(death_mean=20,death_sd=5,dev=0)
     print(s.dead)

"""

import numpy
from scipy.stats import gamma, nbinom
from numpy.random import binomial, choice

prob_list = ['gamma','nbinom']

hash = { p:{} for p in prob_list }

def nbinom_dist_prob(array,mean,sd):
  p = numpy.float64(mean) / (sd * sd)
  r = mean * p / (1.0 - p)
  return 1.0 - (nbinom.sf(array+1,r,p))/(nbinom.sf(array,r,p))

def nbinom_dist_hash(array,mean,sd):
  ret = []
  for a in array:
      if not (mean in hash['nbinom']):
          hash['nbinom'][mean] = {}
      if not (sd in hash['nbinom'][mean]):
          hash['nbinom'][mean][sd] = nbinom_dist_prob(a,mean,sd)
      ret.append(hash['nbinom'][mean][sd])
  return numpy.array(ret)

def gamma_dist_prob(array,mean,sd):
  theta = numpy.float64(sd) * sd / mean
  k = mean / theta
  return 1.0 - (gamma.sf(array+1,k,0,theta)/gamma.sf(array,k,0,theta))

def gamma_dist_hash(array,mean,sd):
  ret = []
  for a in array:
      if not (mean in hash['gamma']):
          hash['gamma'][mean] = {}
      if not (sd in hash['gamma'][mean]):
          hash['gamma'][mean][sd] = gamma_dist_prob(a,mean,sd)
      ret.append(hash['gamma'][mean][sd])
  return numpy.array(ret)

class spop:
    def __init__(self,stochastic=True,prob='gamma'):
        if not (prob in prob_list):
            print("ERROR: prob can be one of the following:")
            for p in prob_list:
                print("-> %s" %(p))
            return None
        self.prob = prob
        self.prob_fun = None
        if self.prob == 'gamma':
            self.prob_fun = gamma_dist_prob
        elif self.prob == 'nbinom':
            self.prob_fun = nbinom_dist_prob
        self.stochastic = stochastic
        self.pop = numpy.ndarray((0,4),dtype=numpy.int32 if self.stochastic else numpy.float64)
        self.developed = 0
        self.dead = 0
        self.size = 0
        self.devtable = numpy.ndarray((0,4),dtype=numpy.int32 if self.stochastic else numpy.float64)
    #
    def add(self,pop):
        if self.stochastic:
            pop = numpy.array(pop,dtype=numpy.int32)
        else:
            pop = numpy.array(pop,dtype=numpy.float64)
        if (not isinstance(pop, numpy.ndarray) or pop.shape[1] != 4):
            print("ERROR: Wrong array dimension (4 columns needed for age, devcycle, development, and age)")
            return
        if (pop.shape[0] == 0):
            return
        for row in pop:
            i = (self.pop[:,0] == row[0]) & (self.pop[:,1] == row[1]) & (self.pop[:,2] == row[2])
            if any(i):
                self.pop[i,3] += row[3]
            else:
                self.pop = numpy.append(self.pop,[row],axis=0)
        self.size = numpy.sum(self.pop[:,3])
    #
    def checkNaN(self,d):
        if isinstance(d,numpy.float64) and numpy.isnan(d):
            return 1.0
        if isinstance(d,numpy.ndarray):
            d[numpy.isnan(d)] = 1.0
        return d
    #
    def extract(self,num):
        if num <= 0 or self.size <= 0:
            return []
        if self.stochastic:
            ret = self.pop.copy()
            ret[:,3] = 0
            for j in numpy.arange(num):
                i = choice(numpy.arange(self.pop.shape[0]), p=self.pop[:,3] / self.size)
                self.pop[i,3] -= 1
                self.size -= 1
                ret[i,3] += 1
            ret = ret[ret[:,3]>0,:]
        else:
            ret = self.pop.copy()
            ret[:,3] *= num / self.size
            self.pop[:,3] -= ret[:,3]
            self.size -= num
        return ret
    #
    def iterate(self,dev=None,death=None,dev_mean=None,dev_sd=None,death_mean=None,death_sd=None,pause=False):
        if dev is None:
            if (not (dev_mean is None)) and ((dev_sd is None) or numpy.all(dev_sd == 0)):
                dev = numpy.float64(self.pop[:,2] >= numpy.array(dev_mean) - 1.0)
            elif (not (dev_mean is None)) and (not ((dev_sd is None) or numpy.all(dev_sd == 0))):
                dev = self.prob_fun(self.pop[:,2],dev_mean,dev_sd)
            else:
                print("Wrong probability: %s" %(dev))
                return numpy.nan
        dev = self.checkNaN(dev)
        #
        if death is None:
            if (not (death_mean is None)) and ((death_sd is None) or numpy.all(death_sd == 0)):
                death = numpy.float64(self.pop[:,0] >= numpy.array(death_mean) - 1.0)
            elif (not (death_mean is None)) and (not ((death_sd is None) or numpy.all(death_sd == 0))):
                death = self.prob_fun(self.pop[:,0],death_mean,death_sd)
            else:
                print("Wrong probability: %s" %(death))
                return numpy.nan
        death = self.checkNaN(death)
        #
        if self.stochastic:
            k = binomial(self.pop[:,3],death)
            self.pop[:,3] -= k
            # 
            d = binomial(self.pop[:,3],dev)
            self.pop[:,3] -= d
        else:
            k = self.pop[:,3] * death
            self.pop[:,3] -= k
            # 
            d = self.pop[:,3] * dev
            self.pop[:,3] -= d
        #
        if not pause:
            self.pop[:,0] += 1
            self.pop[:,2] += 1
        #
        self.devtable = self.pop.copy()
        if not pause:
            self.devtable[:,1] += 1
            self.devtable[:,2] = 0
        self.devtable[:,3] = d
        self.devtable = self.devtable[self.devtable[:,3]>0,:]
        # 
        self.pop = self.pop[self.pop[:,3]>0,:]
        self.dead = numpy.sum(k)
        self.developed = numpy.sum(d)
        self.size = numpy.sum(self.pop[:,3])
    #
    def flush(self):
        self.developed = 0
        self.dead = 0
        self.devtable = numpy.ndarray((0,4),dtype=numpy.int32 if self.stochastic else numpy.float64)
