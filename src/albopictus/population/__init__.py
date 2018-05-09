import numpy
from scipy.stats import gamma
from numpy.random import binomial

def gamma_dist_prob(array,mean,sd):
  theta = numpy.float64(sd) * sd / mean
  k = mean / theta
  return 1.0 - (gamma.sf(array+1,k,0,theta)/gamma.sf(array,k,0,theta))

class spop:
    def __init__(self,stochastic=True):
        self.stochastic = stochastic
        self.pop = numpy.ndarray((0,2),dtype=numpy.int32 if self.stochastic else numpy.float64)
    def add(self,number,age=0):
        i = self.pop[:,0] == age
        if any(i):
            self.pop[i,1] += number
        else:
            self.pop = numpy.append(self.pop,[[age,number]],axis=0)
    def iterate(self,dev=None,death=None,dev_mean=None,dev_sd=None,death_mean=None,death_sd=None):
        if dev==None:
            if dev_mean!=None and dev_sd!=None: # variable development probability
                dev = gamma_dist_prob(self.pop[:,0],dev_mean,dev_sd)
            else:
                print "Wrong development probability:",dev
                return
        if death==None:
            if death_mean!=None and death_sd!=None: # variable probability of death
                death = gamma_dist_prob(self.pop[:,0],death_mean,death_sd)
            else:
                print "Wrong probability of death:",death
                return
        #
        if self.stochastic:
            k = binomial(self.pop[:,1],death)
            self.pop[:,1] -= k
            # 
            d = binomial(self.pop[:,1],dev)
            self.pop[:,1] -= d
        else:
            self.pop[:,1] *= (1.0 - death)
            # 
            d = self.pop[:,1] * dev
            self.pop[:,1] *= 1.0 - dev
        # 
        self.pop[:,0] += 1
        self.pop = self.pop[self.pop[:,1]>0,:]
        return numpy.sum(d)

