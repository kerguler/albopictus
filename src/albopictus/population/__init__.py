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
        self.pop = numpy.ndarray((0,4),dtype=numpy.int32 if self.stochastic else numpy.float64)
        self.developed = 0
        self.dead = 0
        self.size = 0
        self.devtable = numpy.ndarray((0,4),dtype=numpy.int32 if self.stochastic else numpy.float64)
    #
    def add(self,pop):
        pop = numpy.array(pop)
        if (not isinstance(pop, numpy.ndarray) or pop.shape[1] != 4):
            print "ERROR: Wrong array dimension (4 columns needed for age, devcycle, development, and age)"
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
    def iterate(self,dev=None,death=None,dev_mean=None,dev_sd=None,death_mean=None,death_sd=None):
        if dev==None:
            if dev_mean!=None and dev_sd!=None: # variable development probability
                if dev_sd == 0:
                    dev = numpy.float64(self.pop[:,2] >= dev_mean - 1.0)
                else:
                    dev = gamma_dist_prob(self.pop[:,2],dev_mean,dev_sd)
            else:
                print "Wrong development probability:",dev
                return
        if death==None:
            if death_mean!=None and death_sd!=None: # variable probability of death
                if death_sd == 0:
                    death = numpy.float64(self.pop[:,0] >= death_mean - 1.0)
                else:
                    death = gamma_dist_prob(self.pop[:,0],death_mean,death_sd)
            else:
                print "Wrong probability of death:",death
                return
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
        self.pop[:,0] += 1
        self.pop[:,2] += 1
        #
        self.devtable = self.pop.copy()
        self.devtable[:,3] = d
        self.devtable[:,1] += 1
        self.devtable[:,2] = 0
        self.devtable = self.devtable[self.devtable[:,3]>0,:]
        # 
        self.pop = self.pop[self.pop[:,3]>0,:]
        self.dead = numpy.sum(k)
        self.developed = numpy.sum(d)
        self.size = numpy.sum(self.pop[:,3])

