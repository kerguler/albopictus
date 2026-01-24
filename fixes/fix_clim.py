# This script has been applied to adapt "clim" to vector08c

import numpy
import json

from matplotlib import pylab as plt

import albopictus as aa

vector = aa.vector08c
vecpar = aa.param['Q4.a100+1']

climfile = "src/albopictus/data/climate.json"

clim = json.load(open(climfile,"r"))

for prv in clim:
    for i in range(len(clim[prv])):
        l = len(clim[prv][i]['dates'])
        clim[prv][i]['envar'] += [clim[prv][i]['envar'][-1] for j in range(l-1)]
        for j in numpy.arange(l*2,l*3):
            clim[prv][i]['envar'][j] /= 1000.0

json.dump(clim,open(climfile,'w'))

sim = vector.simPar(clim['PR'][0],vecpar[0])

plt.plot(sim['colegg'][1:])
plt.show()
