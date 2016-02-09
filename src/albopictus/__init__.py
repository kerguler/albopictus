"""

Environmentally-driven population dynamics model of Aedes albopictus

This package implements the stage- and age-structured population dynamics model of temperate Aedes albopictus (Skuse) as described in Erguler et.al.. The package includes a representative environmental dataset for the Emilia-Romagna region. The package also includes three sets of parameter vectors sampled from different posterior modes offering explanation to the data obtained from surveillance.

References

     Erguler K, Smith-Unna SE, Waldock J, Proestos Y, Christophides GK, Lelieveld J, Parham PE. Large-scale modelling of the environmentally-driven population dynamics of temperate Aedes albopictus (Skuse). PLOS ONE (in press)

See also

     hoppMCMC - an adaptive basin-hopping Markov-chain Monte Carlo algorithm for Bayesian optimisation

vector.simPar

    Simulate the vector model at a specific environmental grid point with a given parameter vector

    Parameters
    ----------

         clim:
              environmental dataset to be used
              it must be a dictionary with keys 'dates' and 'envar'
              clim['dates'] holds and array of dates corresponding to each data point with a length of final_time
              clim['envar'] should hold an array with the following data in the given order and size
                   photoperiod [final_time]
                   mean air temperature [final_time]
                   daily precipitation [final_time]
                   human population density [1]
              the object 'clim' holds an example dataset for the Emilia-Romagna region

         pr:
              an array of size vector.numpar holding the parameter vector
     
Examples

     import albopictus as aa

     # The following line simulates the model at the first environmental grid point
     # of Bologna using the first parameter vector of the posterior mode Q1.
     sim = aa.vector.simPar(aa.clim['BO'][0],aa.param['Q1'][0])

     import pylab

     pylab.ylim([0,50])
     pylab.plot(aa.clim['BO'][0]['dates'],sim['colegg'])
     pylab.show()

"""

# modelAalbopictus - climateData ------------------------- //
from distutils.sysconfig import get_python_lib
import pkg_resources
import numpy
import json
from datetime import date

param = json.load(open(pkg_resources.resource_filename(__name__, "data/posterior.json"),"r"))
prior = json.load(open(pkg_resources.resource_filename(__name__, "data/prior.json"),"r"))
clim = json.load(open(pkg_resources.resource_filename(__name__, "data/climate.json"),"r"))
for pr in clim:
    for clm in clim[pr]:
        clm['dates'] = [date.fromordinal(d) for d in clm['dates']]

provinces = ["Bologna","Ferrara","Modena","Piacenza","Parma","Ravenna","Reggio Emilia"]
prvn = ["BO","FE","MO","PC","PR","RA","RE"]

# modelAalbopictus --------------------------------------- //

from readModel import prepareModel

vector = prepareModel(pkg_resources.resource_filename(__name__, "modelAalbopictus.so"))

