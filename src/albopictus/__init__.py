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
     
         cpr: (optional)
              an array of vector control parameters of the following order:
              [day0, dayf, daily fraction of breeding site reduction,
               day0, dayf, daily fraction of egg reduction,
               day0, dayf, daily fraction of larva reduction,
               day0, dayf, daily fraction of pupa reduction,
               day0, dayf, daily fraction of adult vector reduction]
     
Examples

     import albopictus as aa

     # The following line simulates the model at the first environmental grid point
     # of Bologna using the first parameter vector of the posterior mode Q4.
     sim = aa.vector.simPar(aa.clim['BO'][0],aa.param['Q4'][0])

     import pylab

     pylab.ylim([0,50])
     pylab.plot(aa.clim['BO'][0]['dates'],sim['colegg'])
     pylab.show()

chikv.simSpread

    Simulate the chikungunya transmission model for a given/simulated vector abundance dataset

    Parameters
    ----------

         clim:
              vector abundance dataset to be used
              it must be a dictionary with keys 'dates' and 'envar'
              clim['dates'] holds and array of dates corresponding to each data point with a length of final_time
              clim['envar'] should hold an array with the following data in the given order and size
                   adult abundance (coln4f) [final_time]
                   expected adult lifespan (cold4) [final_time]
                   standard deviation of adult lifespan (cold4s)  [final_time]
                   expected daily fecundity (colF4)  [final_time]
                   number of people [1]

         probs:
              n x n matrix of transmission probabilities between patches 
              (n: the number of patches)
              diagonal elements represent the daily probability of staying in a patch
              non diagonal entries will be normalised per row to yield a row sum of 1

         pr:
              an array of size vector.numpar holding the parameter vector

         cpr: (optional)
              an array of vector control parameters of the following order:
              [day0, dayf, daily fraction of breeding site reduction,
               day0, dayf, daily fraction of egg reduction,
               day0, dayf, daily fraction of larva reduction,
               day0, dayf, daily fraction of pupa reduction,
               day0, dayf, daily fraction of adult vector reduction]
     
Examples

     import numpy as np
     import albopictus as aa

     # The following line simulates the model at the first environmental grid point
     # of Bologna using the first parameter of the default posterior mode.
     sim = aa.vector.simPar(aa.clim['BO'][0],aa.vecpar[0])

     # This will generate the vector abundance dataset required for transmission simulation
     clim = [{
          'dates': aa.clim['BO'][0]['dates'],
          'envar': sim['coln4f'].tolist()+sim['cold4'].tolist()+sim['cold4s'].tolist()+sim['colF4'].tolist()+[100]
     }]

     # The following lines simulate chikv transmission for the simulated abundance data
     # using the first parameter of the default posterior mode.
     pr = np.array(aa.chikp[0]).copy()
     pr[aa.chikv.parids['introduce_time']] = 200
     simc = aa.chikv.simSpread(clim,[[1.0]],pr)

     import pylab

     pylab.xlim([aa.clim['BO'][0]['dates'][0],aa.clim['BO'][0]['dates'][500]])
     pylab.plot(aa.clim['BO'][0]['dates'],sim['coln4f'])
     pylab.plot(aa.clim['BO'][0]['dates'],simc[0]['colnInew'])
     pylab.show()

"""

__version__ = '2.0.4'

# modelAalbopictus - climateData ------------------------- //
from importlib_resources import files, as_file
import numpy
import json
from datetime import date

with (files(__package__) / "data" / "posterior.json").open("r", encoding="utf-8") as f:
    param = json.load(f)

with (files(__package__) / "data" / "prior03.json").open("r", encoding="utf-8") as f:
    prior03 = json.load(f)

with (files(__package__) / "data" / "prior08.json").open("r", encoding="utf-8") as f:
    prior08 = json.load(f)

with (files(__package__) / "data" / "prior13.json").open("r", encoding="utf-8") as f:
    prior13 = json.load(f)

with (files(__package__) / "data" / "prior18.json").open("r", encoding="utf-8") as f:
    prior18 = json.load(f)

with (files(__package__) / "data" / "priorAeae.json").open("r", encoding="utf-8") as f:
    priorae = json.load(f)

with (files(__package__) / "data" / "climate.json").open("r", encoding="utf-8") as f:
    clim = json.load(f)

for pr in clim:
    for clm in clim[pr]:
        clm['dates'] = [date.fromordinal(d) for d in clm['dates']]

provinces = ["Bologna","Ferrara","Modena","Piacenza","Parma","Ravenna","Reggio Emilia"]
prvn = ["BO","FE","MO","PC","PR","RA","RE"]

# modelAalbopictus --------------------------------------- //

from albopictus.readModel import prepareModel

with as_file(files(__package__) / "modelAalbopictus03.so") as model_path:
    vector03 = prepareModel(str(model_path))

with as_file(files(__package__) / "modelAalbopictus08.so") as model_path:
    vector08 = prepareModel(str(model_path))

with as_file(files(__package__) / "modelAalbopictus08b.so") as model_path:
    vector08b = prepareModel(str(model_path))

with as_file(files(__package__) / "modelAalbopictus08c.so") as model_path:
    vector08c = prepareModel(str(model_path))

with as_file(files(__package__) / "modelAalbopictus13.so") as model_path:
    vector13 = prepareModel(str(model_path))

with as_file(files(__package__) / "modelAalbopictus18.so") as model_path:
    vector18 = prepareModel(str(model_path))

with as_file(files(__package__) / "modelStochAalbopictus.so") as model_path:
    vectorst = prepareModel(str(model_path), "vectorst")

with as_file(files(__package__) / "modelStochAaegypti.so") as model_path:
    vectorae = prepareModel(str(model_path), "vectorae")


with as_file(files(__package__) / "modelStochCHIKV.so") as model_path:
    chikv = prepareModel(str(model_path), "chikv")

with as_file(files(__package__) / "modelStochCDZ.so") as model_path:
    cdz = prepareModel(str(model_path), "cdz")

# modelCulex --------------------------------------------- //

with as_file(files(__package__) / "modelCulex.so") as model_path:
    culex = prepareModel(str(model_path))

# modelStochSand - climateData --------------------------- //

with (files(__package__) / "data" / "priorSand.json").open("r", encoding="utf-8") as f:
    priorSand = json.load(f)

with as_file(files(__package__) / "modelStochSand.so") as model_path:
    sand = prepareModel(str(model_path), "sandfly")

# Set defaults ------------------------------------------- //

vector = vector08c
vecpar = param['Q4.a100+1']

chikv = chikv
chikp = param['QI']

prior = prior08

# Accessories -------------------------------------------- //

from albopictus.accessory import *

from albopictus.plotPos import plotPosterior

