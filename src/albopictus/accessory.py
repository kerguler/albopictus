import numpy

def calcEnsemble(sim):
    if len(sim)==0:
        return None
    p = numpy.percentile(sim,[2.5,50,97.5],axis=0)
    return {
        'lower':p[0],
        'median':p[1],
        'higher':p[2],
        'mean':numpy.mean(sim,axis=0),
        'std':numpy.std(sim,axis=0)
        }
