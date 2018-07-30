import ctypes
import numpy
import numpy.ctypeslib as npct
array_1d_double = npct.ndpointer(dtype=numpy.float64, ndim=1, flags='CONTIGUOUS')
array_1d_int = npct.ndpointer(dtype=numpy.int32, ndim=1, flags='CONTIGUOUS')

from scipy.stats import gamma
from albopictus.accessory import calcEnsemble, calc_ss

class prepareModel:
    def __init__(self,modelname,label):
        global exit_gamma
        global exit_rng_setup
        global exit_rng_destroy
        #
        import atexit
        self.modelname = modelname
        self.model = ctypes.cdll.LoadLibrary(self.modelname)
        #
        try:
            self.set_sim_mode = self.model.set_sim_mode
            self.set_sim_mode.restype = None
            self.set_sim_mode.argtypes = [npct.c_char]
        except:
            pass
        #
        try:
            self.set_gamma_mode = self.model.set_gamma_mode
            self.set_gamma_mode.restype = None
            self.set_gamma_mode.argtypes = [npct.c_char]
        except:
            pass
        #
        try:
            atexit.register(self.model.gamma_mean_destroy)
        except:
            pass
        #
        try:
            self.model.rng_setup.restype = None
            self.model.rng_setup.argtypes = [ctypes.c_char_p]
            self.model.rng_setup(label)
        except:
            pass
        try:
            atexit.register(self.model.rng_destroy)
        except:
            pass
        #
        self.numparModel = self.model.numparModel
        self.numparModel.restype = None
        self.numparModel.argtypes = [array_1d_int,array_1d_int]
        self.numpar = numpy.arange(1,dtype=numpy.int32)
        self.nummet = numpy.arange(1,dtype=numpy.int32)
        ret = self.numparModel(self.numpar,self.nummet)
        self.numpar = self.numpar[0]
        self.nummet = self.nummet[0]
        #
        try:
            self.param_model = self.model.param_model
            self.param_model.restype = None
            self.param_model.argtypes = [ctypes.POINTER(ctypes.c_char_p),array_1d_double]
            temp = (ctypes.c_char_p * (self.nummet+self.numpar))(256)
            self.param = numpy.array([0 for n in range(self.numpar)],dtype=numpy.float64)
            ret = self.param_model(temp,self.param)
            temp = numpy.array([str(elm,'utf-8') for elm in temp])
            self.metnames = numpy.copy(temp[:self.nummet])
            self.parnames = numpy.copy(temp[-self.numpar:])
        except:
            print("Skipping default parameters")
            self.metnames = numpy.array(["coln%d" %(n) for n in range(self.nummet)])
            self.parnames = numpy.array(["par%d" %(n) for n in range(self.numpar)])
        #
        self.metids = {}
        for elm in self.metnames:
            self.metids[elm] = numpy.where(elm==self.metnames)[0][0]
        #
        self.parids = {}
        for elm in self.parnames:
            self.parids[elm] = numpy.where(elm==self.parnames)[0][0]
        #
        try:
            self.sim_model = self.model.sim_model
            self.sim_model.restype = None
            self.sim_model.argtypes = [array_1d_double,
                                    array_1d_double,
                                    array_1d_int,
                                    array_1d_int,
                                    array_1d_double,
                                    array_1d_int]
        except:
            pass
        #
        try:
            self.sim_spread = self.model.sim_spread
            self.sim_spread.restype = None
            self.sim_spread.argtypes = [array_1d_double,
                                        array_1d_double,
                                        array_1d_double,
                                        array_1d_int,
                                        array_1d_int,
                                        array_1d_int,
                                        array_1d_double,
                                        array_1d_int]
        except:
            pass
        #
        try:
            self.model.gamma_pdf
            self.model.gamma_pdf.restype = ctypes.c_double
            self.model.gamma_pdf.argtypes = [ctypes.c_double,
                                             ctypes.c_double,
                                             ctypes.c_double]
        except:
            pass
        #
        try:
            self.model.gamma_dist_prob
            self.model.gamma_dist_prob.restype = ctypes.c_double
            self.model.gamma_dist_prob.argtypes = [ctypes.c_double,
                                                   ctypes.c_double,
                                                   ctypes.c_double]
        except:
            pass
    #
    def gamma_pdf(self,n,mean,sd):
        """
        Probability density function of the gamma distribution with given mean and standard deviation
        """
        if not hasattr(n, "__iter__"):
            return self.model.gamma_pdf(n,mean,sd)
        return [self.model.gamma_pdf(nn,mean,sd) for nn in n]
    #
    def gamma_prob(self,n,mean,sd):
        """
        Probability of daily development/survival with given mean and standard deviation
        """
        if not hasattr(n, "__iter__"):
            return self.model.gamma_dist_prob(mean,sd,n)
        return [self.model.gamma_dist_prob(mean,sd,nn) for nn in n]
    #
    def simPar(self,clim,pr,cpr=[]):
        """
        Main simulation routine
        """
        fT = numpy.array(len(clim['dates']),dtype=numpy.int32,ndmin=1)
        envar = numpy.array(clim['envar'],dtype=numpy.float64)
        param = numpy.array(pr,dtype=numpy.float64)
        param = numpy.hstack([param,cpr])
        control = numpy.array(len(cpr)!=0,dtype=numpy.int32,ndmin=1)
        result = numpy.ndarray((self.nummet+1)*fT[0],dtype=numpy.float64)
        success = numpy.array(0,dtype=numpy.int32,ndmin=1)
        ret = self.sim_model(envar,
                             param,
                             fT,
                             control,
                             result,
                             success)
        ret = {
            'colT':result[0:fT[0]],
            'success':success
            }
        for n in range(self.nummet):
            ret[self.metnames[n]] = result[((n+1)*fT[0]):((n+2)*fT[0])]
        return ret
    # 
    def simParMat(self,clim,parmat,boil=None,control=[],ensemble=False):
        """
        Wrapper routine for simPar with a list of parameters sampled from the posterior
        """
        if boil is None:
            boil = {}
            for elm in self.metnames:
                boil[elm] = []
        for pr in parmat:
            s = self.simPar(clim,pr,control)
            if (s.__class__.__name__!='dict' or not s['success']):
                return {'success':False}
            for elm in boil.keys():
                boil[elm].append(s[elm].copy())
        if ensemble:
            for elm in boil.keys():
                boil[elm] = calcEnsemble(boil[elm])
        else:
            for elm in boil.keys():
                boil[elm] = numpy.mean(boil[elm],axis=0)
        boil['colT'] = s['colT'].copy()
        boil['success'] = True
        return boil
    # for the calculation of HSI
    def simParPosterior(self,clim,parmat,weeksBefore=52):
        """
        Wrapper routine to calculate RSI and HSI
        """
        sims = [sim['coln4f'][(7*weeksBefore):] for sim in [self.simPar(clim,prm) for prm in parmat] if sim['success'][0]==1]
        if len(sims) == 0:
            return {'sim':numpy.nan,
                    'mean':numpy.nan,
                    'sindex':numpy.nan}
        elif len(sims) == 1:
            return {'sim':numpy.array(sims),
                    'mean':numpy.array(sims),
                    'sindex':numpy.mean(sims)}
        ret = {}
        ret['sim'] = numpy.array(sims)
        ret['mean'] = numpy.mean(ret['sim'],axis=0)
        ret['sindex'] = numpy.mean(ret['mean'])
        return ret
    # for Chikungunya
    def calcTProbs(self,probs,numreg):
        """
        Prepares the propensity matrix for transmissions between regions
        """
        tprobs = numpy.array(probs,dtype=numpy.float64).copy()
        mask = numpy.diag(numpy.repeat(True,numreg))
        for i in range(numreg):
            tprobs[i][~mask[i]] = (1.0-tprobs[i][mask[i]]) * tprobs[i][~mask[i]] / numpy.sum(tprobs[i][~mask[i]])
        return numpy.concatenate(tprobs)
    #
    def simSpread(self,clim,probs,pr,cpr=[]):
        """
        Main spatiotemporal simulation routine
        """
        numreg = numpy.array(len(clim),dtype=numpy.int32,ndmin=1)
        tprobs = numpy.concatenate(probs)
        envar = numpy.concatenate([clm['envar'] for clm in clim])
        fT = numpy.array(len(clim[0]['dates']),dtype=numpy.int32,ndmin=1)
        param = numpy.array(pr,dtype=numpy.float64)
        param = numpy.hstack([param,cpr])
        control = numpy.array(len(cpr)!=0,dtype=numpy.int32,ndmin=1)
        result = numpy.ndarray((self.nummet+1)*fT[0]*numreg[0],dtype=numpy.float64)
        success = numpy.ndarray(numreg[0],dtype=numpy.int32)
        ret = self.sim_spread(envar,
                              tprobs,
                              param,
                              fT,
                              control,
                              numreg,
                              result,
                              success)
        l = fT[0]*(self.nummet+1)
        ret = [{
            'colT':result[(m*l):(m*l+fT[0])],
            'success':success[m]
            } for m in range(numreg[0])]
        for m in range(numreg[0]):
            for n in range(self.nummet):
                ret[m][self.metnames[n]] = result[(m*l+(n+1)*fT[0]):(m*l+(n+2)*fT[0])]
        return ret
    # Prior distribution score function for all
    def scorePar(self,pr,prior,verbose=False):
        """
        Calculates prior probability for a list of parameter values and a given prior definition

        Parameters
        ----------

              pr:
                    List of parameter values for the model (must have numpar elements)
              prior:
                    Dictionary of prior probability definitions. An example is given below.
                    {'Critical temperature': {
                         'parids': ['PP.ta.thr'], 
                         'mean': [21.0], 
                         'var': [[9.0]], 
                         'var.inv': [[1.0/9.0]],
                         'min': 0.0,
                         'max': 100.0
                    } }
              verbose:
                    Warn about unacceptable values (Default: False)
        """
        scr = 0
        for key in prior:
            if 'min' in prior[key] and any([pr[self.parids[pid]] < prior[key]['min'] for pid in prior[key]['parids']]):
                if verbose:
                    print("Warning: Parameter %s breached minimum allowance!" %(key))
                return numpy.Inf
            if 'max' in prior[key] and any([pr[self.parids[pid]] > prior[key]['max'] for pid in prior[key]['parids']]):
                if verbose:
                    print("Warning: Parameter %s breached maximum allowance!" %(key))
                return numpy.Inf
            if 'mean' in prior[key] and 'var' in prior[key] and 'var.inv' in prior[key]:
                scr += 0.5*calc_ss(numpy.array([pr[self.parids[x]] for x in prior[key]['parids']]),
                                       prior[key]['mean'],
                                       prior[key]['var'],
                                       prior[key]['var.inv'])
        return scr
    # Simulate a parameter vector form the prior distribution
    def simPrior(self,pr,prior):
        """
        Simulates a parameter vector from a given prior definition

        Parameters
        ----------

              pr:
                    The initial list of parameter values (must have numpar elements)
              prior:
                    Dictionary of prior probability definitions. An example is given below.
                    {'Critical temperature': {
                         'parids': ['PP.ta.thr'], 
                         'mean': [21.0], 
                         'var': [[9.0]], 
                         'var.inv': [[1.0/9.0]]
                    } }
        """
        scr = numpy.Inf
        while scr<0 or scr==numpy.Inf:
            for key in prior:
                # simulate
                ids = [self.parids[x] for x in prior[key]['parids']]
                r = []
                if 'mean' in prior[key] and 'var' in prior[key] and 'var.inv' in prior[key]:
                    r = numpy.random.multivariate_normal(prior[key]['mean'],prior[key]['var'],1)[0]
                elif 'min' in prior[key] and 'max' in prior[key]:
                    r = numpy.random.uniform(prior[key]['min'],prior[key]['max'],len(prior[key]['parids']))
                if len(r)>0:
                    for i in numpy.arange(len(ids)):
                        pr[ids[i]] = r[i]
            # check
            scr = self.scorePar(pr,prior)
