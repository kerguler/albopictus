import ctypes
import numpy
import numpy.ctypeslib as npct
array_1d_double = npct.ndpointer(dtype=numpy.float64, ndim=1, flags='CONTIGUOUS')
array_1d_int = npct.ndpointer(dtype=numpy.int32, ndim=1, flags='CONTIGUOUS')

from scipy.stats import gamma

exit_gamma = False
exit_rng_setup = False
exit_rng_destroy = False

class prepareModel:
    def __init__(self,modelname):
        global exit_gamma
        global exit_rng_setup
        global exit_rng_destroy
        #
        import atexit
        self.modelname = modelname
        self.model = ctypes.cdll.LoadLibrary(self.modelname)
        #
        try:
            self.set_gamma_mode = self.model.set_gamma_mode
            self.set_gamma_mode.restype = None
            self.set_gamma_mode.argtypes = [npct.c_char]
        except:
            pass
        #
        try:
            if not exit_gamma:
                atexit.register(self.model.gamma_mean_destroy)
                exit_gamma = True
        except:
            pass
        #
        try:
            if not exit_rng_setup:
                self.model.rng_setup()
                exit_rng_setup = True
        except:
            pass
        try:
            if not exit_rng_destroy:
                atexit.register(self.model.rng_destroy)
                exit_rng_destroy = True
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
            temp = numpy.array([elm for elm in temp])
            self.metnames = numpy.copy(temp[:self.nummet])
            self.parnames = numpy.copy(temp[-self.numpar:])
        except:
            print "Skipping default parameters"
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
                                        array_1d_double,
                                        array_1d_int]
        except:
            pass
    #
    def simPar(self,clim,pr):
        fT = numpy.array(len(clim['dates']),dtype=numpy.int32,ndmin=1)
        envar = numpy.array(clim['envar'],dtype=numpy.float64)
        param = numpy.array(pr,dtype=numpy.float64)
        control = numpy.array(0,dtype=numpy.int32,ndmin=1)
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
    # for albopictus
    def simControl(self,clim,pr,cpr):
        fT = numpy.array(len(clim['dates']),dtype=numpy.int32,ndmin=1)
        envar = numpy.array(clim['envar'],dtype=numpy.float64)
        param = numpy.array(pr,dtype=numpy.float64)
        param = numpy.hstack([param,cpr])
        control = numpy.array(1,dtype=numpy.int32,ndmin=1)
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
    # for Chikungunya
    def calcTProbs(self,probs,numreg):
        tprobs = numpy.array(probs,dtype=numpy.float64).copy()
        mask = numpy.diag(numpy.repeat(True,numreg))
        for i in range(numreg):
            tprobs[i][~mask[i]] = (1.0-tprobs[i][mask[i]]) * tprobs[i][~mask[i]] / numpy.sum(tprobs[i][~mask[i]])
        for i in range(numreg):
            tprobs[:,i] /= numpy.sum(tprobs[:,i:],axis=1)
        return numpy.concatenate(tprobs)
    #
    def simSpread(self,clim,probs,pr):
        numreg = numpy.array(len(clim),dtype=numpy.int32,ndmin=1)
        tprobs = self.calcTProbs(probs,numreg)
        envar = numpy.concatenate([clm['envar'] for clm in clim])
        fT = numpy.array(len(clim[0]['dates']),dtype=numpy.int32,ndmin=1)
        param = numpy.array(pr,dtype=numpy.float64)
        result = numpy.ndarray((self.nummet+1)*fT[0]*numreg[0],dtype=numpy.float64)
        success = numpy.ndarray(numreg[0],dtype=numpy.int32)
        ret = self.sim_spread(envar,
                              tprobs,
                              param,
                              fT,
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
