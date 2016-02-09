import ctypes
import numpy
import numpy.ctypeslib as npct
array_1d_double = npct.ndpointer(dtype=numpy.float64, ndim=1, flags='CONTIGUOUS')
array_1d_int = npct.ndpointer(dtype=numpy.int32, ndim=1, flags='CONTIGUOUS')

class prepareModel:
    def __init__(self,modelname):
        self.modelname = modelname
        self.model = ctypes.cdll.LoadLibrary(self.modelname)
        #
        import atexit
        try:
            self.model.rng_setup()
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
        self.sim_model = self.model.sim_model
        self.sim_model.restype = None
        self.sim_model.argtypes = [array_1d_double,
                                   array_1d_double,
                                   array_1d_int,
                                   array_1d_double,
                                   array_1d_int]
    #
    def simPar(self,clim,pr):
        fT = numpy.array(len(clim['dates']),dtype=numpy.int32,ndmin=1)
        envar = numpy.array(clim['envar'],dtype=numpy.float64)
        param = numpy.array(pr,dtype=numpy.float64)
        result = numpy.ndarray((self.nummet+1)*fT[0],dtype=numpy.float64)
        success = numpy.array(0,dtype=numpy.int32,ndmin=1)
        ret = self.sim_model(envar,
                             param,
                             fT,
                             result,
                             success)
        ret = {
            'colT':result[0:fT[0]],
            'success':success
            }
        for n in range(self.nummet):
            ret[self.metnames[n]] = result[((n+1)*fT[0]):((n+2)*fT[0])]
        return ret
