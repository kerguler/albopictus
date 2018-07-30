from distutils.sysconfig import get_python_lib
import pkg_resources
import numpy

def ssig(x,a1,a2,a3):
    tmp = (a1)/(1.0+numpy.exp((a2)*((a3)-(x))))
    tmp[tmp<0.0] = 0.0
    tmp[tmp>1.0] = 1.0
    return tmp
def dsig(x,a1,a2,a3):
    tmp = (a1)/((1.0+numpy.exp((a2)-(x)))*(1.0+numpy.exp((x)-(a3))))
    tmp[tmp<0.0] = 0.0
    tmp[tmp>1.0] = 1.0
    return tmp
def dsig2(x,a1,a2,a3):
    tmp = (a1)/((1.0+numpy.exp((a2)-(x)))*(1.0+numpy.exp((x)-(a3))))
    tmp[tmp<0.0] = 0.0
    return tmp
def poly(x,a1,a2,a3):
    tmp = (a1) + (a2)*(x) + (a3)*pow((x),2)
    tmp[tmp<0.0] = 0.0
    return tmp
def expd(x,k):
    return numpy.exp(-k*x)
def flin(x,a1,a2):
    tmp = (a1) + (a2)*(x)
    tmp[tmp<0.0] = 0.0
    tmp[tmp>1.0] = 1.0
    return tmp
def fpow(x,a1,a2):
    tmp = (a1)*((x)**(a2))
    tmp[tmp<1.0] = 1.0
    return tmp

data = {'n23.comb.dens': numpy.genfromtxt(pkg_resources.resource_filename(__name__,"data/albopictus_n23_dev.txt"),delimiter="\t",skip_header=1,usecols=[1,2]),
        'n23.comb.temp': numpy.genfromtxt(pkg_resources.resource_filename(__name__,"data/albopictus_n23_temp.txt"),delimiter="\t",skip_header=1,usecols=[1,2]),
        'n23.surv': numpy.genfromtxt(pkg_resources.resource_filename(__name__,"data/albopictus_n23_surv.txt"),delimiter="\t",skip_header=1,usecols=[1,2]),
        'tbm.temp': numpy.genfromtxt(pkg_resources.resource_filename(__name__,"data/albopictus_tbm_temp.txt"),delimiter="\t",skip_header=1,usecols=[1,2]),
        'F4.temp': numpy.genfromtxt(pkg_resources.resource_filename(__name__,"data/albopictus_F4_temp.txt"),delimiter="\t",skip_header=1,usecols=[1,2]),
        'd1.temp': numpy.genfromtxt(pkg_resources.resource_filename(__name__,"data/albopictus_d1_temp.txt"),delimiter="\t",skip_header=1,usecols=[1,2]),
        'd2.temp': numpy.genfromtxt(pkg_resources.resource_filename(__name__,"data/albopictus_d2_temp.txt"),delimiter="\t",skip_header=1,usecols=[1,2]),
        'd3.temp': numpy.genfromtxt(pkg_resources.resource_filename(__name__,"data/albopictus_d3_temp.txt"),delimiter="\t",skip_header=1,usecols=[1,2]),
        'p0.temp': numpy.genfromtxt(pkg_resources.resource_filename(__name__,"data/albopictus_p0_temp.txt"),delimiter="\t",skip_header=1,usecols=[1,2]),
        'p1.temp': numpy.genfromtxt(pkg_resources.resource_filename(__name__,"data/albopictus_p1_temp.txt"),delimiter="\t",skip_header=1,usecols=[1,2]),
        'p2.temp': numpy.genfromtxt(pkg_resources.resource_filename(__name__,"data/albopictus_p2_temp.txt"),delimiter="\t",skip_header=1,usecols=[1,2]),
        'p3.temp': numpy.genfromtxt(pkg_resources.resource_filename(__name__,"data/albopictus_p3_temp.txt"),delimiter="\t",skip_header=1,usecols=[1,2]),
        'd4.temp': numpy.genfromtxt(pkg_resources.resource_filename(__name__,"data/albopictus_d4_temp.txt"),delimiter="\t",skip_header=1,usecols=[1,2])}

def readDataType(filename,begin=0):
    dat = numpy.genfromtxt(filename,delimiter="\t",skip_header=1,usecols=[0],dtype='S')
    count = begin
    for r in range(len(dat)):
        if dat[r] not in dat[:r]:
            dat[numpy.where(dat[r]==dat)[0]] = count
            count += 1
    return numpy.array(dat,dtype=numpy.int)

priors = {
    'n23.comb': ['n23.1', 'n23.2', 'n23.3', 'n23.4', 'n23.5'],
    'n23.surv': ['n23.surv'],
    'tbm.temp': ['tbm.1','tbm.2','tbm.3'],
    'F4.temp': ['F4.1', 'F4.2', 'F4.3'],
    'd1.temp': ['d1.1', 'd1.2', 'd1.3'],
    'd2.temp': ['d2.1', 'd2.2', 'd2.3'],
    'd3.temp': ['d3.1', 'd3.2', 'd3.3'],
    'p0.temp': ['p0.1', 'p0.2'],
    'p1.temp': ['p1.1', 'p1.2', 'p1.3'],
    'p2.temp': ['p2.1', 'p2.2', 'p2.3'],
    'p3.temp': ['p3.1', 'p3.2', 'p3.3'],
    'd4.temp': ['d4.1', 'd4.2', 'd4.3']
    }

def getFits(model):
    return {
    'n23.comb.dens': {'best':[],'stats':[],'pop':[],'pri':numpy.array([model.parids[x] for x in ['n23.1', 'n23.2', 'n23.3', 'n23.4', 'n23.5']])},
    'n23.comb.temp': {'best':[],'stats':[],'pop':[],'pri':numpy.array([model.parids[x] for x in ['n23.1', 'n23.2', 'n23.3', 'n23.4', 'n23.5']])},
    'n23.surv': {'best':[],'stats':[],'pop':[],'pri':numpy.array([model.parids[x] for x in ['n23.surv']])},
    'tbm.temp': {'best':[],'stats':[],'pop':[],'pri':numpy.array([model.parids[x] for x in ['tbm.1','tbm.2','tbm.3']])},
    'F4.temp': {'best':[],'stats':[],'pop':[],'pri':numpy.array([model.parids[x] for x in ['F4.1', 'F4.2', 'F4.3']])},
    'd1.temp': {'best':[],'stats':[],'pop':[],'pri':numpy.array([model.parids[x] for x in ['d1.1', 'd1.2', 'd1.3']])},
    'd2.temp': {'best':[],'stats':[],'pop':[],'pri':numpy.array([model.parids[x] for x in ['d2.1', 'd2.2', 'd2.3']])},
    'd3.temp': {'best':[],'stats':[],'pop':[],'pri':numpy.array([model.parids[x] for x in ['d3.1', 'd3.2', 'd3.3']])},
    'p0.temp': {'best':[],'stats':[],'pop':[],'pri':numpy.array([model.parids[x] for x in ['p0.1', 'p0.2']])},
    'p1.temp': {'best':[],'stats':[],'pop':[],'pri':numpy.array([model.parids[x] for x in ['p1.1', 'p1.2', 'p1.3']])},
    'p2.temp': {'best':[],'stats':[],'pop':[],'pri':numpy.array([model.parids[x] for x in ['p2.1', 'p2.2', 'p2.3']])},
    'p3.temp': {'best':[],'stats':[],'pop':[],'pri':numpy.array([model.parids[x] for x in ['p3.1', 'p3.2', 'p3.3']])},
    'd4.temp': {'best':[],'stats':[],'pop':[],'pri':numpy.array([model.parids[x] for x in ['d4.1', 'd4.2', 'd4.3']])}
    }

datatype = {'n23.comb.dens': readDataType(pkg_resources.resource_filename(__name__,"data/albopictus_n23_dev.txt")),
            'n23.comb.temp': readDataType(pkg_resources.resource_filename(__name__,"data/albopictus_n23_temp.txt"),1+max(readDataType(pkg_resources.resource_filename(__name__,"data/albopictus_n23_dev.txt")))),
            'n23.surv': readDataType(pkg_resources.resource_filename(__name__,"data/albopictus_n23_surv.txt")),
            'tbm.temp': readDataType(pkg_resources.resource_filename(__name__,"data/albopictus_tbm_temp.txt")),
            'F4.temp': readDataType(pkg_resources.resource_filename(__name__,"data/albopictus_F4_temp.txt")),
            'd1.temp': readDataType(pkg_resources.resource_filename(__name__,"data/albopictus_d1_temp.txt")),
            'd2.temp': readDataType(pkg_resources.resource_filename(__name__,"data/albopictus_d2_temp.txt")),
            'd3.temp': readDataType(pkg_resources.resource_filename(__name__,"data/albopictus_d3_temp.txt")),
            'p0.temp': readDataType(pkg_resources.resource_filename(__name__,"data/albopictus_p0_temp.txt")),
            'p1.temp': readDataType(pkg_resources.resource_filename(__name__,"data/albopictus_p1_temp.txt")),
            'p2.temp': readDataType(pkg_resources.resource_filename(__name__,"data/albopictus_p2_temp.txt")),
            'p3.temp': readDataType(pkg_resources.resource_filename(__name__,"data/albopictus_p3_temp.txt")),
            'd4.temp': readDataType(pkg_resources.resource_filename(__name__,"data/albopictus_d4_temp.txt"))}
    
fun = { 'n23.comb.dens': fpow,
        'n23.comb.temp': poly,
        'n23.surv': expd,
        'tbm.temp': poly,
        'F4.temp': poly,
        'd1.temp': poly,
        'd2.temp': poly,
        'd3.temp': poly,
        'p0.temp': flin,
        'p1.temp': dsig,
        'p2.temp': dsig,
        'p3.temp': dsig,
        'd4.temp': dsig2 }
        
tasks = {
    'n23.comb.dens': {
        'prior': 'n23.comb',
        'fitness': lambda pr: (numpy.Inf if any(pr<-2) or any(pr>2) or pr[1]<0 or pr[4]>0 else sum(numpy.hstack([((data['n23.comb.dens'][:,1]-fun['n23.comb.dens'](data['n23.comb.dens'][:,0],pr[0],pr[1]*fun['n23.comb.temp'](numpy.array([28.0]),pr[2],pr[3],pr[4])[0]))/0.1)**2,((data['n23.comb.temp'][:,1]-fun['n23.comb.temp'](data['n23.comb.temp'][:,0],pr[2],pr[3],pr[4]))/0.1)**2]))),
        'fun': lambda x,a1,a2,a3,a4,a5: (fun['n23.comb.dens'](x,a1,a2*fun['n23.comb.temp'](numpy.array([28.0]),a3,a4,a5)[0])),
        'numpar': 5,
        'domain': {'x':numpy.arange(0,5000,1),'y':[0,4]},
        'xlabel': "Density of immature stages (per Litre)",
        'ylabel': "Fractional change in development"
        },
    'n23.comb.temp': {
        'prior': 'n23.comb',
        'fitness': lambda pr: (numpy.Inf if any(pr<-2) or any(pr>2) or pr[1]<0 or pr[4]>0 else sum(numpy.hstack([((data['n23.comb.dens'][:,1]-fun['n23.comb.dens'](data['n23.comb.dens'][:,0],pr[0],pr[1]*fun['n23.comb.temp'](numpy.array([28.0]),pr[2],pr[3],pr[4])[0]))/0.1)**2,((data['n23.comb.temp'][:,1]-fun['n23.comb.temp'](data['n23.comb.temp'][:,0],pr[2],pr[3],pr[4]))/0.1)**2]))),
        'fun': lambda x,a1,a2,a3,a4,a5: (fun['n23.comb.temp'](x,a3,a4,a5)),
        'numpar': 5,
        'domain': {'x':numpy.arange(-20,60,1),'y':[0,3]},
        'xlabel': "Water temperature ($^\circ$C)",
        'ylabel': "Fractional change in development"
        },
    'n23.surv': {
        'prior': 'n23.surv',
        'fitness': lambda pr: (numpy.Inf if pr[0]<0 or pr[0]>1 else sum(((data['n23.surv'][:,1]-fun['n23.surv'](data['n23.surv'][:,0],pr[0]))/0.01)**2)),
        'numpar': 1,
        'domain': {'x':numpy.arange(0,5000,1),'y':[-0.05,1.05]},
        'xlabel': "Density of immature stages (per Litre)",
        'ylabel': "Larva and pupa survival"
        },
    'tbm.temp': {
        'prior': 'tbm.temp',
        'fitness': lambda pr: (numpy.Inf if pr[2]<0.0 or any(pr<-30) or any(pr>30) else sum(((data['tbm.temp'][:,1]-fun['tbm.temp'](data['tbm.temp'][:,0],pr[0],pr[1],pr[2]))/1.0)**2)),
        'numpar': 3,
        'domain': {'x':numpy.arange(-20,60,1),'y':[-1,20]},
        'xlabel': "Air temperature ($^\circ$C)",
        'ylabel': "Days to first blood meal"
        },
    'F4.temp': {
        'prior': 'F4.temp',
        'fitness': lambda pr: (numpy.Inf if pr[2]>-0.01 or any(pr<-30) or any(pr>30) else sum(((data['F4.temp'][:,1]-fun['F4.temp'](data['F4.temp'][:,0],pr[0],pr[1],pr[2]))/1.0)**2)),
        'numpar': 3,
        'domain': {'x':numpy.arange(-20,60,1),'y':[-1,30]},
        'xlabel': "Air temperature ($^\circ$C)",
        'ylabel': "Fecundity (eggs/female/day)"
        },
    'd1.temp': {
        'prior': 'd1.temp',
        'fitness': lambda pr: (numpy.Inf if pr[2]<0 or any(pr<-25) or any(pr>25) else sum(((data['d1.temp'][:,1]-fun['d1.temp'](data['d1.temp'][:,0],pr[0],pr[1],pr[2]))/1.0)**2)),
        'numpar': 3,
        'domain': {'x':numpy.arange(-20,60,1),'y':[-1,15]},
        'xlabel': "Water temperature ($^\circ$C)",
        'ylabel': "Egg development time (days)"
        },
    'd2.temp': {
        'prior': 'd2.temp',
        'fitness': lambda pr: (numpy.Inf if pr[2]<0 or any(pr<-100) or any(pr>100) else sum(((data['d2.temp'][:,1]-fun['d2.temp'](data['d2.temp'][:,0],pr[0],pr[1],pr[2]))/1.0)**2)),
        'numpar': 3,
        'domain': {'x':numpy.arange(-20,60,1),'y':[-1,40]},
        'xlabel': "Water temperature ($^\circ$C)",
        'ylabel': "Larva development time (days)"
        },
    'd3.temp': {
        'prior': 'd3.temp',
        'fitness': lambda pr: (numpy.Inf if pr[2]<0 or any(pr<-25) or any(pr>25) else sum(((data['d3.temp'][:,1]-fun['d3.temp'](data['d3.temp'][:,0],pr[0],pr[1],pr[2]))/1.0)**2)),
        'numpar': 3,
        'domain': {'x':numpy.arange(-20,60,1),'y':[-1,20]},
        'xlabel': "Water temperature ($^\circ$C)",
        'ylabel': "Pupa development time (days)"
        },
    'p0.temp': {
        'prior': 'p0.temp',
        'fitness': lambda pr: (numpy.Inf if pr[0]<-1 or pr[0]>1 or pr[1]<-1 or pr[1]>1 else sum(((data['p0.temp'][:,1]-fun['p0.temp'](data['p0.temp'][:,0],pr[0],pr[1]))/0.05)**2)),
        'numpar': 2,
        'domain': {'x':numpy.arange(-20,60,1),'y':[-0.05,1.05]},
        'xlabel': "Air temperature ($^\circ$C)",
        'ylabel': "Daily survival proportion (diapausing egg)"
        },
    'p1.temp': {
        'prior': 'p1.temp',
        'fitness': lambda pr: (numpy.Inf if pr[0]<0 or pr[0]>1 or pr[1]<-100 or pr[1]>100 or pr[2]<-100 or pr[2]>100 else sum(((data['p1.temp'][:,1]-fun['p1.temp'](data['p1.temp'][:,0],pr[0],pr[1],pr[2]))/0.1)**2)),
        'numpar': 3,
        'domain': {'x':numpy.arange(-20,60,1),'y':[-0.05,1.05]},
        'xlabel': "Water temperature ($^\circ$C)",
        'ylabel': "Daily survival proportion (non-diapausing egg)"
        },
    'p2.temp': {
        'prior': 'p2.temp',
        'fitness': lambda pr: (numpy.Inf if pr[0]<0 or pr[0]>1 or pr[1]<-100 or pr[1]>100 or pr[2]<-100 or pr[2]>100 else sum(((data['p2.temp'][:,1]-fun['p2.temp'](data['p2.temp'][:,0],pr[0],pr[1],pr[2]))/0.1)**2)),
        'numpar': 3,
        'domain': {'x':numpy.arange(-20,60,1),'y':[-0.05,1.05]},
        'xlabel': "Water temperature ($^\circ$C)",
        'ylabel': "Daily survival proportion (larva)"
        },
    'p3.temp': {
        'prior': 'p3.temp',
        'fitness': lambda pr: (numpy.Inf if pr[0]<0 or pr[0]>1 or pr[1]<-100 or pr[1]>100 or pr[2]<-100 or pr[2]>100 else sum(((data['p3.temp'][:,1]-fun['p3.temp'](data['p3.temp'][:,0],pr[0],pr[1],pr[2]))/0.1)**2)),
        'numpar': 3,
        'domain': {'x':numpy.arange(-20,60,1),'y':[-0.05,1.05]},
        'xlabel': "Water temperature ($^\circ$C)",
        'ylabel': "Daily survival proportion (pupa)"
        },
    'd4.temp': {
        'prior': 'd4.temp',
        'fitness': lambda pr: (numpy.Inf if pr[0]<0 or pr[0]>100 or pr[1]<-100 or pr[1]>100 or pr[2]<-100 or pr[2]>100 else sum(((data['d4.temp'][:,1]-fun['d4.temp'](data['d4.temp'][:,0],pr[0],pr[1],pr[2]))/5.0)**2)),
        'numpar': 3,
        'domain': {'x':numpy.arange(-20,60,1),'y':[0,80]},
        'xlabel': "Air temperature ($^\circ$C)",
        'ylabel': "Adult survival time (days)"
        }
}
