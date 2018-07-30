import albopictus.setPrior
import numpy

markers = numpy.array(["o","v","^","<",">","1","2","3","4","8","s","p","*","h","H","+","x","D","d","|","_"])
rnorm = numpy.random.multivariate_normal
extension = 'png'

############################################################
# ---
############################################################

def pplot(model,prior,parmat,options):
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    import matplotlib.patches as patches
    from scipy.stats import norm
    # ---
    parmat = numpy.array(parmat)
    print("Plotting %s..." %(options['name']))
    if not options['xlim']:
        if len(parmat)>0:
            options['xlim'] = [numpy.min(parmat[:,model.parids[options['par_name']]]*options['scale']),numpy.max(parmat[:,model.parids[options['par_name']]]*options['scale'])]
        else:
            options['xlim'] = [0,1]
        if options['xlim'][0]==options['xlim'][1]:
            options['xlim'] = [options['xlim'][0]-1,options['xlim'][1]+1]
    x = numpy.arange(options['xlim'][0],options['xlim'][1],(options['xlim'][1]-options['xlim'][0])/1000.0)
    matplotlib.rcParams.update({'font.size': 20})
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left()
    ax.set_xlim(options['xlim'])
    if options['prior_name']:
        ax.plot(x,norm.pdf(x,loc=prior[options['prior_name']]['mean'][0],scale=numpy.sqrt(prior[options['prior_name']]['var'][0][0])),color='blue',alpha=1,lw=3,zorder=0)
    if len(parmat)>0:
        a = ax.hist(parmat[:,model.parids[options['par_name']]]*options['scale'],normed=True,color='red',alpha=0.5,lw=1,zorder=0,label=options['name'])
    else:
        a = [0,1]
    if not options['ylim']:
        options['ylim'] = [numpy.min(a[0]),numpy.max(a[0])]
        if options['ylim'][0]==options['ylim'][1]:
            options['ylim'] = [options['ylim'][0]-1,options['ylim'][1]+1]
    ax.set_ylim(options['ylim'])
    ax.set_aspect(numpy.float(65.0*(options['xlim'][1]-options['xlim'][0])/(options['ylim'][1]-options['ylim'][0])/80.0))
    plt.xlabel(options['xlab'])
    plt.ylabel(options['ylab'])
    plt.savefig("%s/figures.%s.png" %(options['directory'],options['name']),bbox_inches="tight",transparent=True)
    plt.close()

############################################################

def plot_fit(parmat,plts,fts,tsk):
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    import matplotlib.patches as patches
    from scipy.stats import norm
    # ---
    parmat = numpy.array(parmat)
    f = plts['fun'] if 'fun' in plts.keys() else setPrior.fun[tsk]
    x = plts['domain']['x']
    #
    ret = []
    for pr in parmat[:,fts['pri']]:
        if plts['numpar']==1:
            ret.append(f(x,pr[0]))
        elif plts['numpar']==2:
            ret.append(f(x,pr[0],pr[1]))
        elif plts['numpar']==3:
            ret.append(f(x,pr[0],pr[1],pr[2]))
        elif plts['numpar']==4:
            ret.append(f(x,pr[0],pr[1],pr[2],pr[3]))
        else:
            ret.append(f(x,pr[0],pr[1],pr[2],pr[3],pr[4]))
    fts['pop'] = ret

def plot_pos(prior,plts,dat,typ,directory,name,fts,tsk):
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    import matplotlib.patches as patches
    from scipy.stats import norm
    # ---
    f = plts['fun'] if 'fun' in plts.keys() else setPrior.fun[tsk]
    x = plts['domain']['x']
    ret = []
    for pr in rnorm(prior[plts['prior']]['mean'],prior[plts['prior']]['var'],1000):
        if plts['numpar']==1:
            ret.append(f(x,pr[0]))
        elif plts['numpar']==2:
            ret.append(f(x,pr[0],pr[1]))
        elif plts['numpar']==3:
            ret.append(f(x,pr[0],pr[1],pr[2]))
        elif plts['numpar']==4:
            ret.append(f(x,pr[0],pr[1],pr[2],pr[3]))
        else:
            ret.append(f(x,pr[0],pr[1],pr[2],pr[3],pr[4]))
    res = numpy.percentile(ret,[2.5,50,97.5],axis=0)
    matplotlib.rcParams.update({'font.size': 20})
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left()
    ylim = plts['domain']['y']
    xlim = [min(plts['domain']['x']),max(plts['domain']['x'])]
    ax.set_ylim(ylim)
    ax.set_aspect(numpy.float(65.0*(xlim[1]-xlim[0])/(ylim[1]-ylim[0])/80.0))
    a = numpy.vstack([x,res[0,:]]).T
    b = numpy.vstack([x[::-1],res[2,:][::-1]]).T
    path = Path(numpy.vstack([a,b]))
    patch = patches.PathPatch(path, facecolor='blue', lw=0, alpha=0.35)
    ax.add_patch(patch)
    ax.plot(x,res[1,:],color='blue',alpha=1,lw=3)# ,zorder=0)
    if False and name=='d4.temp':
        gamma_matrix_sd = 0.375
        ax.plot(x,res[0,:]*(1.0-gamma_matrix_sd),color='blue',lw=1.5,alpha=1.0)
        ax.plot(x,res[2,:]*(1.0+gamma_matrix_sd),color='blue',lw=1.5,alpha=1.0)
    for p in fts['pop']:
        ax.plot(x,p,color='red',alpha=0.5,lw=1)# ,zorder=0)
    for mk in numpy.unique(typ):
        ax.plot(dat[typ==mk,0],dat[typ==mk,1],marker="o",markersize=15,markeredgecolor='black',markerfacecolor='white',linestyle='None')
        ax.plot(dat[typ==mk,0],dat[typ==mk,1],marker=r"$%s$" %(mk+1),markersize=12,markeredgecolor='black',markerfacecolor='black',linestyle='None')
    plt.xlabel(plts['xlabel'])
    plt.ylabel(plts['ylabel'])
    plt.savefig("%s/figures.%s.png" %(directory,name),bbox_inches="tight",transparent=True)
    plt.close()

############################################################

def plotPosterior(directory,model,prior,parmat=[]):
    """
    Plots prior and posterior distributions as a series a PNG files in the given directory
    """
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    import matplotlib.patches as patches
    from scipy.stats import norm
    # ---
    pplot(model, prior, parmat, {'directory': directory, 'name': 'deltaT', 'xlab': "Temperature difference ($^\circ$C)", 'ylab': "Probability", 'ylim': [0, 0.75], 'xlim': [-15, 15], 'scale': 1.0, 'prior_name': 'deltaT', 'par_name': 'new.deltaT'})
    pplot(model, prior, parmat, {'directory': directory, 'name': 'PP.ta.thr', 'xlab': "Temperature ($^\circ$C)", 'ylab': "Probability", 'ylim': [0, 0.75], 'xlim': [10, 30], 'scale': 1.0, 'prior_name': 'alpha.ta.thr', 'par_name': 'PP.ta.thr'})
    pplot(model, prior, parmat, {'directory': directory, 'name': 'PP.thr', 'xlab': "Photoperiod (hours)", 'ylab': "Probability", 'ylim': [], 'xlim': [0, 24], 'scale': 24.0, 'prior_name': '', 'par_name': 'PP.thr'})
    pplot(model, prior, parmat, {'directory': directory, 'name': 'PP.init', 'xlab': "Diapausing eggs (t=0)", 'ylab': "Probability", 'ylim': [], 'xlim': [], 'scale': 1.0, 'prior_name': '', 'par_name': 'PP.init'})
    pplot(model, prior, parmat, {'directory': directory, 'name': 'PP.strong', 'xlab': "p$_s$", 'ylab': "Probability", 'ylim': [], 'xlim': [], 'scale': 1.0, 'prior_name': '', 'par_name': 'PP.strong'})
    pplot(model, prior, parmat, {'directory': directory, 'name': 'PP.normal', 'xlab': "p$_n$", 'ylab': "Probability", 'ylim': [], 'xlim': [], 'scale': 1.0, 'prior_name': '', 'par_name': 'PP.normal'})
    fits = setPrior.getFits(model)
    for tsk in setPrior.tasks:
        if len(parmat)>0:
            print("Fitting %s..." %(tsk))
            plot_fit(parmat,setPrior.tasks[tsk],fits[tsk],tsk)
        print("Plotting %s..." %(tsk))
        plot_pos(prior,setPrior.tasks[tsk],setPrior.data[tsk],setPrior.datatype[tsk],directory,tsk,fits[tsk],tsk)
