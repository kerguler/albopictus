import numpy

def calc_ss(pr,mn,ss,ss_inv=None):
    if (ss_inv is None):
        ss_inv = linalg.inv(ss)
    pr_mn = pr-mn
    if len(pr_mn)>1:
        return pr_mn.dot(ss_inv).dot(pr_mn.T)
    else:
        return (pr_mn*ss_inv*pr_mn)[0][0]

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

def getLonLat(times):
    """
    Python wrapper to obtain the longitude latitude grid for any climate dataset on http://vbd.cyi.ac.cy
    Parameters
    ----------
         times: climate dataset
                'local'        - E-OBS dataset with 1/4 deg resolution (2007-2013)
                'global'       - EMAC dataset with 1/2 deg resolution  (2000-2010)
                'future'       - EMAC dataset with 1/2 deg resolution  (2045-2055)
                'farfuture'    - EMAC dataset with 1/2 deg resolution  (2080-2100)
                'NASA_45_2018' - NASA ACCESS1-0 global downscaled NEX CMIP5 Climate Projection, r1i1p1, RCP4.5 (2017-2020)
                'NASA_85_2018' - NASA ACCESS1-0 global downscaled NEX CMIP5 Climate Projection, r1i1p1, RCP8.5 (2017-2020)
    """
    import json
    import requests
    #
    if not (times in ['local','global','future','farfuture','NASA_45_2018','NASA_85_2018']):
        print("Error: 'times' argument should be one of ['local','global','future','farfuture','NASA_45_2018','NASA_85_2018']")
        return None
    res = {
        'local':True,
        'global':False,
        'future':False,
        'farfuture':False,
        'NASA_45_2018':False,
        'NASA_85_2018':False
        }[times]
    #
    req = json.dumps({
        "lonlat":1,
        "lon":0,
        "lat":0,
        "times":times,
        "res":res,
        "sim":0,
        "version":""
        })
    r = requests.post("http://vbd.cyi.ac.cy/python", data={'request':req})
    if r.status_code!=200:
        print(r.status_code, r.reason)
        return None
    rt = json.loads(r.text)
    if len(rt['lonlat'])==0:
        print("Error: We do not have the climate dataset you requested in our database")
    return rt

def getClimate(lon,lat,times='global'):
    """
    Python wrapper to obtain climate variables from http://vbd.cyi.ac.cy
    Parameters
    ----------
         lon: longitude
         lat: latitude
         times: climate dataset (default: 'global')
                'local'        - E-OBS dataset with 1/4 deg resolution (2007-2013)
                'global'       - EMAC dataset with 1/2 deg resolution  (2000-2010)
                'future'       - EMAC dataset with 1/2 deg resolution  (2045-2055)
                'farfuture'    - EMAC dataset with 1/2 deg resolution  (2080-2100)
                'NASA_45_2018' - NASA ACCESS1-0 global downscaled NEX CMIP5 Climate Projection, r1i1p1, RCP4.5 (2017-2020)
                'NASA_85_2018' - NASA ACCESS1-0 global downscaled NEX CMIP5 Climate Projection, r1i1p1, RCP8.5 (2017-2020)
    """
    import json
    import requests
    #
    if not (times in ['local','global','future','farfuture','NASA_45_2018','NASA_85_2018']):
        print("Error: 'times' argument should be one of ['local','global','future','farfuture','NASA_45_2018','NASA_85_2018']")
        return None
    res = {
        'local':True,
        'global':False,
        'future':False,
        'farfuture':False,
        'NASA_45_2018':False,
        'NASA_85_2018':False
        }[times]
    #
    req = json.dumps({
        "lonlat":0,
        "lon":lon,
        "lat":lat,
        "times":times,
        "res":res,
        "sim":0,
        "version":'PDv1.1.1'
        })
    r = requests.post("http://vbd.cyi.ac.cy/python", data={'request':req})
    if r.status_code!=200:
        print(r.status_code, r.reason)
        return None
    rt = json.loads(r.text)
    if len(rt['mean_air_temp'])==0:
        print("Error: We do not have the climate variables you requested in our database")
    return rt
