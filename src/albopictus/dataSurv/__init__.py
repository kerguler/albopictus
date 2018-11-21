"""

dataset

    Parameters
    ----------

    filename: Name of the surveillance data file (format: CSV)

    drange: Surveillance duration ('custom','daily','weekly','biweekly','monthly','annual')

    dtype: Manually sets the data type of each column (default: None)

    Methods
    -------

    readData
    getPosMat
    getPosMatAll
    getDays

Examples

    surv.csv
    --------
    year,location,isoweek,samplings,mean,sd
    2008,LOC1,1,10,12.3,1.0
    2008,LOC1,2,8,13.4,1.1
    2008,LOC2,1,9,10.5,0.8
    2008,LOC2,2,10,10.2,0.85
    2008,LOC2,3,10,11.9,0.79
    --------

    from albopictus.dataSurv import dataset

    data = dataset("surv.csv","weekly")

    print(data.surv)
    
Warnings

    - Handling of dates have not been fully tested

"""

import numpy
from datetime import date, datetime, timedelta
from dateutil.relativedelta import relativedelta
from calendar import monthrange, isleap
from isoweek import Week

class dataset:
    def __init__(self,filename,drange,dtype=None):
        self.dranges = ['custom','daily','weekly','biweekly','monthly','annual']
        self.drange = drange
        if self.drange not in self.dranges:
            print("Data must be presented in any of these forms:")
            print(self.dranges)
            return
        self.filename = filename
        self.readData(dtype)
    #
    def readData(self,dtype=None):
        # Week 1 of a year starts with Monday
        data = numpy.genfromtxt(self.filename,names=True,delimiter=",",dtype=dtype)
        keys = data.dtype.fields.keys()
        loc = ['default']
        if 'location' in keys:
            # ---
            if dtype == None:
                dt = data.dtype.descr
                for i in numpy.arange(len(dt)):
                    if dt[i][0]=='location':
                        dt[i] = ('location',dt[i][1].replace('S','U').replace('s','u'))
                data = data.astype(numpy.dtype(dt))
            # ---
            loc = numpy.unique(data['location'])
        self.surv = {}
        for pr in loc:
            self.surv[pr] = {}
            for key in keys:
                self.surv[pr][key] = data[data['location']==pr][key] if 'location' in keys else data[key]
        for pr in self.surv.keys():
            srv = self.surv[pr]
            if 'day' in keys and 'month' in keys and 'year' in keys:
                self.surv[pr]['dates'] = numpy.array([date(int(srv['year'][i]),int(srv['month'][i]),int(srv['day'][i])) for i in range(srv['day'].shape[0])])
            elif 'month' in keys and 'year' in keys:
                self.surv[pr]['dates'] = numpy.array([date(int(srv['year'][i]),int(srv['month'][i]),1)+timedelta(days=monthrange(int(srv['year'][i]),int(srv['month'][i]))[1]-1) for i in range(srv['month'].shape[0])])
            # simple weeks
            elif 'week' in keys and 'year' in keys:
                self.surv[pr]['dates'] = numpy.array([datetime.strptime("%4d %02d 0" %(srv['year'][i],srv['week'][i]),'%Y %W %w').date() for i in range(srv['week'].shape[0])])
            # ISO weeks
            elif 'isoweek' in keys and 'year' in keys:
                self.surv[pr]['dates'] = numpy.array([Week(int(srv['year'][i]),int(srv['isoweek'][i])).sunday() for i in range(srv['isoweek'].shape[0])])
            # Excel weeks
            elif 'xweek' in keys and 'year' in keys:
                self.surv[pr]['dates'] = numpy.array([datetime.strptime("%4d %02d 0" %(srv['year'][i],srv['xweek'][i]-1),'%Y %W %w').date() for i in range(srv['xweek'].shape[0])])
            else:
                raise(ValueError, 'Cannot determine the dates of collections!')
            if 'duration' in keys:
                self.surv[pr]['duration'] = numpy.array(srv['duration'])
    #
    def getPosMat(self,dates):
        td = 0
        if self.drange == 'weekly':
            td = 6
        elif self.drange == 'biweekly':
            td = 13
        pos_mat = {}
        for pr in self.surv.keys():
            ret=[]
            for i in range(self.surv[pr]['dates'].shape[0]):
                if self.drange == 'monthly':
                    td = monthrange(self.surv[pr]['dates'][i].year,self.surv[pr]['dates'][i].month)[1]-1
                elif self.drange == 'annual':
                    td = 365+isleap(self.surv[pr]['dates'][i].year-1)-1
                elif self.drange == 'custom':
                    td = self.surv[pr]['duration'][i]-1
                b=numpy.repeat(0.0,dates.shape[0])
                b[(dates>=(self.surv[pr]['dates'][i]-timedelta(days=int(td)))) & (dates<=self.surv[pr]['dates'][i])]=1.0
                if (numpy.sum(b)>0):
                    ret.append(b.tolist())
            pos_mat[pr] = numpy.array(ret).T
        return pos_mat
    #
    def iterMonth(self,d,dr=1):
        return d + relativedelta(months=dr)
    #
    def iterDate(self,d,dr=1,custom=0):
        if self.drange == 'daily':
            d += dr*timedelta(days=1)
        elif self.drange == 'weekly':
            d += dr*timedelta(days=7)
        elif self.drange == 'biweekly':
            d += dr*timedelta(days=14)
        elif self.drange == 'monthly':
            d += relativedelta(months=dr)
        elif self.drange == 'annual':
            d += dr*timedelta(365+isleap(d.year+(dr if dr>0 else 0)))
        elif self.drange == 'custom':
            d += dr*timedelta(days=int(custom))
#            raise(ValueError, "iterDate is not defined for 'custom' surveillance!")
#            return None
        return d
    #
    def getPosMatAll(self,dates):
#        if self.drange == 'custom':
#            raise(ValueError, "getPosMatAll is not defined for 'custom' surveillance!")
#            return None
        pos_mat = {}
        pos_dates = {}
        for pr in self.surv.keys():
            ret = []
            dts = []
            ref = self.surv[pr]['dates'][0]
            custom = self.surv[pr]['duration'][0] if 'duration' in self.surv[pr] else 0
            if self.drange == 'monthly':
                ref = ref.replace(day=1)
            else:
                while ref>=dates[0]:
                    ref = self.iterDate(ref,dr=-1,custom=custom)
                ref = d = self.iterDate(ref,dr=1,custom=custom)
            while True:
                d = self.iterDate(ref,dr=1,custom=custom)
                if d>(dates[-1]+timedelta(days=1)):
                    break
                b=numpy.repeat(0.0,dates.shape[0])
                b[(dates>=ref) & (dates<d)]=1.0
                if (numpy.sum(b)>0):
                    ret.append(b.tolist())
                    dts.append(d)
                ref = d
            pos_mat[pr] = numpy.array(ret).T
            pos_dates[pr] = numpy.array(dts)
        return {'pos_mat':pos_mat,
                'pos_dates':pos_dates}
    #
    def getDays(self,dates):
        days = {}
        for pr in self.surv.keys():
            days[pr] = numpy.array([numpy.where(dates==i)[0][0] for i in self.surv[pr]['dates']])
        return days
