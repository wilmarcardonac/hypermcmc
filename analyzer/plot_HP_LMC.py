import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py
import math

def Mwt(A,b,Pe):
    Mwt = A + b*(math.log10(Pe)-1.)
    return Mwt

host = ['LMC']

Pe = np.linspace(2.,300.,num=50.)

marker = ['o']#'o','v','^','<','>','D','+','x','*']

alpha_eff = np.dtype([('Period',np.float32),('observed_w',np.float32),('residual',np.float32),('error',np.float32),('alpha',np.float32),('name',np.str_,5)])

P,mw,re,er,HP,galaxy = np.loadtxt('../output/chains/previous_runs/fit_c/effective_hyperparameters_cepheids.txt',unpack=True,usecols=[0,1,2,3,4,5],dtype=alpha_eff)

si = np.loadtxt('../output/chains/previous_runs/fit_c/means.txt')

bi = np.loadtxt('../output/chains/previous_runs/fit_c/bestfit.txt')

bd = np.loadtxt('../output/chains/previous_runs/fit_d/bestfit.txt')

be =  np.ones(len(Pe))

bd2 =  np.ones(len(Pe))

for index in range(len(Pe)):
    be[index] = Mwt(bi[0],bi[1],Pe[index])
    bd2[index] = Mwt(bd[0],bd[1],Pe[index])

er = np.sqrt(er**2 + (10.**(si[2]))**2)

fig = py.figure()

#fig.subplots_adjust(bottom=0.1, left=0.06, top=0.85, right=0.97,hspace=0.3)

py.subplot(2,1,1)

indexs = 0

indexf = 0

for index in range(len(P)):

    if galaxy[index] == host[0]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],mw[index],yerr=er[index]/np.sqrt(HP[index]),xerr=None,ecolor='k',marker=marker[0],mfc='k',fmt='',ls='None',label=host[0])

    elif HP[index] < 1. and HP[index] >= .1:

        py.errorbar(P[index],mw[index],yerr=er[index]/np.sqrt(HP[index]),xerr=None,ecolor='k',marker=marker[0],mfc='b',fmt='',ls='None',label=host[0])

    elif HP[index] < .1 and HP[index] >= 1.e-2:

        py.errorbar(P[index],mw[index],yerr=er[index]/np.sqrt(HP[index]),xerr=None,ecolor='k',marker=marker[0],mfc='g',fmt='',ls='None',label=host[0])

    elif HP[index] < 1.e-2 and HP[index] >= 1.e-3:

        py.errorbar(P[index],mw[index],yerr=er[index]/np.sqrt(HP[index]),xerr=None,ecolor='k',marker=marker[0],mfc='r',fmt='',ls='None',label=host[0])

    elif HP[index] < 1.e-3 :

        py.errorbar(P[index],mw[index],yerr=er[index]/np.sqrt(HP[index]),xerr=None,ecolor='k',mfc='y',marker=marker[0],fmt='',ls='None',label=host[0])

#py.scatter(P,mw)

py.plot(Pe,be,markersize='small',color='k')

py.plot(Pe,bd2,color='r',ls='dotted',linewidth=2)

py.xscale('log')

py.xlim(2,2.e2)

py.title('LMC Cepheid variables')

#py.xlabel('Period [days]')

py.ylabel('W [mag]')

py.gca().invert_yaxis()

py.subplot(2,1,2)

indexs = 0

indexf = 0

for index in range(len(P)):

    if galaxy[index] == host[0]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[0],mfc='k',fmt='',ls='None',label=host[0])

    elif HP[index] < 1. and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[0],mfc='b',fmt='',ls='None',label=host[0])

    elif HP[index] < .1 and HP[index] >= 1.e-2:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[0],mfc='g',fmt='',ls='None',label=host[0])

    elif HP[index] < 1.e-2 and HP[index] >= 1.e-3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[0],mfc='r',fmt='',ls='None',label=host[0])

    elif HP[index] < 1.e-3 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',mfc='y',marker=marker[0],fmt='',ls='None',label=host[0])


py.xscale('log')

py.yscale('linear')

py.hlines(0.,1.,2.e2,color='k',linestyles='dotted',linewidth=2)

py.plot(Pe,bd2-be,color='r',ls='dotted',linewidth=2)

#py.ylim(5.e-3,1.e1)

py.xlim(2.,2.e2)

py.xlabel('Period [days]')

py.ylabel('W residual [mag]')



#py.legend(loc=0,numpoints=1,ncol=4)


py.savefig('../output/chains/previous_runs/effective_HP_cepheids_LMC.pdf')

exit()





