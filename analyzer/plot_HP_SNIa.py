import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

#host = ['LMC']

marker = ['o']#'o','v','^','<','>','D','+','x','*']

alpha_eff = np.dtype([('name',np.str_,5),('mu0_i',np.float32),('observed_mvi5av',np.float32),('error',np.float32),('theory',np.float32),('alpha',np.float32)])

snia,mui,mvi,er,th,HP = np.loadtxt('../output/chains/effective_hyperparameters_SNIa.txt',unpack=True,usecols=[0,1,2,3,4,5],dtype=alpha_eff)

fig = py.figure()

#fig.subplots_adjust(bottom=0.1, left=0.06, top=0.85, right=0.97,hspace=0.3)

py.plot(th,mui,markersize='small',color='k')

#py.xscale('log')

py.xlim(13.4,17.5)

py.ylim(29.,33.)

#py.title('LMC Cepheid variables')

#py.xlabel('Period [days]')

#py.ylabel('W [mag]')

#exit()
#py.subplot(2,1,1)

#indexs = 0

#indexf = 0

#for index in range(len(P)):

#    if galaxy[index] == host[0]:

#        indexf = index

for index in range(len(snia)):

    if HP[index] == 1.:

        py.errorbar(mvi[index],mui[index],xerr=er[index]/np.sqrt(HP[index]),yerr=None,ecolor='k',marker=marker[0],mfc='k',fmt='',ls='None')

    elif HP[index] < 1. and HP[index] >= .1:

        py.errorbar(mvi[index],mui[index],xerr=er[index]/np.sqrt(HP[index]),yerr=None,ecolor='k',marker=marker[0],mfc='b',fmt='',ls='None')

    elif HP[index] < .1 and HP[index] >= 1.e-2:

        py.errorbar(mvi[index],mui[index],xerr=er[index]/np.sqrt(HP[index]),yerr=None,ecolor='k',marker=marker[0],mfc='g',fmt='',ls='None')

    elif HP[index] < 1.e-2 and HP[index] >= 1.e-3:

        py.errorbar(mvi[index],mui[index],xerr=er[index]/np.sqrt(HP[index]),yerr=None,ecolor='k',marker=marker[0],mfc='r',fmt='',ls='None')

    elif HP[index] < 1.e-3 :

        py.errorbar(mvi[index],mui[index],xerr=er[index]/np.sqrt(HP[index]),yerr=None,ecolor='k',mfc='y',marker=marker[0],fmt='',ls='None')

#py.scatter(P,mw)

py.savefig('../output/chains/effective_HP_SNIa.pdf')

exit()
py.plot(P,mw-re,markersize='small',color='k')

py.xscale('log')

py.xlim(1,2.e2)

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

py.hlines(0.,1.,2.e2,color='k',linestyles='dotted')

#py.ylim(5.e-3,1.e1)

py.xlim(1.e0,2.e2)

py.xlabel('Period [days]')

py.ylabel('W residual [mag]')



#py.legend(loc=0,numpoints=1,ncol=4)


py.savefig('../output/chains/effective_HP_cepheids.pdf')

exit()





