import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

host = ['n4536','n4639','n3982','n3370','n3021','n1309','n4038','n5584','n4258'] 

marker = ['o','v','^','<','>','D','+','x','*']

alpha_eff = np.dtype([('Period',np.float32),('residual',np.float32),('error',np.float32),('alpha',np.float32),('name',np.str_,5)])

P,re,er,HP,galaxy = np.loadtxt('../output/chains/effective_hyperparameters_cepheids.txt',unpack=True,usecols=[0,1,2,3,4],dtype=alpha_eff)

fig = py.figure()

#fig.subplots_adjust(bottom=0.1, left=0.06, top=0.85, right=0.97,hspace=0.3)

#py.subplot(1,2,2)

py.xscale('log')

py.yscale('linear')

#py.ylim(5.e-3,1.e1)

py.xlim(1.e1,3.e2)

py.xlabel('Period [days]')

py.ylabel('F160W residual[mag]')

py.title('SNe hosts')

py.hlines(0.,1.,3.e2,color='k',linestyles='dotted')

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


indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[1]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[1],mfc='k',fmt='',ls='None',label=host[1])

    elif HP[index] < 1. and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[1],mfc='b',fmt='',ls='None',label=host[1])

    elif HP[index] < .1 and HP[index] >= 1.e-2:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[1],mfc='g',fmt='',ls='None',label=host[1])

    elif HP[index] < 1.e-2 and HP[index] >= 1.e-3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[1],mfc='r',fmt='',ls='None',label=host[1])

    elif HP[index] < 1.e-3 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[1],mfc='y',fmt='',ls='None',label=host[1])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[2]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[2],mfc='k',fmt='',ls='None',label=host[2])

    elif HP[index] < 1. and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[2],mfc='b',fmt='',ls='None',label=host[2])

    elif HP[index] < .1 and HP[index] >= 1.e-2:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[2],mfc='g',fmt='',ls='None',label=host[2])

    elif HP[index] < 1.e-2 and HP[index] >= 1.e-3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[2],mfc='r',fmt='',ls='None',label=host[2])

    elif HP[index] < 1.e-3 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[2],mfc='y',fmt='',ls='None',label=host[2])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[3]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[3],mfc='k',fmt='',ls='None',label=host[3])

    elif HP[index] < 1. and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[3],mfc='b',fmt='',ls='None',label=host[3])

    elif HP[index] < .1 and HP[index] >= 1.e-2:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[3],mfc='g',fmt='',ls='None',label=host[3])

    elif HP[index] < 1.e-2 and HP[index] >= 1.e-3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[3],mfc='r',fmt='',ls='None',label=host[3])

    elif HP[index] < 1.e-3 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[3],mfc='y',fmt='',ls='None',label=host[3])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[4]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[4],mfc='k',fmt='',ls='None',label=host[4])

    elif HP[index] < 1. and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[4],mfc='b',fmt='',ls='None',label=host[4])

    elif HP[index] < .1 and HP[index] >= 1.e-2:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[4],mfc='g',fmt='',ls='None',label=host[4])

    elif HP[index] < 1.e-2 and HP[index] >= 1.e-3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[4],mfc='r',fmt='',ls='None',label=host[4])

    elif HP[index] < 1.e-3 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[4],mfc='y',fmt='',ls='None',label=host[4])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[5]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[5],mfc='k',fmt='',ls='None',label=host[5])

    elif HP[index] < 1. and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[5],mfc='b',fmt='',ls='None',label=host[5])

    elif HP[index] < .1 and HP[index] >= 1.e-2:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[5],mfc='g',fmt='',ls='None',label=host[5])

    elif HP[index] < 1.e-2 and HP[index] >= 1.e-3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[5],mfc='r',fmt='',ls='None',label=host[5])

    elif HP[index] < 1.e-3 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[5],mfc='y',fmt='',ls='None',label=host[5])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[6]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[6],mfc='k',fmt='',ls='None',label=host[6])

    elif HP[index] < 1. and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[6],mfc='b',fmt='',ls='None',label=host[6])

    elif HP[index] < .1 and HP[index] >= 1.e-2:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[6],mfc='g',fmt='',ls='None',label=host[6])

    elif HP[index] < 1.e-2 and HP[index] >= 1.e-3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[6],mfc='r',fmt='',ls='None',label=host[6])

    elif HP[index] < 1.e-3 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[6],mfc='y',fmt='',ls='None',label=host[6])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[7]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[7],mfc='k',fmt='',ls='None',label=host[7])

    elif HP[index] < 1. and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[7],mfc='b',fmt='',ls='None',label=host[7])

    elif HP[index] < .1 and HP[index] >= 1.e-2:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[7],mfc='g',fmt='',ls='None',label=host[7])

    elif HP[index] < 1.e-2 and HP[index] >= 1.e-3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[7],mfc='r',fmt='',ls='None',label=host[7])

    elif HP[index] < 1.e-3 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[7],mfc='y',fmt='',ls='None',label=host[7])

indexs = indexf

py.savefig('../output/chains/effective_HP_F160W_cepheids_SNe.pdf')

py.close(fig)

fig = py.figure()

#py.subplot(1,2,1)

py.title('NGC4258')

py.xscale('log')

py.yscale('linear')

#py.ylim(5.e-3,1.e1)

py.xlim(2.e0,3.e2)

py.xlabel('Period [days]')

py.ylabel('F160W residual [mag]')

py.hlines(0.,2.,3.e2,color='k',linestyles='dotted')

for index in range(len(P)):

    if galaxy[index] == host[8]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[8],mfc='k',fmt='',ls='None',label=host[8])

    elif HP[index] < 1. and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[8],mfc='b',fmt='',ls='None',label=host[8])

    elif HP[index] < .1 and HP[index] >= 1.e-2:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[8],mfc='g',fmt='',ls='None',label=host[8])

    elif HP[index] < 1.e-2 and HP[index] >= 1.e-3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[8],mfc='r',fmt='',ls='None',label=host[8])

    elif HP[index] < 1.e-3 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[8],mfc='y',fmt='',ls='None',label=host[8])


#py.legend(loc=0,numpoints=1,ncol=4)

py.savefig('../output/chains/effective_HP_F160W_cepheids_NGC4258.pdf')

py.close(fig)

exit()


