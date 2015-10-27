import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

host = ['n4536','n4639','n3982','n3370','n3021','n1309','n4038','n5584','n4258'] 

marker = ['o','v','^','<','>','D','+','x','*']

alpha_eff = np.dtype([('Period',np.float32),('residual',np.float32),('error',np.float32),('alpha',np.float32),('name',np.str_,5)])

P,re,er,HP,galaxy = np.loadtxt('../output/effective_hyperparameters_cepheids.txt',unpack=True,usecols=[0,1,2,3,4],dtype=alpha_eff)

indexs = 0

indexf = 0

for index in range(len(P)):

    if galaxy[index] == host[0]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[0],mfc='k',fmt='',ls='None',label=host[0])

    elif HP[index] < 1. and HP[index] >= .5:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[0],mfc='b',fmt='',ls='None',label=host[0])

    elif HP[index] < .5 and HP[index] >= .3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[0],mfc='g',fmt='',ls='None',label=host[0])

    elif HP[index] < .3 and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[0],mfc='r',fmt='',ls='None',label=host[0])

    elif HP[index] < .1 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',mfc='y',marker=marker[0],fmt='',ls='None',label=host[0])


indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[1]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[1],mfc='k',fmt='',ls='None',label=host[1])

    elif HP[index] < 1. and HP[index] >= .5:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[1],mfc='b',fmt='',ls='None',label=host[1])

    elif HP[index] < .5 and HP[index] >= .3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[1],mfc='g',fmt='',ls='None',label=host[1])

    elif HP[index] < .3 and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[1],mfc='r',fmt='',ls='None',label=host[1])

    elif HP[index] < .1 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[1],mfc='y',fmt='',ls='None',label=host[1])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[2]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[2],mfc='k',fmt='',ls='None',label=host[2])

    elif HP[index] < 1. and HP[index] >= .5:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[2],mfc='b',fmt='',ls='None',label=host[2])

    elif HP[index] < .5 and HP[index] >= .3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[2],mfc='g',fmt='',ls='None',label=host[2])

    elif HP[index] < .3 and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[2],mfc='r',fmt='',ls='None',label=host[2])

    elif HP[index] < .1 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[2],mfc='y',fmt='',ls='None',label=host[2])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[3]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[3],mfc='k',fmt='',ls='None',label=host[3])

    elif HP[index] < 1. and HP[index] >= .5:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[3],mfc='b',fmt='',ls='None',label=host[3])

    elif HP[index] < .5 and HP[index] >= .3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[3],mfc='g',fmt='',ls='None',label=host[3])

    elif HP[index] < .3 and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[3],mfc='r',fmt='',ls='None',label=host[3])

    elif HP[index] < .1 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[3],mfc='y',fmt='',ls='None',label=host[3])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[4]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[4],mfc='k',fmt='',ls='None',label=host[4])

    elif HP[index] < 1. and HP[index] >= .5:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[4],mfc='b',fmt='',ls='None',label=host[4])

    elif HP[index] < .5 and HP[index] >= .3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[4],mfc='g',fmt='',ls='None',label=host[4])

    elif HP[index] < .3 and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[4],mfc='r',fmt='',ls='None',label=host[4])

    elif HP[index] < .1 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[4],mfc='y',fmt='',ls='None',label=host[4])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[5]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[5],mfc='k',fmt='',ls='None',label=host[5])

    elif HP[index] < 1. and HP[index] >= .5:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[5],mfc='b',fmt='',ls='None',label=host[5])

    elif HP[index] < .5 and HP[index] >= .3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[5],mfc='g',fmt='',ls='None',label=host[5])

    elif HP[index] < .3 and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[5],mfc='r',fmt='',ls='None',label=host[5])

    elif HP[index] < .1 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[5],mfc='y',fmt='',ls='None',label=host[5])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[6]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[6],mfc='k',fmt='',ls='None',label=host[6])

    elif HP[index] < 1. and HP[index] >= .5:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[6],mfc='b',fmt='',ls='None',label=host[6])

    elif HP[index] < .5 and HP[index] >= .3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[6],mfc='g',fmt='',ls='None',label=host[6])

    elif HP[index] < .3 and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[6],mfc='r',fmt='',ls='None',label=host[6])

    elif HP[index] < .1 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[6],mfc='y',fmt='',ls='None',label=host[6])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[7]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[7],mfc='k',fmt='',ls='None',label=host[7])

    elif HP[index] < 1. and HP[index] >= .5:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[7],mfc='b',fmt='',ls='None',label=host[7])

    elif HP[index] < .5 and HP[index] >= .3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[7],mfc='g',fmt='',ls='None',label=host[7])

    elif HP[index] < .3 and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[7],mfc='r',fmt='',ls='None',label=host[7])

    elif HP[index] < .1 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[7],mfc='y',fmt='',ls='None',label=host[7])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[8]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[8],mfc='k',fmt='',ls='None',label=host[8])

    elif HP[index] < 1. and HP[index] >= .5:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[8],mfc='b',fmt='',ls='None',label=host[8])

    elif HP[index] < .5 and HP[index] >= .3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[8],mfc='g',fmt='',ls='None',label=host[8])

    elif HP[index] < .3 and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[8],mfc='r',fmt='',ls='None',label=host[8])

    elif HP[index] < .1 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[8],mfc='y',fmt='',ls='None',label=host[8])

py.xscale('log')

py.yscale('linear')

#py.ylim(5.e-3,1.e1)

py.xlim(1.e0,5.e2)

py.xlabel('Period [days]')

py.ylabel('W [mag]')

#py.legend(loc=0,numpoints=1,ncol=4)

py.savefig('effective_HP_cepheids.pdf')

exit()


