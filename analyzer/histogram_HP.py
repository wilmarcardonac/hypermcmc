import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py
from matplotlib.colors import LogNorm

host = ['m101 ','n1015','n1309','n1365','n1448','n2442','n3021','n3370','n3447','n3972','n3982','n4038','n4424','n4536','n4639','n5584','n5917','n7250','u9391','n4258'] # HOST GALAXIES IN SAME ORDER LISTED IN TABLE 5 OF R16

marker = ['o','v','^','<','>','D','+','x','*']

alpha_eff = np.dtype([('Period',np.float32),('residual',np.float32),('error',np.float32),('alpha',np.float32),('name',np.str_,5)])

bins = [0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.]

P1,re1,er1,HP1,galaxy1 = np.loadtxt('../output/chains/effective_hyperparameters_cepheids_R11.txt',unpack=True,usecols=[0,2,3,4,5],dtype=alpha_eff)

P2,re2,er2,HP2,galaxy2 = np.loadtxt('../output/chains/effective_hyperparameters_cepheids_LMC.txt',unpack=True,usecols=[0,2,3,4,5],dtype=alpha_eff)

P3,re3,er3,HP3,galaxy3 = np.loadtxt('../output/chains/effective_hyperparameters_cepheids_MW.txt',unpack=True,usecols=[0,2,3,4,5],dtype=alpha_eff)

P4,re4,er4,HP4,galaxy4 = np.loadtxt('../output/chains/effective_hyperparameters_cepheids_M31.txt',unpack=True,usecols=[0,2,3,4,5],dtype=alpha_eff)

py.hist(HP1,bins=bins,alpha=0.2,label='R16')

py.hist(HP2,bins=bins,alpha=0.4,label='LMC')

py.hist(HP3,bins=bins,alpha=0.5,label='MW')

py.hist(HP4,bins=bins,alpha=0.6,label='M31')

py.legend(loc=0)

py.savefig('../output/chains/effective_HP_histogram_R16.pdf')

py.close()

#py.hist2d(P1,HP1,bins=[[0,60,205],bins],norm=LogNorm())

#py.colorbar()

#py.savefig('../output/chains/effective_HP_2D_histogram_R16.pdf')

print 'TOTAL NUMBER OF CEPHEIDS IN R16 SAMPLE', len(HP1)

counter = 0

for index in range(len(HP1)):
    if HP1[index] == 1.:
        counter = counter + 1

print 'NUMBER OF CEPHEIDS WITH HP EQUAL TO ONE IN R16 SAMPLE :', counter, 'WHICH IS % ',counter/float(len(HP1))*100

counter = 0

for index in range(len(HP1)):
    if (HP1[index] < 1. and HP1[index]>=0.1):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN R16 SAMPLE 0.1 >= HP < 1:', counter, 'WHICH IS % ',counter/float(len(HP1))*100

counter = 0

for index in range(len(HP1)):
    if (HP1[index] < 0.1 and HP1[index]>=0.01):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN R16 SAMPLE 0.01 >= HP < 0.1:', counter, 'WHICH IS % ',counter/float(len(HP1))*100

counter = 0

for index in range(len(HP1)):
    if (HP1[index] < 0.01 and HP1[index]>=0.001):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN R16 SAMPLE 0.001 >= HP < 0.01:', counter, 'WHICH IS % ',counter/float(len(HP1))*100

counter = 0

for index in range(len(HP1)):
    if (HP1[index] < 0.001):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN R16 SAMPLE HP < 0.001:', counter, 'WHICH IS % ',counter/float(len(HP1))*100

print 'TOTAL NUMBER OF CEPHEIDS IN LMC SAMPLE', len(HP2)

counter = 0

for index in range(len(HP2)):
    if HP2[index] == 1.:
        counter = counter + 1

print 'NUMBER OF CEPHEIDS WITH HP EQUAL TO ONE IN LMC SAMPLE :', counter, 'WHICH IS % ',counter/float(len(HP2))*100

counter = 0

for index in range(len(HP2)):
    if (HP2[index] < 1. and HP2[index]>=0.1):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN LMC 0.1>= HP <1. :', counter, 'WHICH IS % ',counter/float(len(HP2))*100

counter = 0

for index in range(len(HP2)):
    if (HP2[index] < 0.1 and HP2[index]>=0.01):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN LMC SAMPLE 0.01 >= HP < 0.1:', counter, 'WHICH IS % ',counter/float(len(HP2))*100

counter = 0

for index in range(len(HP2)):
    if (HP2[index] < 0.01 and HP2[index]>=0.001):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN LMC SAMPLE 0.001 >= HP < 0.01:', counter, 'WHICH IS % ',counter/float(len(HP2))*100

counter = 0

for index in range(len(HP2)):
    if (HP2[index] < 0.001):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN LMC SAMPLE HP < 0.001:', counter, 'WHICH IS % ',counter/float(len(HP2))*100

print 'TOTAL NUMBER OF CEPHEIDS IN MW SAMPLE', len(HP3)

counter = 0

for index in range(len(HP3)):
    if HP3[index] == 1:
        counter = counter + 1

print 'NUMBER OF CEPHEIDS WITH HP EQUAL TO ONE IN MW :', counter, 'WHICH IS % ',counter/float(len(HP3))*100

counter = 0

for index in range(len(HP3)):
    if (HP3[index] < 1. and HP3[index]>=0.1):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN MW 0.1>= HP <1. :', counter, 'WHICH IS % ',counter/float(len(HP3))*100

counter = 0

for index in range(len(HP3)):
    if (HP3[index] < 0.1 and HP3[index]>=0.01):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN MW SAMPLE 0.01 >= HP < 0.1:', counter

counter = 0

for index in range(len(HP3)):
    if (HP3[index] < 0.01 and HP3[index]>=0.001):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN MW SAMPLE 0.001 >= HP < 0.01:', counter

counter = 0

for index in range(len(HP3)):
    if (HP3[index] < 0.001):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN MW SAMPLE HP < 0.001:', counter

print 'TOTAL NUMBER OF CEPHEIDS IN M31 SAMPLE', len(HP4)

counter = 0

for index in range(len(HP4)):
    if HP4[index] == 1:
        counter = counter + 1

print 'NUMBER OF CEPHEIDS WITH HP EQUAL TO ONE IN M31 :', counter, 'WHICH IS % ',counter/float(len(HP4))*100

counter = 0

for index in range(len(HP4)):
    if (HP4[index] < 1. and HP4[index]>=0.1):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN M31 0.1>= HP <1. :', counter, 'WHICH IS % ',counter/float(len(HP4))*100

counter = 0

for index in range(len(HP4)):
    if (HP4[index] < 0.1 and HP4[index]>=0.01):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN M31 SAMPLE 0.01 >= HP < 0.1:', counter

counter = 0

for index in range(len(HP4)):
    if (HP4[index] < 0.01 and HP4[index]>=0.001):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN M31 SAMPLE 0.001 >= HP < 0.01:', counter

counter = 0

for index in range(len(HP4)):
    if (HP4[index] < 0.001):
        counter = counter + 1

print 'NUMBER OF OUTLIERS IN M31 SAMPLE HP < 0.001:', counter

exit()

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

py.savefig('../output/chains/effective_HP_cepheids.pdf')

exit()


