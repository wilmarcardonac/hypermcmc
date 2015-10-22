import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

host = ['n4536','n4639','n3982','n3370','n3021','n1309','n4038','n5584','n4258'] 

marker = ['bo','gv','r^','k<','y>','mD','b+','gx','r*']

alpha_eff = np.dtype([('Period',np.float32),('alpha',np.float32),('name',np.str_,5)])

P,HP,galaxy = np.loadtxt('../output/effective_hyperparameters.txt',unpack=True,usecols=[0,1,2],dtype=alpha_eff)

indexs = 0

indexf = 0

for index in range(len(P)):

    if galaxy[index] == host[0]:

        indexf = index

py.plot(P[indexs:indexf+1],HP[indexs:indexf+1],marker[0],label=host[0])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[1]:

        indexf = index

py.plot(P[indexs:indexf+1],HP[indexs:indexf+1],marker[1],label=host[1])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[2]:

        indexf = index

py.plot(P[indexs:indexf+1],HP[indexs:indexf+1],marker[2],label=host[2])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[3]:

        indexf = index

py.plot(P[indexs:indexf+1],HP[indexs:indexf+1],marker[3],label=host[3])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[4]:

        indexf = index

py.plot(P[indexs:indexf+1],HP[indexs:indexf+1],marker[4],label=host[4])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[5]:

        indexf = index

py.plot(P[indexs:indexf+1],HP[indexs:indexf+1],marker[5],label=host[5])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[6]:

        indexf = index

py.plot(P[indexs:indexf+1],HP[indexs:indexf+1],marker[6],label=host[6])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[7]:

        indexf = index

py.plot(P[indexs:indexf+1],HP[indexs:indexf+1],marker[7],label=host[7])

indexs = indexf

for index in range(len(P)):

    if galaxy[index] == host[8]:

        indexf = index

py.plot(P[indexs:indexf+1],HP[indexs:indexf+1],marker[8],label=host[8])

py.xscale('log')

py.yscale('log')

py.ylim(5.e-3,1.e1)

py.xlim(1.e0,5.e2)

py.xlabel('Period [days]')

py.ylabel(r'$\alpha_{eff}$')

py.legend(loc=0,numpoints=1,ncol=4)

py.savefig('effective_HP.pdf')

exit()


