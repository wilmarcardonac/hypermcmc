import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

host = ['n4536','n4639','n3982','n3370','n3021','n1309','n4038','n5584','n4258'] 

marker = ['bo','gv','r^','k<','y>','mD','b+','gx','r*']

alpha_eff = np.dtype([('N',np.int16),('alpha',np.float32),('name',np.str_,5)])

N,HP,galaxy = np.loadtxt('../output/effective_hyperparameters_hosts.txt',unpack=True,usecols=[0,1,2],dtype=alpha_eff)

for index in range(len(N)):

    if galaxy[index] == host[index]:

        py.plot(N[index],HP[index],marker[index],label=host[index])

py.xscale('linear')

py.yscale('log')

#py.ylim(5.e-3,1.e1)

py.xlim(0,10)

py.xlabel('Host')

py.ylabel(r'$\alpha_{eff}$')

py.legend(loc=0,numpoints=1,ncol=4)

py.savefig('effective_HP_hosts.pdf')

exit()


