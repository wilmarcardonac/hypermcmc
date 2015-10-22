import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

host = ['n4536','n4639','n3982','n3370','n3021','n1309','n4038','n5584','n4258'] 

P,HP,galaxy = np.loadtxt('../output/effective_hyperparameters.txt',unpack=True,usecols=[0,1,2])

for index in range(len(P)):

    if galaxy[index] == host[0]:

        py.plot(P[index],HP[index],'bo')

    elif galaxy[index] == host[1]:

        py.plot(P[index],HP[index],'gv')

    elif galaxy[index] == host[2]:

        py.plot(P[index],HP[index],'r^')

    elif galaxy[index] == host[3]:

        py.plot(P[index],HP[index],'k<')

    elif galaxy[index] == host[4]:

        py.plot(P[index],HP[index],'y>')

    elif galaxy[index] == host[5]:

        py.plot(P[index],HP[index],'mD')

    elif galaxy[index] == host[6]:

        py.plot(P[index],HP[index],'b+')

    elif galaxy[index] == host[7]:

        py.plot(P[index],HP[index],'gx')

    elif galaxy[index] == host[8]:

        py.plot(P[index],HP[index],'r*')


py.xscale('log')

py.yscale('log')

py.ylim(5.e-3,2.e0)

py.xlabel('Period [days]')

py.ylabel(r'$\alpha_{eff}$')

py.legend(["bo","gv","r^","k<","y>","mD","b+","gx","r*"],host,loc=4,numpoints=1)

py.savefig('effective_HP.pdf')

exit()
