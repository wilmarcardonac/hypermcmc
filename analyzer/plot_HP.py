import numpy as np
import pylab as py

P,HP = np.loadtxt('../output/effective_hyperparameters.txt',unpack=True,usecols=[0,1])

py.plot(P,HP)

py.xscale('log')

py.yscale('linear')

py.ylim(0.,1.)

py.xlabel('Period [days]')

py.ylabel(r'$\alpha_{eff}$')

py.savefig('effective_HP.pdf')

exit()
