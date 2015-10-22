import numpy as np
import pylab as py

P,HP = np.loadtxt('../output/effective_hyperparameters.txt',unpack=True,usecols=[0,1])

py.plot(P,HP,'bo')

py.xscale('log')

py.yscale('log')

py.ylim(5.e-3,2.e0)

py.xlabel('Period [days]')

py.ylabel(r'$\alpha_{eff}$')

py.savefig('effective_HP.pdf')

exit()
