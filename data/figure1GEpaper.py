import numpy as np
import pylab as py

# A new data type 

fulldata = np.dtype([('Name',np.str_,20),('Period',np.float32),('H',np.float32),('Sigma_m',np.float32),('V',np.float32),('I',np.float32)])

N,P,H,Sm,V,I = np.loadtxt('data.txt',unpack=True,usecols=[0,1,2,3,4,5],dtype=fulldata)

mw = H - 0.41*(V-I)

py.scatter(P,mw)

py.xscale('log')

py.xlim(2,200)

py.xlabel('Period [days]')

py.ylabel('W [mag]')

py.gca().invert_yaxis()

py.show()

exit()
