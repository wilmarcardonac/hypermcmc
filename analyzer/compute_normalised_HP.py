import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

alpha_eff = np.dtype([('Period',np.float32),('observed_w',np.float32),('residual',np.float32),('error',np.float32),('alpha',np.float32),('name',np.str_,5),('Z_w',np.float32)])

P,mw,re,er,HP,galaxy,zw = np.loadtxt('../output/chains/effective_hyperparameters_cepheids_R11.txt',unpack=True,usecols=[0,1,2,3,4,5,6],dtype=alpha_eff)

P,mw,re,er,HP1,galaxy,zw = np.loadtxt('../output/chains/effective_hyperparameters_cepheids_LMC.txt',unpack=True,usecols=[0,1,2,3,4,5,6],dtype=alpha_eff)

P,mw,re,er,HP2,galaxy,zw = np.loadtxt('../output/chains/effective_hyperparameters_cepheids_MW.txt',unpack=True,usecols=[0,1,2,3,4,5,6],dtype=alpha_eff)

P,mw,re,er,HP3,galaxy,zw = np.loadtxt('../output/chains/effective_hyperparameters_cepheids_M31.txt',unpack=True,usecols=[0,1,2,3,4,5,6],dtype=alpha_eff)

alpha_eff1 = np.dtype([('name',np.str_,5),('mu0_i',np.float32),('observed_mvi5av',np.float32),('error',np.float32),('theory',np.float32),('alpha',np.float32)])

snia,mui,mvi,er3,th,HPsnia = np.loadtxt('../output/chains/effective_hyperparameters_SNIa.txt',unpack=True,usecols=[0,1,2,3,4,5],dtype=alpha_eff1)

#print 'HP R11 Cepheids ', np.sum(HP[:])/len(HP[:])

#print 'HP R11 LMC Cepheids ', np.sum(HP1[:])/len(HP1[:])

print 'HP Cepheids ', (np.sum(HP[:]) + np.sum(HP1[:]) + np.sum(HP2[:]) + np.sum(HP3[:]) )/(len(HP[:]) + len(HP1[:]) + len(HP2[:]) + len(HP3[:]) )

print 'HP R11 SN Ia ', np.sum(HPsnia[1:])/len(HPsnia[1:])

exit()


