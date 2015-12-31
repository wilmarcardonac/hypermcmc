import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

#host = ['LMC']

marker = ['o']#'o','v','^','<','>','D','+','x','*']

alpha_eff = np.dtype([('name',np.str_,5),('mu0_i',np.float32),('observed_mvi5av',np.float32),('error',np.float32),('theory',np.float32),('alpha',np.float32)])

snia,mui,mvi,er,th,HP = np.loadtxt('../output/chains/effective_hyperparameters_SNIa.txt',unpack=True,usecols=[0,1,2,3,4,5],dtype=alpha_eff)

mean,error = np.loadtxt('../output/chains/1Dstatistics_HP_R11_W.txt',unpack=True,usecols=[1,2],skiprows=17)

fig = py.figure()

#py.plot(mvi,mui,markersize='small',color='k')

py.xlim(13.4,17.5)

py.ylim(-0.3,3.5)

py.xlabel(r'SN Ia $m_v^0 + 5a_v$ [mag]')

py.ylabel(r'Cepheid $(\mu_0-\mu_{0,4258})$ [mag]')

py.errorbar(mean[8],0.,xerr=error[8],yerr=None,ecolor='k',marker=marker[0],mfc='k',fmt='',ls='None')

for index in range(1,len(snia)):

    if HP[index] == 1.:

        py.errorbar(mvi[index],mean[index-1],xerr=er[index]/np.sqrt(HP[index]),yerr=error[index-1],ecolor='k',marker=marker[0],mfc='k',fmt='',ls='None')

    elif HP[index] < 1. and HP[index] >= .1:

        py.errorbar(mvi[index],mean[index-1],xerr=er[index]/np.sqrt(HP[index]),yerr=error[index-1],ecolor='k',marker=marker[0],mfc='b',fmt='',ls='None')

    elif HP[index] < .1 and HP[index] >= 1.e-2:

        py.errorbar(mvi[index],mean[index-1],xerr=er[index]/np.sqrt(HP[index]),yerr=error[index-1],ecolor='k',marker=marker[0],mfc='g',fmt='',ls='None')

    elif HP[index] < 1.e-2 and HP[index] >= 1.e-3:

        py.errorbar(mvi[index],mean[index-1],xerr=er[index]/np.sqrt(HP[index]),yerr=error[index-1],ecolor='k',marker=marker[0],mfc='r',fmt='',ls='None')

    elif HP[index] < 1.e-3 :

        py.errorbar(mvi[index],mean[index-1],xerr=er[index]/np.sqrt(HP[index]),yerr=error[index-1],ecolor='k',mfc='y',marker=marker[0],fmt='',ls='None')

py.savefig('../output/chains/effective_HP_SNIa.pdf')

exit()






