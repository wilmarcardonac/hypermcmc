import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

def y(m,b,x):
    y = m*x - b
    return y

marker = ['o','s']#'o','v','^','<','>','D','+','x','*']

alpha_eff = np.dtype([('name',np.str_,5),('mu0_i',np.float32),('observed_mvi5av',np.float32),('error',np.float32),('theory',np.float32),('alpha',np.float32)])

snia,mui,mvi,er,th,HP = np.loadtxt('../output/chains/previous_runs/fit_44/effective_hyperparameters_SNIa.txt',unpack=True,usecols=[0,1,2,3,4,5],dtype=alpha_eff)

mean,error = np.loadtxt('../output/chains/previous_runs/fit_44/1Dstatistics_HP_R11_W.txt',unpack=True,usecols=[1,2],skiprows=31)

#A = mean[2] 

#B = error[2]

#mean[2] = mean[3]

#error[2] = error[3]

#mean[3] = A

#error[3] = B

#A = mean[6] 

#B = error[6]

#mean[6] = mean[7]

#error[6] = error[7]

#mean[7] = A

#error[7] = B

fig = py.figure()

xe =  np.linspace(13.,17.5)

ye = np.ones(len(xe))

for index in range(len(xe)):
    ye[index] = y(1.,th[0],xe[index])

py.plot(xe,ye,color='k')

#py.plot(mvi,mui,markersize='small',color='k')

py.xlim(13.,17.5)

py.ylim(-0.3,4.)

py.title('Relative distances from Cepheid variables and SNe Ia')

py.xlabel(r'$\mathrm{SN\, Ia} \,m_B^0+5a_B \,\mathrm{[mag]}$',fontsize=18)

py.ylabel(r'$\mathrm{Cepheid}\, (\mu_0 - \mu_{0,4258})\, \mathrm{[mag]}$',fontsize=18)

#exit()
#py.subplot(2,1,1)

#indexs = 0

#indexf = 0

#for index in range(len(P)):

#    if galaxy[index] == host[0]:

#        indexf = index

py.errorbar(mean[19],0.,xerr=error[19],yerr=None,ecolor='k',marker=marker[1],mfc='k',fmt='',ls='None')

for index in range(1,len(snia)):

    print snia[index], mvi[index], 'deltamu ', mean[index-1], 'error', error[index-1], 'HP ', HP[index]

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

py.savefig('../output/chains/previous_runs/fit_44/effective_HP_SNIa_R16.pdf')

exit()
py.plot(P,mw-re,markersize='small',color='k')

py.xscale('log')

py.xlim(1,2.e2)

py.title('LMC Cepheid variables')

#py.xlabel('Period [days]')

py.ylabel('W [mag]')

py.gca().invert_yaxis()

py.subplot(2,1,2)

indexs = 0

indexf = 0

for index in range(len(P)):

    if galaxy[index] == host[0]:

        indexf = index

for index in range(indexs,indexf+1):

    if HP[index] == 1.:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[0],mfc='k',fmt='',ls='None',label=host[0])

    elif HP[index] < 1. and HP[index] >= .1:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[0],mfc='b',fmt='',ls='None',label=host[0])

    elif HP[index] < .1 and HP[index] >= 1.e-2:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[0],mfc='g',fmt='',ls='None',label=host[0])

    elif HP[index] < 1.e-2 and HP[index] >= 1.e-3:

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',marker=marker[0],mfc='r',fmt='',ls='None',label=host[0])

    elif HP[index] < 1.e-3 :

        py.errorbar(P[index],re[index],yerr=er[index],xerr=None,ecolor='k',mfc='y',marker=marker[0],fmt='',ls='None',label=host[0])


py.xscale('log')

py.yscale('linear')

py.hlines(0.,1.,2.e2,color='k',linestyles='dotted')

#py.ylim(5.e-3,1.e1)

py.xlim(1.e0,2.e2)

py.xlabel('Period [days]')

py.ylabel('W residual [mag]')



#py.legend(loc=0,numpoints=1,ncol=4)


py.savefig('../output/chains/effective_HP_cepheids.pdf')

exit()





