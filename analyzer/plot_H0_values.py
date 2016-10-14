import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

#host = ['LMC']

# Indirect determinations 

wmap09  = [70.0,2.2]

planck13 = [67.9,1.5]

planck15 = [67.81,0.92]

# Riess 2011

ngc4258_R11 = [74.8,3.0]

lmc_R11 = [71.3,3.8]

mw_R11 = [75.7,2.6]

ngc4258_lmc_R11 = [72.3,2.3]

ngc4258_mw_R11 = [74.5,2.3]

lmc_mw_R11 = [74.4,2.4]

ngc4258_lmc_mw_R11 = [73.8,2.4]

# Efstathiou 2015

ngc4258_E15 = [70.6,3.3]

lmc_E15 = [73.4,3.1]

mw_E15 = [75.3,3.5]

ngc4258_lmc_E15 = [71.8,2.6]

ngc4258_mw_E15 = [72.2,2.8]

lmc_mw_E15 = [73.9,2.7]

ngc4258_lmc_mw_E15 = [72.5,2.5]

# Riess 2016

ngc4258_R16 = [72.25,2.51]

lmc_R16 = [72.04,2.67]

mw_R16 = [76.18,2.37]

#ngc4258_lmc_ = [,]

ngc4258_mw_R16 = [74.04,1.93]

#lmc_mw_ = [,]

ngc4258_lmc_mw_R16 = [73.24,1.74]

# HP

ngc4258_HP = [72.1,4.2]

lmc_HP = [74.4,4.0]

mw_HP = [78.1,4.4]

ngc4258_lmc_HP = [71.1,4.0]

ngc4258_mw_HP = [76.9,4.0]

lmc_mw_HP = [76.1,3.8]

ngc4258_lmc_mw_HP_M1a = [75.0,3.9]

ngc4258_lmc_mw_HP_M1ag = [73.2,2.5]

ngc4258_lmc_mw_HP_M1ah = [74.1,3.7]

ngc4258_lmc_mw_HP_M1aj = [72.4,2.2]

marker = ['o','D','v','>','<','*','^']

fig = py.figure()

#py.errorbar(70.5,2009,yerr=None,xerr=1.3,ecolor='k',marker=marker[0],mfc='k',fmt='',ls='None',label='WMAP09')

#py.errorbar(74.2,2009,yerr=None,xerr=3.6,ecolor='b',marker=marker[1],mfc='b',fmt='',ls='None',label='Riess09')

#py.errorbar(73.8,2011,yerr=None,xerr=2.4,ecolor='r',marker=marker[1],mfc='r',fmt='',ls='None',label='Riess11')

#py.errorbar(67.9,2013,yerr=None,xerr=1.5,ecolor='k',marker=marker[0],mfc='k',fmt='',ls='None',label='PLANCK13')

#py.errorbar(70.6,2015,yerr=None,xerr=3.3,ecolor='y',mfc='y',marker=marker[1],fmt='',ls='None',label='Efstathiou15')

#py.errorbar(67.81,2015,yerr=None,xerr=0.92,ecolor='k',marker=marker[0],mfc='k',fmt='',ls='None',label='PLANCK15')

#py.errorbar(73.0,2016,yerr=None,xerr=1.8,ecolor='g',marker=marker[1],mfc='g',fmt='',ls='None',label='Riess16')

#py.scatter(P,mw)

#py.plot(P,mw-re,markersize='small',color='k')

#py.xscale('log')

#py.xlim(60.,80.)

#py.ylim(2007,2017)

#py.title(' Cepheid variables')

#py.xlabel(r'$H_0$ [km/s/Mpc]')

#py.ylabel('Year')

#py.legend(loc=0,numpoints=1,ncol=1)

#py.gca().invert_yaxis()

#py.subplot(2,1,2)

#py.savefig('./H0_values.pdf')

# NGC4258 + LMC + MW

fig = py.figure()

# Indirect determinations

#py.errorbar(wmap09[0],wmap09[1]/wmap09[0]*100.,yerr=None,xerr=wmap09[1],ecolor='blue',marker=marker[0],mfc='blue',fmt='',ls='None',label='WMAP 2009')

#py.errorbar(planck13[0],planck13[1]/planck13[0]*100.,yerr=None,xerr=planck13[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK13')

py.errorbar(planck15[0],1,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK 2015',elinewidth=3,ms=10)
py.annotate('Planck 2015 (1.4%)', xy=(66,3), xycoords='data', size=14, color='black')

# Riess 2011

py.errorbar(ngc4258_lmc_mw_R11[0],18,yerr=None,xerr=ngc4258_lmc_mw_R11[1],ecolor='#f91fb4',marker=marker[0],mfc='#f91fb4',fmt='',ls='None',label='Riess et al. 2011',elinewidth=3,ms=10)
py.annotate('Riess et al. 2011 (3.3%)', xy=(66,20), xycoords='data', size=14, color='#f91fb4')

# Efstathiou 2015 

py.errorbar(ngc4258_lmc_mw_E15[0],25,yerr=None,xerr=ngc4258_lmc_mw_E15[1],ecolor='dodgerblue',marker=marker[0],mfc='dodgerblue',fmt='',ls='None',label='Efstathiou 2014',elinewidth=3,ms=10)
py.annotate('Efstathiou 2014 (3 anchors) (3.4%)', xy=(66,27), xycoords='data', size=14, color='dodgerblue')

# Riess 2016

#py.errorbar(ngc4258_lmc_mw_R16[0],ngc4258_lmc_mw_R16[1]/ngc4258_lmc_mw_R16[0]*100.,yerr=None,xerr=ngc4258_lmc_mw_R16[1],ecolor='green',marker=marker[2],mfc='green',fmt='',ls='None',label='Riess et al. 2016')

# HP

py.errorbar(ngc4258_lmc_mw_HP_M1a[0],46,yerr=None,xerr=ngc4258_lmc_mw_HP_M1a[1],ecolor='#0c54c0',marker=marker[0],mfc='#0c54c0',fmt='',ls='None',label=r'HPs in Cepheids, SN Ia, and distance modulus',elinewidth=3,ms=10)
py.annotate('Cardona et al. 2016 (5.2%) [R11 data with 722 HP]', xy=(66,48), xycoords='data', size=14, color='#0c54c0')

print 'HPs in Cepheids, SN Ia, and distance modulus', ngc4258_lmc_mw_HP_M1a[1]/ngc4258_lmc_mw_HP_M1a[0]*100.

py.errorbar(ngc4258_lmc_mw_HP_M1ag[0],32,yerr=None,xerr=ngc4258_lmc_mw_HP_M1ag[1],ecolor='#0c54c0',marker=marker[0],mfc='#0c54c0',fmt='',ls='None',label=r'HPs in Cepheids and distance modulus',elinewidth=3,ms=10)
py.annotate('Cardona et al. 2016 (3.4%) [R11 data with 714 HP]', xy=(66,34), xycoords='data', size=14, color='#0c54c0')

print 'HPs in Cepheids and distance modulus', ngc4258_lmc_mw_HP_M1ag[1]/ngc4258_lmc_mw_HP_M1ag[0]*100.

py.errorbar(ngc4258_lmc_mw_HP_M1ah[0],39,yerr=None,xerr=ngc4258_lmc_mw_HP_M1ah[1],ecolor='#0c54c0',marker=marker[0],mfc='#0c54c0',fmt='',ls='None',label=r'HPs in Cepheids and SN Ia',elinewidth=3,ms=10)
py.annotate('Cardona et al. 2016 (5.0%) [R11 data with 720 HP]', xy=(66,41), xycoords='data', size=14, color='#0c54c0')

print 'HPs in Cepheids and SN Ia', ngc4258_lmc_mw_HP_M1ah[1]/ngc4258_lmc_mw_HP_M1ah[0]*100.

py.errorbar(ngc4258_lmc_mw_HP_M1aj[0],11,yerr=None,xerr=ngc4258_lmc_mw_HP_M1aj[1],ecolor='#0c54c0',marker=marker[0],mfc='#0c54c0',fmt='',ls='None',label=r'HPs in Cepheids',elinewidth=3,ms=10)
py.annotate('Cardona et al. 2016 (3.0%) [R11 data with 712 HP]', xy=(66,13), xycoords='data', size=14, color='#0c54c0')

print 'HPs in Cepheids', ngc4258_lmc_mw_HP_M1aj[1]/ngc4258_lmc_mw_HP_M1aj[0]*100.

py.xlim(65.,81.)

py.ylim(-2.,55.)

#py.title('NGC 4258 + LMC + MW')

py.xlabel(r'$H_0$ [km/s/Mpc]',fontsize=18)

#py.ylabel(r'$ \frac{\mathrm{Standard\, deviation}}{\mathrm{Mean\, value}}*100$')

py.tick_params(axis='both',which='major',labelsize=15)

minor_ticks = np.arange(65,81,1)

ax = fig.add_subplot(111)

ax.tick_params(axis='both',which='minor',labelsize=12)

ax.set_xticks(minor_ticks, minor = True)

py.gca().axes.get_yaxis().set_visible(False)

#py.legend(loc=0,numpoints=1,ncol=1)

py.tight_layout(pad=0, h_pad=0, w_pad=0)
py.savefig('./H0_values_three_anchors.pdf')

###### NGC4258 ###### 

fig = py.figure()

fig.subplots_adjust(bottom=0.1, left=0.06, top=0.85, right=0.97, hspace=0.2)

py.subplot(3,2,1)

print 'NGC :'
# Indirect determinations

py.errorbar(wmap09[0],2.25,yerr=None,xerr=wmap09[1],ecolor='#a710f9',marker=marker[0],mfc='#a710f9',fmt='',ls='None',label='WMAP 2009 (3.1%)')

print 'WMAP 2009 = ', wmap09[1]/wmap09[0]*100.
#py.errorbar(planck13[0],planck13[1]/planck13[0]*100.,yerr=None,xerr=planck13[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK13')

py.errorbar(planck15[0],1.5,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None')#,label='PLANCK 2015')
#py.annotate('1.4%',xy=(70,1.5),xycoords='data',size=12,color='black')
print 'Planck 2015 = ', planck15[1]/planck15[0]*100.

# Riess 2011

py.errorbar(ngc4258_R11[0],3.75,yerr=None,xerr=ngc4258_R11[1],ecolor='#f91fb4',marker=marker[0],mfc='#f91fb4',fmt='',ls='None')#,label='Riess et al. 2011')
py.annotate('4.0%',xy=(78.2,3.62),xycoords='data',size=10,color='#f91fb4')

print 'R11 = ', ngc4258_R11[1]/ngc4258_R11[0]*100.
# Efstathiou 2015 

py.errorbar(ngc4258_E15[0],4.5,yerr=None,xerr=ngc4258_E15[1],ecolor='dodgerblue',marker=marker[0],mfc='dodgerblue',fmt='',ls='None')#,label='Efstathiou 2015')
py.annotate('4.7%',xy=(74.2,4.37),xycoords='data',size=10,color='dodgerblue')

print 'E15 = ', ngc4258_E15[1]/ngc4258_E15[0]*100.
# Riess 2016

py.errorbar(ngc4258_R16[0],3.,yerr=None,xerr=ngc4258_R16[1],ecolor='#f91fb4',marker=marker[2],mfc='#f91fb4',fmt='',ls='None')#,label='Riess et al. 2016')
py.annotate('3.5%',xy=(75.,2.87),xycoords='data',size=10,color='#f91fb4')

print 'R16 = ', ngc4258_R16[1]/ngc4258_R16[0]*100.
# HP

py.errorbar(ngc4258_HP[0],5.25,yerr=None,xerr=ngc4258_HP[1],ecolor='#0c54c0',marker=marker[5],mfc='#0c54c0',fmt='',ls='None',label='HPs in Cepheids and SN Ia')
py.annotate('5.8%',xy=(76.5,5.12),xycoords='data',size=10,color='#0c54c0')

print 'HP = ', ngc4258_HP[1]/ngc4258_HP[0]*100.

py.xlim(66.,83.)

py.gca().axes.get_yaxis().set_visible(False)

py.ylim(1.,7.)

py.title('NGC 4258')

#py.xlabel(r'$H_0$ [km/s/Mpc]')

#py.ylabel(r'$ \frac{\mathrm{Standard\, deviation}}{\mathrm{Mean\, value}}*100$')

py.legend(prop={'size':8},numpoints=1,ncol=2,loc=0)

#py.savefig('./H0_values_ngc4258.pdf')

###### LMC #####

#fig = py.figure()

py.subplot(3,2,2)

# Indirect determinations
print 'LMC :'
py.errorbar(wmap09[0],2.25,yerr=None,xerr=wmap09[1],ecolor='#a710f9',marker=marker[0],mfc='#a710f9',fmt='',ls='None')#,label='WMAP 2009')

#py.errorbar(planck13[0],planck13[1]/planck13[0]*100.,yerr=None,xerr=planck13[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK13')

py.errorbar(planck15[0],1.5,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK 2015 (1.4%)')

# Riess 2011

py.errorbar(lmc_R11[0],4.5,yerr=None,xerr=lmc_R11[1],ecolor='#f91fb4',marker=marker[0],mfc='#f91fb4',fmt='',ls='None')#,label='Riess et al. 2011')
py.annotate('5.3%',xy=(75.2,4.37),xycoords='data',size=10,color='#f91fb4')
print 'R11 = ',lmc_R11[1]/lmc_R11[0]*100.

# Efstathiou 2015 

py.errorbar(lmc_E15[0],3.75,yerr=None,xerr=lmc_E15[1],ecolor='dodgerblue',marker=marker[0],mfc='dodgerblue',fmt='',ls='None')#,label='Efstathiou 2015')
py.annotate('4.2%',xy=(76.7,3.62),xycoords='data',size=10,color='dodgerblue')
print 'E15 =',lmc_E15[1]/lmc_E15[0]*100.

# Riess 2016

py.errorbar(lmc_R16[0],3.,yerr=None,xerr=lmc_R16[1],ecolor='#f91fb4',marker=marker[2],mfc='#f91fb4',fmt='',ls='None')#,label='Riess et al. 2016')
py.annotate('3.7%',xy=(75.,2.87),xycoords='data',size=10,color='#f91fb4')
print 'R16 =', lmc_R16[1]/lmc_R16[0]*100.
# HP

py.errorbar(lmc_HP[0],5.25,yerr=None,xerr=lmc_HP[1],ecolor='#0c54c0',marker=marker[5],mfc='#0c54c0',fmt='',ls='None')#,label='HPs in Cepheids, SN Ia, and distance modulus')
py.annotate('5.4%',xy=(78.5,5.12),xycoords='data',size=10,color='#0c54c0')
print 'HP =',lmc_HP[1]/lmc_HP[0]*100.

py.xlim(66.,83.)

py.gca().axes.get_yaxis().set_visible(False)

py.ylim(1.,7.)

py.title('LMC')

#py.xlabel(r'$H_0$ [km/s/Mpc]')

#py.ylabel(r'$ \frac{\mathrm{Standard\, deviation}}{\mathrm{Mean\, value}}*100$')

py.legend(numpoints=1,ncol=1,prop={'size':8},loc=2)
#py.legend(bbox_to_anchor=(-1.2,2.4,2.2,.102),loc=0,mode='expand',borderaxespad=0.,numpoints=1,ncol=1)

#py.savefig('./H0_values_lmc.pdf',bbox_inches='tight')

###### MW #####

#fig = py.figure()

py.subplot(3,2,3)

# Indirect determinations
print 'MW :'
py.errorbar(wmap09[0],2.25,yerr=None,xerr=wmap09[1],ecolor='#a710f9',marker=marker[0],mfc='#a710f9',fmt='',ls='None')#,label='WMAP 2009')

#py.errorbar(planck13[0],planck13[1]/planck13[0]*100.,yerr=None,xerr=planck13[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK13')

py.errorbar(planck15[0],1.5,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None')#,label='PLANCK 2015')

# Riess 2011

py.errorbar(mw_R11[0],3.75,yerr=None,xerr=mw_R11[1],ecolor='#f91fb4',marker=marker[0],mfc='#f91fb4',fmt='',ls='None',label='Riess et al. 2011')
py.annotate('3.4%',xy=(78.5,3.62),xycoords='data',size=10,color='#f91fb4')
print 'R11 = ',mw_R11[1]/mw_R11[0]*100.

# Efstathiou 2015 

py.errorbar(mw_E15[0],4.5,yerr=None,xerr=mw_E15[1],ecolor='dodgerblue',marker=marker[0],mfc='dodgerblue',fmt='',ls='None')#,label='Efstathiou 2015')
py.annotate('4.6%',xy=(79.,4.37),xycoords='data',size=10,color='dodgerblue')
print 'E15 =',mw_E15[1]/mw_E15[0]*100.

# Riess 2016

py.errorbar(mw_R16[0],3.,yerr=None,xerr=mw_R16[1],ecolor='#f91fb4',marker=marker[2],mfc='#f91fb4',fmt='',ls='None')#,label='Riess et al. 2016')
py.annotate('3.1%',xy=(79.,2.87),xycoords='data',size=10,color='#f91fb4')
print 'R16 =',mw_R16[1]/mw_R16[0]*100.
# HP

py.errorbar(mw_HP[0],5.25,yerr=None,xerr=mw_HP[1],ecolor='#0c54c0',marker=marker[5],mfc='#0c54c0',fmt='',ls='None')#,label='HPs in Cepheids, SN Ia, and distance modulus')
py.annotate('5.6%',xy=(72.,5.12),xycoords='data',size=10,color='#0c54c0')
print 'HP =',mw_HP[1]/mw_HP[0]*100.

py.xlim(66.,83.)

py.gca().axes.get_yaxis().set_visible(False)

py.ylim(1.,7.)

py.title('MW')

#py.xlabel(r'$H_0$ [km/s/Mpc]')

#py.ylabel(r'$ \frac{\mathrm{Standard\, deviation}}{\mathrm{Mean\, value}}*100$')

py.legend(loc=2,numpoints=1,ncol=1,prop={'size':8})

#py.savefig('./H0_values_mw.pdf')

# NGC 4258 + LMC

#fig = py.figure()

py.subplot(3,2,4)

# Indirect determinations
print 'NGC+LMC :'
py.errorbar(wmap09[0],2.25,yerr=None,xerr=wmap09[1],ecolor='#a710f9',marker=marker[0],mfc='#a710f9',fmt='',ls='None')#,label='WMAP 2009')

#py.errorbar(planck13[0],planck13[1]/planck13[0]*100.,yerr=None,xerr=planck13[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK13')

py.errorbar(planck15[0],1.5,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None')#,label='PLANCK 2015')

# Riess 2011

py.errorbar(ngc4258_lmc_R11[0],3.,yerr=None,xerr=ngc4258_lmc_R11[1],ecolor='#f91fb4',marker=marker[0],mfc='#f91fb4',fmt='',ls='None')#,label='Riess et al. 2011')
py.annotate('3.2%',xy=(75.,2.87),xycoords='data',size=10,color='#f91fb4')
print 'R11 =',ngc4258_lmc_R11[1]/ngc4258_lmc_R11[0]*100.

# Efstathiou 2015 

py.errorbar(ngc4258_lmc_E15[0],3.75,yerr=None,xerr=ngc4258_lmc_E15[1],ecolor='dodgerblue',marker=marker[0],mfc='dodgerblue',fmt='',ls='None',label='Efstathiou 2015')
py.annotate('3.6%',xy=(75.,3.62),xycoords='data',size=10,color='dodgerblue')
print 'E15 =',ngc4258_lmc_E15[1]/ngc4258_lmc_E15[0]*100.
# HPs

py.errorbar(ngc4258_lmc_HP[0],4.5,yerr=None,xerr=ngc4258_lmc_HP[1],ecolor='#0c54c0',marker=marker[0],mfc='#0c54c0',fmt='',ls='None')#,label='HPs in Cepheids, SN Ia, and distance modulus')
py.annotate('5.6%',xy=(76.,4.37),xycoords='data',size=10,color='#0c54c0')
print 'HP =',ngc4258_lmc_HP[1]/ngc4258_lmc_HP[0]*100.

py.xlim(66.,83.)

py.gca().axes.get_yaxis().set_visible(False)

py.ylim(1.,7.)

py.title('NGC 4258 + LMC')

#py.xlabel(r'$H_0$ [km/s/Mpc]')

#py.ylabel(r'$ \frac{\mathrm{Standard\, deviation}}{\mathrm{Mean\, value}}*100$')

py.legend(loc=2,numpoints=1,ncol=1,prop={'size':8})

#py.savefig('./H0_values_ngc4258_lmc.pdf')


# NGC 4258 + MW
print 'NGC + MW:'
#fig = py.figure()

py.subplot(3,2,5)

# Indirect determinations

py.errorbar(wmap09[0],2.25,yerr=None,xerr=wmap09[1],ecolor='#a710f9',marker=marker[0],mfc='#a710f9',fmt='',ls='None')#,label='WMAP 2009')

#py.errorbar(planck13[0],planck13[1]/planck13[0]*100.,yerr=None,xerr=planck13[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK13')

py.errorbar(planck15[0],1.5,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None')#,label='PLANCK 2015')

# Riess 2011

py.errorbar(ngc4258_mw_R11[0],3.75,yerr=None,xerr=ngc4258_mw_R11[1],ecolor='#f91fb4',marker=marker[0],mfc='#f91fb4',fmt='',ls='None')#,label='Riess et al. 2011')
py.annotate('3.1%',xy=(77.,3.62),xycoords='data',size=10,color='#f91fb4')
print 'R11 = ',ngc4258_mw_R11[1]/ngc4258_mw_R11[0]*100.

# Efstathiou 2015 

py.errorbar(ngc4258_mw_E15[0],4.5,yerr=None,xerr=ngc4258_mw_E15[1],ecolor='dodgerblue',marker=marker[0],mfc='dodgerblue',fmt='',ls='None')#,label='Efstathiou 2015')
py.annotate('3.9%',xy=(75.2,4.37),xycoords='data',size=10,color='dodgerblue')
print 'E15 = ',ngc4258_mw_E15[1]/ngc4258_mw_E15[0]*100.

# Riess 2016

py.errorbar(ngc4258_mw_R16[0],3.,yerr=None,xerr=ngc4258_mw_R16[1],ecolor='#f91fb4',marker=marker[2],mfc='#f91fb4',fmt='',ls='None',label='Riess et al. 2016')
py.annotate('2.6%',xy=(76.2,2.87),xycoords='data',size=10,color='#f91fb4')
print 'R16 = ',ngc4258_mw_R16[1]/ngc4258_mw_R16[0]*100.
# HP

py.errorbar(ngc4258_mw_HP[0],5.25,yerr=None,xerr=ngc4258_mw_HP[1],ecolor='#0c54c0',marker=marker[0],mfc='#0c54c0',fmt='',ls='None')#,label='HPs in Cepheids, SN Ia, and distance modulus')
py.annotate('5.2%',xy=(81.2,5.12),xycoords='data',size=10,color='#0c54c0')
print 'HP = ',ngc4258_mw_HP[1]/ngc4258_mw_HP[0]*100.

py.xlim(66.,83.)

py.gca().axes.get_yaxis().set_visible(False)

py.ylim(1.,7.)

py.title('NGC 4258 + MW')

py.xlabel(r'$H_0$ [km/s/Mpc]')

#py.ylabel(r'$ \frac{\mathrm{Standard\, deviation}}{\mathrm{Mean\, value}}*100$')

py.legend(loc=2,numpoints=1,ncol=1,prop={'size':8})

#py.savefig('./H0_values_ngc4258_mw.pdf')

# LMC + MW
print 'LMC + Mw:'
#fig = py.figure()

py.subplot(3,2,6)

# Indirect determinations

py.errorbar(wmap09[0],2.25,yerr=None,xerr=wmap09[1],ecolor='#a710f9',marker=marker[0],mfc='#a710f9',fmt='',ls='None')#,label='WMAP 2009')

#py.errorbar(planck13[0],planck13[1]/planck13[0]*100.,yerr=None,xerr=planck13[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK13')

py.errorbar(planck15[0],1.5,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None')#,label='PLANCK 2015')

# Riess 2011

py.errorbar(lmc_mw_R11[0],3.,yerr=None,xerr=lmc_mw_R11[1],ecolor='#f91fb4',marker=marker[0],mfc='#f91fb4',fmt='',ls='None')#,label='Riess et al. 2011')
py.annotate('3.2%',xy=(77.5,2.87),xycoords='data',size=10,color='#f91fb4')
print 'R11 = ',lmc_mw_R11[1]/lmc_mw_R11[0]*100.
# Efstathiou 2015 

py.errorbar(lmc_mw_E15[0],3.75,yerr=None,xerr=lmc_mw_E15[1],ecolor='dodgerblue',marker=marker[0],mfc='dodgerblue',fmt='',ls='None')#,label='Efstathiou 2015')
py.annotate('3.7%',xy=(77.,3.62),xycoords='data',size=10,color='dodgerblue')
print 'E15 = ', lmc_mw_E15[1]/lmc_mw_E15[0]*100.

# Riess 2016

#py.errorbar(ngc4258_lmc_[0],ngc4258_lmc_[1]/ngc4258_lmc_[0],yerr=None,xerr=ngc4258_lmc_[1],ecolor='red',marker=marker[4],mfc='red',fmt='',ls='None',label='NGC 4258 + LMC')

#py.errorbar(lmc_mw_[0],lmc_mw_[1]/lmc_mw_[0],yerr=None,xerr=lmc_mw_[1],ecolor='red',marker=marker[6],mfc='red',fmt='',ls='None',label='LMC + MW')

# HP

py.errorbar(lmc_mw_HP[0],4.5,yerr=None,xerr=lmc_mw_HP[1],ecolor='#0c54c0',marker=marker[0],mfc='#0c54c0',fmt='',ls='None',label='HPs in Cepheids, SN Ia, and distance modulus')
py.annotate('5.0%',xy=(80.2,4.37),xycoords='data',size=10,color='#0c54c0')
print 'HP = ',lmc_mw_HP[1]/lmc_mw_HP[0]*100.

py.xlim(66.,83.)

py.gca().axes.get_yaxis().set_visible(False)

py.ylim(1.,7.)

py.title('LMC + MW')

py.xlabel(r'$H_0$ [km/s/Mpc]')

#py.ylabel(r'$ \frac{\mathrm{Standard\, deviation}}{\mathrm{Mean\, value}}*100$')

py.legend(loc=2,numpoints=1,ncol=1,prop={'size':8})

#py.savefig('./H0_values_lmc_mw.pdf')
#py.legend(bbox_to_anchor=(-1.2,2.4,2.2,.102),loc=3,mode='expand',borderaxespad=0.,numpoints=1,ncol=1)

fig.tight_layout()

py.savefig('./H0_values_anchor_combination.pdf',bbox_inches='tight')

exit()

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



py.xlim(1.e0,2.e2)

py.xlabel('Period [days]')

py.ylabel('W residual [mag]')






py.savefig('../output/chains/effective_HP_cepheids.pdf')

exit()





