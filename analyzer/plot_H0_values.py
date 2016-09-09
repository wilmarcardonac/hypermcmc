import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

#host = ['LMC']

# Indirect determinations 

wmap09  = [70.0,2.2]

#planck13 = [67.9,1.5]

planck15 = [67.81,0.92]

# Riess 2011

ngc4258_R11 = [74.8,3.0]

lmc_R11 = [71.3,3.8]

mw_R11 = [75.7,2.6]

ngc4258_lmc_R11 = [72.3,2.3]

ngc4258_mw_R11 = [74.5,2.3]

lmc_mw_R11 = [74.4,2.4]

ngc4258_lmc_mw_R11 = [73.8,2.2]

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

ngc4258_lmc_mw_m31_R16 = [73.46,1.71]

# HP

#ngc4258_HP = [72.1,4.2]

#lmc_HP = [74.4,4.0]

#mw_HP = [78.1,4.4]

#ngc4258_lmc_HP = [71.1,4.0]

#ngc4258_mw_HP = [76.9,4.0]

#lmc_mw_HP = [76.1,3.8]

ngc4258_lmc_mw_HP_29 = [75.0,3.9]

ngc4258_lmc_mw_m31_HP_43 = [73.88,2.16]

#ngc4258_lmc_mw_HP_M1ag = [73.2,2.5]

ngc4258_lmc_mw_m31_HP_56 = [73.54,1.56]

ngc4258_lmc_mw_HP_33 = [72.4,2.2]

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

#a = py.errorbar(wmap09[0],wmap09[1]/wmap09[0]*100.,yerr=None,xerr=wmap09[1],ecolor='blue',marker=marker[0],mfc='blue',fmt='',ls='None',label='WMAP 2009',elinewidth=3,ms=10)

print 'WMAP ', wmap09[1]/wmap09[0]*100.

#a = py.errorbar(wmap09[0],1,yerr=None,xerr=wmap09[1],ecolor='blue',marker=marker[0],mfc='blue',fmt='',ls='None',label='WMAP 2009',elinewidth=3,ms=10)

#l1 = py.legend(bbox_to_anchor=(0,.27,1,.1),loc=3,numpoints=1,handles=[a])

#py.errorbar(planck13[0],planck13[1]/planck13[0]*100.,yerr=None,xerr=planck13[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK13')

#b = py.errorbar(planck15[0],planck15[1]/planck15[0]*100.,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK 2015',elinewidth=3,ms=10)

print 'Planck ', planck15[1]/planck15[0]*100.

#b = py.errorbar(planck15[0],2,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK 2015',elinewidth=3,ms=10)

#l2 = py.legend(bbox_to_anchor=(0.25,.0,1,.1),loc=3,numpoints=1,handles=[b])

#py.gca().add_artist(l1)

# Riess 2011

#c = py.errorbar(ngc4258_lmc_mw_R11[0],ngc4258_lmc_mw_R11[1]/ngc4258_lmc_mw_R11[0]*100.,yerr=None,xerr=ngc4258_lmc_mw_R11[1],ecolor='red',marker=marker[0],mfc='red',fmt='',ls='None',label='Riess et al. 2011',elinewidth=3,ms=10)

print 'R11 ', ngc4258_lmc_mw_R11[1]/ngc4258_lmc_mw_R11[0]*100.

#c = py.errorbar(ngc4258_lmc_mw_R11[0],3,yerr=None,xerr=ngc4258_lmc_mw_R11[1],ecolor='red',marker=marker[0],mfc='red',fmt='',ls='None',label='Riess et al. 2011',elinewidth=3,ms=10)

#l3 = py.legend(bbox_to_anchor=(0.57,.38,1,.1),loc=3,numpoints=1,handles=[c])

#py.gca().add_artist(l2)

# Efstathiou 2015 

#d = py.errorbar(ngc4258_E15[0],ngc4258_E15[1]/ngc4258_E15[0]*100.,yerr=None,xerr=ngc4258_E15[1],ecolor='magenta',marker=marker[1],mfc='magenta',fmt='',ls='None',label='Efstathiou 2015',elinewidth=3,ms=10)

print 'E15 ', ngc4258_E15[1]/ngc4258_E15[0]*100. 

#d = py.errorbar(ngc4258_E15[0],4,yerr=None,xerr=ngc4258_E15[1],ecolor='magenta',marker=marker[1],mfc='magenta',fmt='',ls='None',label='Efstathiou 2015',elinewidth=3,ms=10)

#l4 = py.legend(bbox_to_anchor=(0.55,.7,1,.1),loc=3,numpoints=1,handles=[d])

#py.gca().add_artist(l3)

# Riess 2016

#e = py.errorbar(ngc4258_lmc_mw_R16[0],ngc4258_lmc_mw_R16[1]/ngc4258_lmc_mw_R16[0]*100.,yerr=None,xerr=ngc4258_lmc_mw_R16[1],ecolor='green',marker=marker[2],mfc='green',fmt='',ls='None',label='Riess et al. 2016',elinewidth=3,ms=10)

print 'R16 ', ngc4258_lmc_mw_R16[1]/ngc4258_lmc_mw_R16[0]*100.

#e = py.errorbar(ngc4258_lmc_mw_R16[0],5,yerr=None,xerr=ngc4258_lmc_mw_R16[1],ecolor='green',marker=marker[2],mfc='green',fmt='',ls='None',label='Riess et al. 2016',elinewidth=3,ms=10)

#l5 = py.legend(bbox_to_anchor=(0.,.32,1,.1),loc=3,numpoints=1,handles=[e])

#py.gca().add_artist(l2)

# HP

#f = py.errorbar(ngc4258_lmc_mw_HP_29[0],ngc4258_lmc_mw_HP_29[1]/ngc4258_lmc_mw_HP_29[0]*100.,yerr=None,xerr=ngc4258_lmc_mw_HP_29[1],ecolor='dodgerblue',marker=marker[3],mfc='dodgerblue',fmt='',ls='None',label=r'Cardona et al. 2016 (R11 data with HP)',elinewidth=3,ms=10)

print 'R11 data with HP ', ngc4258_lmc_mw_HP_29[1]/ngc4258_lmc_mw_HP_29[0]*100.

#f = py.errorbar(ngc4258_lmc_mw_HP_29[0],6,yerr=None,xerr=ngc4258_lmc_mw_HP_29[1],ecolor='dodgerblue',marker=marker[3],mfc='dodgerblue',fmt='',ls='None',label=r'Cardona et al. 2016 (R11 data with HP)',elinewidth=3,ms=10)

#l6 = py.legend(bbox_to_anchor=(0,.8,1,.1),loc=3,numpoints=1,handles=[f])

#py.gca().add_artist(l5)

#g = py.errorbar(ngc4258_lmc_mw_HP_33[0],ngc4258_lmc_mw_HP_33[1]/ngc4258_lmc_mw_HP_33[0]*100.,yerr=None,xerr=ngc4258_lmc_mw_HP_33[1],ecolor='dodgerblue',marker=marker[6],mfc='dodgerblue',fmt='',ls='None',label=r'Cardona et al. 2016 (R11 data with HP in Cepheids)',elinewidth=3,ms=10)

print 'R11 data with HP in Cepheids', ngc4258_lmc_mw_HP_33[1]/ngc4258_lmc_mw_HP_33[0]*100.

#g = py.errorbar(ngc4258_lmc_mw_HP_33[0],7,yerr=None,xerr=ngc4258_lmc_mw_HP_33[1],ecolor='dodgerblue',marker=marker[6],mfc='dodgerblue',fmt='',ls='None',label=r'Cardona et al. 2016 (R11 data with HP in Cepheids)',elinewidth=3,ms=10)

#l7 = py.legend(bbox_to_anchor=(0,.7,1,.1),loc=3,numpoints=1,handles=[g])

#py.gca().add_artist(l6)

#h = py.errorbar(ngc4258_lmc_mw_m31_HP_43[0],ngc4258_lmc_mw_m31_HP_43[1]/ngc4258_lmc_mw_m31_HP_43[0]*100.,yerr=None,xerr=ngc4258_lmc_mw_m31_HP_43[1],ecolor='dodgerblue',marker=marker[4],mfc='dodgerblue',fmt='',ls='None',label=r'Cardona et al. 2016 (R16 data with HP)',elinewidth=3,ms=10)

print 'R16 data with HP', ngc4258_lmc_mw_m31_HP_43[1]/ngc4258_lmc_mw_m31_HP_43[0]*100.

#h = py.errorbar(ngc4258_lmc_mw_m31_HP_43[0],8,yerr=None,xerr=ngc4258_lmc_mw_m31_HP_43[1],ecolor='dodgerblue',marker=marker[4],mfc='dodgerblue',fmt='',ls='None',label=r'Cardona et al. 2016 (R16 data with HP)',elinewidth=3,ms=10)

#l8 = py.legend(bbox_to_anchor=(0,.86,1,.1),loc=3,numpoints=1,handles=[f])

#py.gca().add_artist(l5)

#l9 = py.legend(bbox_to_anchor=(0,.5,1,.1),loc=3,numpoints=1,handles=[d,a,g,c,h])

#py.gca().add_artist(l8)

#j = py.errorbar(ngc4258_lmc_mw_m31_HP_56[0],ngc4258_lmc_mw_m31_HP_56[1]/ngc4258_lmc_mw_m31_HP_56[0]*100.,yerr=None,xerr=ngc4258_lmc_mw_m31_HP_56[1],ecolor='dodgerblue',marker=marker[5],mfc='dodgerblue',fmt='',ls='None',label=r'Cardona et al. 2016 (R16 data with HP in Cepheids)',elinewidth=3,ms=10)

print 'R16 data with HP in Cepheids ', ngc4258_lmc_mw_m31_HP_56[1]/ngc4258_lmc_mw_m31_HP_56[0]*100.

#j = py.errorbar(ngc4258_lmc_mw_m31_HP_56[0],9,yerr=None,xerr=ngc4258_lmc_mw_m31_HP_56[1],ecolor='dodgerblue',marker=marker[5],mfc='dodgerblue',fmt='',ls='None',label=r'Cardona et al. 2016 (R16 data with HP in Cepheids)',elinewidth=3,ms=10)

#l10 = py.legend(bbox_to_anchor=(0,.13,1,.1),loc=3,numpoints=1,handles=[j])

#py.gca().add_artist(l9)

#mpl.rcParams['legend.handlelength'] = 0

b = py.errorbar(planck15[0],1,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='Planck 2015',elinewidth=3,ms=10)

#j = py.errorbar(ngc4258_lmc_mw_m31_HP_56[0],10,yerr=None,xerr=ngc4258_lmc_mw_m31_HP_56[1],ecolor='blue',marker=marker[5],mfc='blue',fmt='',ls='None',label=r'Cardona et al. 2016 (R16 data with HP in Cepheids)',elinewidth=3,ms=10)

e = py.errorbar(ngc4258_lmc_mw_R16[0],19,yerr=None,xerr=ngc4258_lmc_mw_R16[1],ecolor='green',marker=marker[2],mfc='green',fmt='',ls='None',label='Riess et al. 2016',elinewidth=3,ms=10)

h = py.errorbar(ngc4258_lmc_mw_m31_HP_43[0],28,yerr=None,xerr=ngc4258_lmc_mw_m31_HP_43[1],ecolor='blue',marker=marker[4],mfc='blue',fmt='',ls='None',label=r'Cardona et al. 2016 (R16 data with HP)',elinewidth=3,ms=10)

c = py.errorbar(ngc4258_lmc_mw_R11[0],37,yerr=None,xerr=ngc4258_lmc_mw_R11[1],ecolor='red',marker=marker[0],mfc='red',fmt='',ls='None',label='Riess et al. 2011',elinewidth=3,ms=10)

#g = py.errorbar(ngc4258_lmc_mw_HP_33[0],46,yerr=None,xerr=ngc4258_lmc_mw_HP_33[1],ecolor='blue',marker=marker[6],mfc='blue',fmt='',ls='None',label=r'Cardona et al. 2016 (R11 data with HP in Cepheids)',elinewidth=3,ms=10)

a = py.errorbar(wmap09[0],55,yerr=None,xerr=wmap09[1],ecolor='dodgerblue',marker=marker[0],mfc='dodgerblue',fmt='',ls='None',label='WMAP 2009',elinewidth=3,ms=10)

d = py.errorbar(ngc4258_E15[0],64,yerr=None,xerr=ngc4258_E15[1],ecolor='magenta',marker=marker[1],mfc='magenta',fmt='',ls='None',label='Efstathiou 2015',elinewidth=3,ms=10)

f = py.errorbar(ngc4258_lmc_mw_HP_29[0],73,yerr=None,xerr=ngc4258_lmc_mw_HP_29[1],ecolor='blue',marker=marker[3],mfc='blue',fmt='',ls='None',label=r'Cardona et al. 2016 (R11 data with HP)',elinewidth=3,ms=10)

py.annotate('Planck 2015 (1.4%)', xy=(66,3), xycoords='data', size=14,color='black')

#py.annotate('Cardona et al. 2016 (2.1%) [R16 data with HP in Cepheids]', xy=(66,12), xycoords='data', size=14,color='blue')

py.annotate('Riess et al. 2016 (2.4%)', xy=(66,21), xycoords='data', size=14,color='green')

py.annotate('Cardona et al. 2016 (2.9%) [R16 data with HP]', xy=(66,30), xycoords='data', size=14,color='blue')

py.annotate('Riess et al. 2011 (3.0%)', xy=(66,39), xycoords='data', size=14,color='red')

#py.annotate('Cardona et al. 2016 (3.0%) [R11 data with HP in Cepheids]', xy=(66,48), xycoords='data', size=14,color='blue')

py.annotate('WMAP 2009 (3.1%)', xy=(66,57), xycoords='data', size=14,color='dodgerblue')

py.annotate('Efstathiou 2015 (4.7%)', xy=(66,66), xycoords='data', size=14,color='magenta')

py.annotate('Cardona et al. 2016 (5.2%) [R11 data with HP]', xy=(66,75), xycoords='data', size=14,color='blue')

#l9 = py.legend(bbox_to_anchor=(0,.4,1,.1),loc=3,numpoints=1,handles=[f,d,a,g,c,h,e,j,b],frameon=False,markerscale=0,handlelength=0,handletextpad=0,fancybox=True)

#for item in l9.legendHandles:
#    item.set_visible(False)

py.xlim(65.,83.)

py.ylim(-2.,81.)

#py.suptitle(r'Best estimates of $H_0$',fontsize=15)

#py.title('ordered by uncertainty (but equidistant)')

py.xlabel(r'$H_0$ [km/s/Mpc]')

py.gca().axes.get_yaxis().set_visible(False)

#py.ylabel(r'$ \frac{\mathrm{Standard\, deviation}}{\mathrm{Mean\, value}}*100$')

#py.legend(loc=0,numpoints=1,ncol=1)

py.savefig('./H0_best_estimates.eps')

exit()
# NGC4258

fig = py.figure()

fig.subplots_adjust(bottom=0.1, left=0.06, top=0.85, right=0.97, hspace=0.2)

py.subplot(3,2,1)

# Indirect determinations

py.errorbar(wmap09[0],wmap09[1]/wmap09[0]*100.,yerr=None,xerr=wmap09[1],ecolor='blue',marker=marker[0],mfc='blue',fmt='',ls='None',label='WMAP 2009')

#py.errorbar(planck13[0],planck13[1]/planck13[0]*100.,yerr=None,xerr=planck13[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK13')

py.errorbar(planck15[0],planck15[1]/planck15[0]*100.,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None')#,label='PLANCK 2015')

# Riess 2011

py.errorbar(ngc4258_R11[0],ngc4258_R11[1]/ngc4258_R11[0]*100.,yerr=None,xerr=ngc4258_R11[1],ecolor='red',marker=marker[0],mfc='red',fmt='',ls='None')#,label='Riess et al. 2011')

# Efstathiou 2015 

py.errorbar(ngc4258_E15[0],ngc4258_E15[1]/ngc4258_E15[0]*100.,yerr=None,xerr=ngc4258_E15[1],ecolor='magenta',marker=marker[1],mfc='magenta',fmt='',ls='None')#,label='Efstathiou 2015')

# Riess 2016

py.errorbar(ngc4258_R16[0],ngc4258_R16[1]/ngc4258_R16[0]*100.,yerr=None,xerr=ngc4258_R16[1],ecolor='green',marker=marker[2],mfc='green',fmt='',ls='None')#,label='Riess et al. 2016')

# HP

py.errorbar(ngc4258_HP[0],ngc4258_HP[1]/ngc4258_HP[0]*100.,yerr=None,xerr=ngc4258_HP[1],ecolor='dodgerblue',marker=marker[5],mfc='dodgerblue',fmt='',ls='None',label='HPs in Cepheids and SN Ia')

py.xlim(66.,79.)

py.ylim(1.,8.)

py.title('NGC 4258')

#py.xlabel(r'$H_0$ [km/s/Mpc]')

py.ylabel(r'$ \frac{\mathrm{Standard\, deviation}}{\mathrm{Mean\, value}}*100$')

py.legend(prop={'size':8},numpoints=1,ncol=2,loc=0)

#py.savefig('./H0_values_ngc4258.pdf')

# LMC

#fig = py.figure()

py.subplot(3,2,2)

# Indirect determinations

py.errorbar(wmap09[0],wmap09[1]/wmap09[0]*100.,yerr=None,xerr=wmap09[1],ecolor='blue',marker=marker[0],mfc='blue',fmt='',ls='None')#,label='WMAP 2009')

#py.errorbar(planck13[0],planck13[1]/planck13[0]*100.,yerr=None,xerr=planck13[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK13')

py.errorbar(planck15[0],planck15[1]/planck15[0]*100.,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK 2015')

# Riess 2011

py.errorbar(lmc_R11[0],lmc_R11[1]/lmc_R11[0]*100.,yerr=None,xerr=lmc_R11[1],ecolor='red',marker=marker[0],mfc='red',fmt='',ls='None')#,label='Riess et al. 2011')

# Efstathiou 2015 

py.errorbar(lmc_E15[0],lmc_E15[1]/lmc_E15[0]*100.,yerr=None,xerr=lmc_E15[1],ecolor='magenta',marker=marker[1],mfc='magenta',fmt='',ls='None')#,label='Efstathiou 2015')

# Riess 2016

py.errorbar(lmc_R16[0],lmc_R16[1]/lmc_R16[0]*100.,yerr=None,xerr=lmc_R16[1],ecolor='green',marker=marker[2],mfc='green',fmt='',ls='None')#,label='Riess et al. 2016')

# HP

py.errorbar(lmc_HP[0],lmc_HP[1]/lmc_HP[0]*100.,yerr=None,xerr=lmc_HP[1],ecolor='dodgerblue',marker=marker[5],mfc='dodgerblue',fmt='',ls='None')#,label='HPs in Cepheids, SN Ia, and distance modulus')

py.xlim(66.,79.)

py.ylim(1.,7.)

py.title('LMC')

#py.xlabel(r'$H_0$ [km/s/Mpc]')

#py.ylabel(r'$ \frac{\mathrm{Standard\, deviation}}{\mathrm{Mean\, value}}*100$')

py.legend(numpoints=1,ncol=1,prop={'size':8},loc=4)
#py.legend(bbox_to_anchor=(-1.2,2.4,2.2,.102),loc=0,mode='expand',borderaxespad=0.,numpoints=1,ncol=1)

#py.savefig('./H0_values_lmc.pdf',bbox_inches='tight')

# MW

#fig = py.figure()

py.subplot(3,2,3)

# Indirect determinations

py.errorbar(wmap09[0],wmap09[1]/wmap09[0]*100.,yerr=None,xerr=wmap09[1],ecolor='blue',marker=marker[0],mfc='blue',fmt='',ls='None')#,label='WMAP 2009')

#py.errorbar(planck13[0],planck13[1]/planck13[0]*100.,yerr=None,xerr=planck13[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK13')

py.errorbar(planck15[0],planck15[1]/planck15[0]*100.,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None')#,label='PLANCK 2015')

# Riess 2011

py.errorbar(mw_R11[0],mw_R11[1]/mw_R11[0]*100.,yerr=None,xerr=mw_R11[1],ecolor='red',marker=marker[0],mfc='red',fmt='',ls='None',label='Riess et al. 2011')

# Efstathiou 2015 

py.errorbar(mw_E15[0],mw_E15[1]/mw_E15[0]*100.,yerr=None,xerr=mw_E15[1],ecolor='magenta',marker=marker[1],mfc='magenta',fmt='',ls='None')#,label='Efstathiou 2015')

# Riess 2016

py.errorbar(mw_R16[0],mw_R16[1]/mw_R16[0]*100.,yerr=None,xerr=mw_R16[1],ecolor='green',marker=marker[2],mfc='green',fmt='',ls='None')#,label='Riess et al. 2016')

# HP

py.errorbar(mw_HP[0],mw_HP[1]/mw_HP[0]*100.,yerr=None,xerr=mw_HP[1],ecolor='dodgerblue',marker=marker[5],mfc='dodgerblue',fmt='',ls='None')#,label='HPs in Cepheids, SN Ia, and distance modulus')

py.xlim(66.,83.)

py.ylim(1.,7.)

py.title('MW')

#py.xlabel(r'$H_0$ [km/s/Mpc]')

py.ylabel(r'$ \frac{\mathrm{Standard\, deviation}}{\mathrm{Mean\, value}}*100$')

py.legend(loc=2,numpoints=1,ncol=1,prop={'size':8})

#py.savefig('./H0_values_mw.pdf')

# NGC 4258 + LMC

#fig = py.figure()

py.subplot(3,2,4)

# Indirect determinations

py.errorbar(wmap09[0],wmap09[1]/wmap09[0]*100.,yerr=None,xerr=wmap09[1],ecolor='blue',marker=marker[0],mfc='blue',fmt='',ls='None')#,label='WMAP 2009')

#py.errorbar(planck13[0],planck13[1]/planck13[0]*100.,yerr=None,xerr=planck13[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK13')

py.errorbar(planck15[0],planck15[1]/planck15[0]*100.,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None')#,label='PLANCK 2015')

# Riess 2011

py.errorbar(ngc4258_lmc_R11[0],ngc4258_lmc_R11[1]/ngc4258_lmc_R11[0]*100.,yerr=None,xerr=ngc4258_lmc_R11[1],ecolor='red',marker=marker[0],mfc='red',fmt='',ls='None')#,label='Riess et al. 2011')

# Efstathiou 2015 

py.errorbar(ngc4258_lmc_E15[0],ngc4258_lmc_E15[1]/ngc4258_lmc_E15[0]*100.,yerr=None,xerr=ngc4258_lmc_E15[1],ecolor='magenta',marker=marker[1],mfc='magenta',fmt='',ls='None',label='Efstathiou 2015')

# HPs

py.errorbar(ngc4258_lmc_HP[0],ngc4258_lmc_HP[1]/ngc4258_lmc_HP[0]*100.,yerr=None,xerr=ngc4258_lmc_HP[1],ecolor='dodgerblue',marker=marker[3],mfc='dodgerblue',fmt='',ls='None')#,label='HPs in Cepheids, SN Ia, and distance modulus')

py.xlim(66.,79.)

py.ylim(1.,7.)

py.title('NGC 4258 + LMC')

#py.xlabel(r'$H_0$ [km/s/Mpc]')

#py.ylabel(r'$ \frac{\mathrm{Standard\, deviation}}{\mathrm{Mean\, value}}*100$')

py.legend(loc=4,numpoints=1,ncol=1,prop={'size':8})

#py.savefig('./H0_values_ngc4258_lmc.pdf')


# NGC 4258 + MW

#fig = py.figure()

py.subplot(3,2,5)

# Indirect determinations

py.errorbar(wmap09[0],wmap09[1]/wmap09[0]*100.,yerr=None,xerr=wmap09[1],ecolor='blue',marker=marker[0],mfc='blue',fmt='',ls='None')#,label='WMAP 2009')

#py.errorbar(planck13[0],planck13[1]/planck13[0]*100.,yerr=None,xerr=planck13[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK13')

py.errorbar(planck15[0],planck15[1]/planck15[0]*100.,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None')#,label='PLANCK 2015')

# Riess 2011

py.errorbar(ngc4258_mw_R11[0],ngc4258_mw_R11[1]/ngc4258_mw_R11[0]*100.,yerr=None,xerr=ngc4258_mw_R11[1],ecolor='red',marker=marker[0],mfc='red',fmt='',ls='None')#,label='Riess et al. 2011')

# Efstathiou 2015 

py.errorbar(ngc4258_mw_E15[0],ngc4258_mw_E15[1]/ngc4258_mw_E15[0]*100.,yerr=None,xerr=ngc4258_mw_E15[1],ecolor='magenta',marker=marker[1],mfc='magenta',fmt='',ls='None')#,label='Efstathiou 2015')

# Riess 2016

py.errorbar(ngc4258_mw_R16[0],ngc4258_mw_R16[1]/ngc4258_mw_R16[0]*100.,yerr=None,xerr=ngc4258_mw_R16[1],ecolor='green',marker=marker[2],mfc='green',fmt='',ls='None',label='Riess et al. 2016')

# HP

py.errorbar(ngc4258_mw_HP[0],ngc4258_mw_HP[1]/ngc4258_mw_HP[0]*100.,yerr=None,xerr=ngc4258_mw_HP[1],ecolor='dodgerblue',marker=marker[3],mfc='dodgerblue',fmt='',ls='None')#,label='HPs in Cepheids, SN Ia, and distance modulus')

py.xlim(66.,82.)

py.ylim(1.,7.)

py.title('NGC 4258 + MW')

py.xlabel(r'$H_0$ [km/s/Mpc]')

py.ylabel(r'$ \frac{\mathrm{Standard\, deviation}}{\mathrm{Mean\, value}}*100$')

py.legend(loc=2,numpoints=1,ncol=1,prop={'size':8})

#py.savefig('./H0_values_ngc4258_mw.pdf')

# LMC + MW

#fig = py.figure()

py.subplot(3,2,6)

# Indirect determinations

py.errorbar(wmap09[0],wmap09[1]/wmap09[0]*100.,yerr=None,xerr=wmap09[1],ecolor='blue',marker=marker[0],mfc='blue',fmt='',ls='None')#,label='WMAP 2009')

#py.errorbar(planck13[0],planck13[1]/planck13[0]*100.,yerr=None,xerr=planck13[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None',label='PLANCK13')

py.errorbar(planck15[0],planck15[1]/planck15[0]*100.,yerr=None,xerr=planck15[1],ecolor='black',marker=marker[0],mfc='black',fmt='',ls='None')#,label='PLANCK 2015')

# Riess 2011

py.errorbar(lmc_mw_R11[0],lmc_mw_R11[1]/lmc_mw_R11[0]*100.,yerr=None,xerr=lmc_mw_R11[1],ecolor='red',marker=marker[0],mfc='red',fmt='',ls='None')#,label='Riess et al. 2011')

# Efstathiou 2015 

py.errorbar(lmc_mw_E15[0],lmc_mw_E15[1]/lmc_mw_E15[0]*100.,yerr=None,xerr=lmc_mw_E15[1],ecolor='magenta',marker=marker[1],mfc='magenta',fmt='',ls='None')#,label='Efstathiou 2015')

# Riess 2016

#py.errorbar(ngc4258_lmc_[0],ngc4258_lmc_[1]/ngc4258_lmc_[0],yerr=None,xerr=ngc4258_lmc_[1],ecolor='red',marker=marker[4],mfc='red',fmt='',ls='None',label='NGC 4258 + LMC')

#py.errorbar(lmc_mw_[0],lmc_mw_[1]/lmc_mw_[0],yerr=None,xerr=lmc_mw_[1],ecolor='red',marker=marker[6],mfc='red',fmt='',ls='None',label='LMC + MW')

# HP

py.errorbar(lmc_mw_HP[0],lmc_mw_HP[1]/lmc_mw_HP[0]*100.,yerr=None,xerr=lmc_mw_HP[1],ecolor='dodgerblue',marker=marker[3],mfc='dodgerblue',fmt='',ls='None',label='HPs in Cepheids, SN Ia, and distance modulus')

py.xlim(66.,81.)

py.ylim(1.,7.)

py.title('LMC + MW')

py.xlabel(r'$H_0$ [km/s/Mpc]')

#py.ylabel(r'$ \frac{\mathrm{Standard\, deviation}}{\mathrm{Mean\, value}}*100$')

py.legend(loc=0,numpoints=1,ncol=1,prop={'size':8})

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





