from getdist import loadMCSamples,plots,covmat
import numpy as np

number_of_parameters = 27 # 14 NGC4258 AS ANCHOR, 16 LMC AS ANCHOR, 15 MW AS ANCHOR, 16 NGC4258+LMC AS ANCHORS, 15 NGC4258+MW AS ANCHORS, 
# 16 LMC+MW AS ANCHORS, 16 NGC4258+LMC+MW AS ANCHORS, 19 NGC4258+LMC+MW AS ANCHORS (VARYING sigma_int), 27 NGC4258+LMC+MW AS ANCHORS AND SIGMA INT PER GALAXY, 26 NGC4258+LMC+MW AS ANCHORS AND SIGMA INT PER GALAXY WITHOUT METALLICITY

samples = loadMCSamples('../output/chains/previous_runs/fit_M1a/mcmc_final_output_HP',settings={'ignore_rows': 0.2 })

samples2 = loadMCSamples('../../../NGC4258-LMC-MW-R16/hypermcmc/output/chains/previous_runs/fit_44/mcmc_final_output_HP',settings={'ignore_rows': 0. }) 


g = plots.getSinglePlotter()

g.settings.rcSizes(axes_fontsize = 7,lab_fontsize = 12)

g.settings.x_label_rotation = -45.

if number_of_parameters == 26 :

    pass

else:

    g.triangle_plot([samples,samples2],['mu04258','muLMC','Mw','Zw','bw','H0'],filled_compare=True,legend_labels=['Cardona et al. 2016 (5.2%) [R11 data with 722 HP]','Cardona et al. 2016 (2.9%) [R16 data with 2298 HP]'],line_args=[{'ls':'dotted','color':'red'},{'lw':2,'color':'blue'}],contour_colors=['red','blue'],legend_loc='upper right',legend_fontsize=2)

    for ax in g.subplots[:,0]:
        ax.axvline(29.25,color='green',ls='--')
        ax.axvline(29.40,color='black',ls='--')
        ax.axvline(29.387,color='red',ls='--')

    
    for ax in g.subplots[1:,1]:
        ax.axvline(18.49,color='black',ls='--')

    for ax in g.subplots[5:,5]:
        ax.axvline(67.81,color='black',ls='--')
        ax.axvline(70.6,color='green',ls='--')
        ax.axvline(73.8,color='red',ls='--')
        ax.axvline(73.24,color='red',ls='dotted')


    g.export('triangle_plot_fit_29_44.pdf')

exit()
f = plots.getSinglePlotter()

f.settings.rcSizes(axes_fontsize = 5,lab_fontsize = 7)

f.settings.x_label_rotation = -45.

f.triangle_plot(samples,['mu01','mu02','mu03','mu04','mu05','mu06','mu07','mu08','Mw','bw','H0'],filled=True)

f.export('../output/chains/previous_runs/fit_M1a/triangle_figure_HP_R11_W_2.pdf')

p = samples.getParams()

samples.addDerived(p.mu01 - p.mu04258, name='mu01_mu04258', label='\mu_{0,1}-\mu_{0,4258}')

samples.addDerived(p.mu02 - p.mu04258, name='mu02_mu04258', label='\mu_{0,2}-\mu_{0,4258}')

samples.addDerived(p.mu03 - p.mu04258, name='mu03_mu04258', label='\mu_{0,3}-\mu_{0,4258}')

samples.addDerived(p.mu04 - p.mu04258, name='mu04_mu04258', label='\mu_{0,4}-\mu_{0,4258}')

samples.addDerived(p.mu05 - p.mu04258, name='mu05_mu04258', label='\mu_{0,5}-\mu_{0,4258}')

samples.addDerived(p.mu06 - p.mu04258, name='mu06_mu04258', label='\mu_{0,6}-\mu_{0,4258}')

samples.addDerived(p.mu07 - p.mu04258, name='mu07_mu04258', label='\mu_{0,7}-\mu_{0,4258}')

samples.addDerived(p.mu08 - p.mu04258, name='mu08_mu04258', label='\mu_{0,8}-\mu_{0,4258}')

samples.addDerived(p.mu04258 + 5.*np.log10(p.H0) - 25., name='mu04258_5av', label='m^0_{v,4258}+5a_v')

if number_of_parameters == 16 :

    pass

else:

    samples.addDerived(np.power(10,p.log10sigma_int_LMC),name='sigma_int_LMC',label='\sigma_{int}^{LMC}')

    samples.addDerived(np.power(10,p.log10sigma_int_MW),name='sigma_int_MW',label='\sigma_{int}^{MW}')

#if number_of_parameters == 27 :


#    samples.addDerived(np.power(10,p.log10sigma_int_1),name='sigma_int_1',label='\sigma_{int,1}')

#    samples.addDerived(np.power(10,p.log10sigma_int_2),name='sigma_int_2',label='\sigma_{int,2}')

#    samples.addDerived(np.power(10,p.log10sigma_int_3),name='sigma_int_3',label='\sigma_{int,3}')

#    samples.addDerived(np.power(10,p.log10sigma_int_4),name='sigma_int_4',label='\sigma_{int,4}')

#    samples.addDerived(np.power(10,p.log10sigma_int_5),name='sigma_int_5',label='\sigma_{int,5}')

#    samples.addDerived(np.power(10,p.log10sigma_int_6),name='sigma_int_6',label='\sigma_{int,6}')

#    samples.addDerived(np.power(10,p.log10sigma_int_7),name='sigma_int_7',label='\sigma_{int,7}')

#    samples.addDerived(np.power(10,p.log10sigma_int_8),name='sigma_int_8',label='\sigma_{int,8}')

#    samples.addDerived(np.power(10,p.log10sigma_int_9),name='sigma_int_9',label='\sigma_{int,9}')

#samples.addDerived(p.mu01 + 5.*np.log10(p.H0) - 25. - 5.*p.av, name='mv1', label='m_{v,1}')

#samples.addDerived(p.mu02 + 5.*np.log10(p.H0) - 25. - 5.*p.av, name='mv2', label='m_{v,2}')

#samples.addDerived(p.mu03 + 5.*np.log10(p.H0) - 25. - 5.*p.av, name='mv3', label='m_{v,3}')

#samples.addDerived(p.mu04 + 5.*np.log10(p.H0) - 25. - 5.*p.av, name='mv4', label='m_{v,4}')

#samples.addDerived(p.mu05 + 5.*np.log10(p.H0) - 25. - 5.*p.av, name='mv5', label='m_{v,5}')

#samples.addDerived(p.mu06 + 5.*np.log10(p.H0) - 25. - 5.*p.av, name='mv6', label='m_{v,6}')

#samples.addDerived(p.mu07 + 5.*np.log10(p.H0) - 25. - 5.*p.av, name='mv7', label='m_{v,7}')

#samples.addDerived(p.mu08 + 5.*np.log10(p.H0) - 25. - 5.*p.av, name='mv8', label='m_{v,8}')

#samples.addDerived(p.mu04258 + 5.*np.log10(p.H0) - 25. - 5.*p.av, name='mv4258', label='m_{v,4258}')

bestfit = samples.getLikeStats()

means = samples.setMeans()

filebestfit = open("../output/chains/previous_runs/fit_M1a/bestfit.txt",'w')

filemeans = open("../output/chains/previous_runs/fit_M1a/means.txt",'w')

for index in range(number_of_parameters) : 

    filebestfit.write(str(bestfit.names[index].bestfit_sample)+"\n")

    filemeans.write(str(means[index])+"\n")

filebestfit.close()

filemeans.close()

stats = samples.getMargeStats()

stats.saveAsText('../output/chains/previous_runs/fit_M1a/1Dstatistics_HP_R11_W.txt')

covariance_matrix = samples.getCov(nparam=number_of_parameters)

covariance_matrix_2 = covmat.CovMat(matrix=covariance_matrix)

covariance_matrix_2.saveToFile('../output/chains/previous_runs/fit_M1a/covariance_matrix.txt')

if number_of_parameters == 27 :

#    h = plots.getSubplotPlotter()

#    h.settings.rcSizes(axes_fontsize = 2, lab_fontsize = 7)

#    h.plots_1d(samples,['sigma_int_LMC','sigma_int_MW','sigma_int_1','sigma_int_2','sigma_int_3','sigma_int_4','sigma_int_5','sigma_int_6','sigma_int_7','sigma_int_8','sigma_int_9'],nx=3)

#    h.export('../output/chains/1D_sigma_int_plots.pdf')

    print 'ANALYZE SCRIPT ENDED'

else:

    print 'ANALYZE SCRIPT ENDED'

exit()
