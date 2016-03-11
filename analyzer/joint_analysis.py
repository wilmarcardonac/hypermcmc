from getdist import loadMCSamples,plots
import numpy as np

number_of_parameters = 3

samples = loadMCSamples('../output/chains/previous_runs/fit_d/mcmc_final_output',settings={'ignore_rows': .2 }) 

samples_HP = loadMCSamples('../output/chains/previous_runs/fit_c/mcmc_final_output_HP',settings={'ignore_rows': .2 })

g = plots.getSubplotPlotter()

g.settings.rcSizes(axes_fontsize = 4,lab_fontsize = 7)

g.triangle_plot(samples,filled=True,legend_labels=[r'$\sum \chi_j^2$',r'$\sum N_j \ln (\chi_j^2)$'])

g.triangle_plot([samples,samples_HP],filled=True,legend_labels=[r'$\sum \chi_j^2$',r'$\sum N_j \ln (\chi_j^2)$'])

g.export('triangle_figure_joint.pdf')

p = samples.getParams()

samples.addDerived(np.power(10,p.log10sigma_int),name='sigma_int',label='\sigma_{int}')

bestfit = samples.getLikeStats()

stats = samples.getMargeStats()

q = samples_HP.getParams()

samples_HP.addDerived(np.power(10,q.log10sigma_int),name='sigma_int',label='\sigma_{int}')

bestfit_HP = samples_HP.getParams()

stats_HP = samples_HP.getMargeStats()

f = plots.getSubplotPlotter()

f.settings.rcSizes(axes_fontsize = 4,lab_fontsize = 7)

f.triangle_plot([samples,samples_HP],['A','bw','sigma_int'],filled=True,legend_labels=[r'$\sum \chi_j^2$',r'$\sum N_j \ln (\chi_j^2)$'])

f.export('triangle_figure_joint_sigma_int.pdf')

print 'Figure has been created '

exit()
