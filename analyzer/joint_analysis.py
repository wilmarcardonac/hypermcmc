from getdist import loadMCSamples,plots

number_of_parameters = 2

samples = loadMCSamples('../output/mcmc_final_output',settings={'ignore_rows': 2 }) 

samples_HP = loadMCSamples('../output/mcmc_final_output_HP',settings={'ignore_rows': 2 })

g = plots.getSubplotPlotter()

g.settings.rcSizes(axes_fontsize = 4,lab_fontsize = 7)

g.triangle_plot([samples,samples_HP],filled=True,legend_labels=[r'$\sum \chi_j^2$',r'$\sum N_j \ln (\chi_j^2)$'])

g.export('triangle_figure_joint.pdf')

print 'Figure has been created '

exit()
