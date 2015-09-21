from getdist import loadMCSamples,plots

number_model_parameters = 2

number_hyperparameters = 53

number_of_parameters = number_model_parameters + number_hyperparameters

samples = loadMCSamples('../output/mcmc_final_output_HP',settings={'ignore_rows':2}) 

g = plots.getSinglePlotter()

g.settings.rcSizes(axes_fontsize = 4,lab_fontsize = 7)

g.plot_2d(samples,'A','bw')

#g.triangle_plot(samples,filled=True)

g.export('2D_plot_HP_as_MCMC.pdf')

k = plots.getSubplotPlotter(width_inch=10)

k.settings.rcSizes(axes_fontsize = 7,lab_fontsize = 7)

k.plots_1d(samples,['alpha_01','alpha_02','alpha_03','alpha_04','alpha_05','alpha_06','alpha_07','alpha_08','alpha_09','alpha_10','alpha_11','alpha_12','alpha_13','alpha_14','alpha_15','alpha_16','alpha_17','alpha_18','alpha_19','alpha_20','alpha_21','alpha_22','alpha_23','alpha_24','alpha_25'],nx=5)

k.export('1D_plot_HP_as_MCMC_1-25.pdf')

k.plots_1d(samples,['alpha_26','alpha_27','alpha_28','alpha_29','alpha_30','alpha_31','alpha_32','alpha_33','alpha_34','alpha_35','alpha_36','alpha_37','alpha_38','alpha_39','alpha_40','alpha_41','alpha_42','alpha_43','alpha_44','alpha_45','alpha_46','alpha_47','alpha_48','alpha_49','alpha_50'],nx=5)

k.export('1D_plot_HP_as_MCMC_25-50.pdf')

k.plots_1d(samples,['alpha_51','alpha_52','alpha_53'],nx=2)

k.export('1D_plot_HP_as_MCMC_51-53.pdf')

bestfit = samples.getLikeStats()

means = samples.setMeans()

filebestfit = open("../output/bestfit.txt",'w')

filemeans = open("../output/means.txt",'w')

#filebestfit.write("-log(Like) = "+str(bestfit.logLike_sample)+"\n")

for index in range(number_of_parameters) : 

    filebestfit.write(str(bestfit.names[index].bestfit_sample)+"\n")

    filemeans.write(str(means[index])+"\n")

filebestfit.close()

filemeans.close()

stats = samples.getMargeStats()

stats.saveAsText('1Dstatistics_HP_as_MCMC.txt')

print 'Figure has been created '

exit()
