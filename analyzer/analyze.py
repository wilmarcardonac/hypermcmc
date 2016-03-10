from getdist import loadMCSamples,plots
import numpy as np

number_of_parameters = 3

samples = loadMCSamples('../output/chains/mcmc_final_output',settings={'ignore_rows': .2 }) 

g = plots.getSinglePlotter()

g.settings.rcSizes(axes_fontsize = 4,lab_fontsize = 7)

g.triangle_plot(samples,filled=True)

g.export('../output/chains/triangle_figure.pdf')

p = samples.getParams()

samples.addDerived(np.power(10,p.log10sigma_int),name='sigma_int',label='\sigma_{int}')

bestfit = samples.getLikeStats()

print 'BEST FIT Log like', bestfit.logLike_sample

means = samples.setMeans()

filebestfit = open("../output/chains/bestfit.txt",'w')

filemeans = open("../output/chains/means.txt",'w')

#filebestfit.write("-log(Like) = "+str(bestfit.logLike_sample)+"\n")

for index in range(number_of_parameters) : 

    filebestfit.write(str(bestfit.names[index].bestfit_sample)+"\n")

    filemeans.write(str(means[index])+"\n")

filebestfit.close()

filemeans.close()

stats = samples.getMargeStats()

stats.saveAsText('../output/chains/1Dstatistics.txt')

f = plots.getSinglePlotter()

f.settings.rcSizes(axes_fontsize = 4,lab_fontsize = 7)

f.triangle_plot(samples,['A','bw','sigma_int'],filled=True)

f.export('../output/chains/triangle_figure_sigma_int.pdf')

print 'Figure has been created '

exit()
