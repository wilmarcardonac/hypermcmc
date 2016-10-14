from getdist import loadMCSamples,plots
import numpy as np

number_of_parameters = 4 # 2 IF LMC ALONE : 3; IF LMC ALONE AND INCLUDING METALLICITY : 4. ,  12, 14

samples = loadMCSamples('../output/chains/mcmc_final_output_HP',settings={'ignore_rows': 0.2 }) 

g = plots.getSinglePlotter()

g.settings.rcSizes(axes_fontsize = 4,lab_fontsize = 7)

g.triangle_plot(samples,filled=True)

g.export('../output/chains/triangle_figure_HP.pdf')

p = samples.getParams()

samples.addDerived(np.power(10,p.log10sigma_int),name='sigma_int',label='\sigma_{int}')

bestfit = samples.getLikeStats()

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

stats.saveAsText('../output/chains/1Dstatistics_HP.txt')

print 'Figure has been created '

exit()
