from getdist import loadMCSamples,plots

number_of_parameters = 14

samples = loadMCSamples('../output/chains/mcmc_final_output_HP',settings={'ignore_rows': 2 }) 

g = plots.getSinglePlotter()

g.settings.rcSizes(axes_fontsize = 2,lab_fontsize = 7)

g.triangle_plot(samples,filled=True)

g.export('../output/chains/triangle_figure_HP_R11_H.pdf')

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

stats.saveAsText('../output/chains/1Dstatistics_HP_R11_H.txt')

print 'Figure has been created '

exit()
