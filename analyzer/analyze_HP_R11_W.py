from getdist import loadMCSamples,plots,covmat

number_of_parameters = 16 # 14 NGC4258 AS ANCHOR, 16 LMC AS ANCHOR, 15 MW AS ANCHOR, 16 NGC4258+LMC AS ANCHORS, 15 NGC4258+MW AS ANCHORS, 
# 16 LMC+MW AS ANCHORS, 16 NGC4258+LMC+MW AS ANCHORS

samples = loadMCSamples('../output/chains/mcmc_final_output_HP',settings={'ignore_rows': 0.2 }) 

g = plots.getSinglePlotter()

g.settings.rcSizes(axes_fontsize = 2,lab_fontsize = 7)

g.triangle_plot(samples,filled=True)

g.export('../output/chains/triangle_figure_HP_R11_W.pdf')

p = samples.getParams()

samples.addDerived(p.mu01 - p.mu04258, name='mu01_mu04258', label='\mu_{0,1}-\mu_{0,4258}')

samples.addDerived(p.mu02 - p.mu04258, name='mu02_mu04258', label='\mu_{0,2}-\mu_{0,4258}')

samples.addDerived(p.mu03 - p.mu04258, name='mu03_mu04258', label='\mu_{0,3}-\mu_{0,4258}')

samples.addDerived(p.mu04 - p.mu04258, name='mu04_mu04258', label='\mu_{0,4}-\mu_{0,4258}')

samples.addDerived(p.mu05 - p.mu04258, name='mu05_mu04258', label='\mu_{0,5}-\mu_{0,4258}')

samples.addDerived(p.mu06 - p.mu04258, name='mu06_mu04258', label='\mu_{0,6}-\mu_{0,4258}')

samples.addDerived(p.mu07 - p.mu04258, name='mu07_mu04258', label='\mu_{0,7}-\mu_{0,4258}')

samples.addDerived(p.mu08 - p.mu04258, name='mu08_mu04258', label='\mu_{0,8}-\mu_{0,4258}')

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

stats.saveAsText('../output/chains/1Dstatistics_HP_R11_W.txt')

covariance_matrix = samples.getCov(nparam=number_of_parameters)

covariance_matrix_2 = covmat.CovMat(matrix=covariance_matrix)

covariance_matrix_2.saveToFile('../output/chains/covariance_matrix.txt')

print 'ANALYZE SCRIPT ENDED'

exit()
