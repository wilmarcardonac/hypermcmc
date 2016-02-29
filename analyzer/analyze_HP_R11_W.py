from getdist import loadMCSamples,plots,covmat
import numpy as np

number_of_parameters = 26 # 14 NGC4258 AS ANCHOR, 16 LMC AS ANCHOR, 15 MW AS ANCHOR, 16 NGC4258+LMC AS ANCHORS, 15 NGC4258+MW AS ANCHORS, 
# 16 LMC+MW AS ANCHORS, 19 NGC4258+LMC+MW AS ANCHORS, 27 NGC4258+LMC+MW AS ANCHORS AND SIGMA INT PER GALAXY, 26 NGC4258+LMC+MW AS ANCHORS AND SIGMA INT PER GALAXY WITHOUT METALLICITY

samples = loadMCSamples('../output/chains/mcmc_final_output_HP',settings={'ignore_rows': 0.2 }) 

g = plots.getSinglePlotter()

g.settings.rcSizes(axes_fontsize = 2,lab_fontsize = 7)

g.triangle_plot(samples,['mu04258','muLMC','Mw','bw','H0','log10sigma_int_LMC','log10sigma_int_MW'],filled=True)

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

samples.addDerived(p.mu04258 + 5.*np.log10(p.H0) - 25., name='mu04258_5av', label='m^0_{v,4258}+5a_v')

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

filebestfit = open("../output/chains/bestfit.txt",'w')

filemeans = open("../output/chains/means.txt",'w')

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
