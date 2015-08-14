import numpy as np

Data = np.loadtxt('mcmc_output.txt')

Cov = np.cov(Data[:,2:],rowvar=0)

np.savetxt('covariance_matrix.txt',Cov,fmt='%.10E')

exit()
