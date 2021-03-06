import numpy as np

data = np.loadtxt('mcmc_output.txt')

#for index in range(len(data[:,4])):
#    data[index,4] = np.log10(data[index,4])

Cov = np.cov(data[:,2:],rowvar=0)

np.savetxt('covariance_matrix.txt',Cov,fmt='%.10E')

exit()
