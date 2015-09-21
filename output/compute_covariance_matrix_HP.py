import numpy as np

data = np.loadtxt('mcmc_output.txt')

#for index1 in range(4,len(data[10,:])):
#    for index in range(len(data[:,index1])):
#        data[index,index1] = np.log10(data[index,index1])

Cov = np.cov(data[:,2:],rowvar=0)

np.savetxt('covariance_matrix.txt',Cov,fmt='%.10E')

exit()
