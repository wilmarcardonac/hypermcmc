import numpy as np
import math

# Title: The Cepheid Period Luminosity relation in the Large Magellanic Cloud  
# Authors: Sebo K.M., Rawson D., Mould J., Madore B.F., Putman M.E., Graham J.A., 
# Freedman W.L., Gibson B.K., Germany L.M. 
# Table: LMC Cepheids  - Properties

# We load data from the paper above

table4 = np.dtype([('name',np.str_,20),('V',np.float32),('I',np.float32)])

n,V,I = np.loadtxt('table4.txt',unpack=True,usecols=[1,11,13],dtype=table4)

# Title: New Cepheid period-luminosity relations for the Large Magellanic Cloud: 92 near-infrared light curves
# Authors: S.E. Persson et al.
# Table: LMC Cepheid Infrared Magnitudes

# We load data from the last paper above

table3 = np.dtype([('name',np.str_,20),('Period',np.float32),('H',np.float32),('Sigma_m',np.float32)])

N,p,h,sm = np.loadtxt('table3_original.txt',unpack=True,usecols=[0,1,2,3],dtype=table3)

print math.log10(p[1])
# We create files to write data 

dataA = open('dataA.txt','w')

dataB = open('dataB.txt','w')

dataAB = open('dataAB.txt','w')

dataC = open('dataC.txt','w')

dataABC = open('dataABC.txt','w')

# An array with header

#data.write("#Name    Period    H    Sigma_m    V    I\n")

# counter = 0

for i in range(len(n)):
    for j in range(len(N)):

        if ( (N[j] == n[i]) and ( math.log10(p[j])<=1.0 )):
            dataA.write(N[j]+" "+str(p[j])+" "+str(h[j])+" "+str(sm[j])+" "+str(V[i])+" "+str(I[i])+"\n" )

        if ( (N[j] == n[i]) and ((math.log10(p[j])>1.0) and (math.log10(p[j])< 1.8))):
            dataB.write(N[j]+" "+str(p[j])+" "+str(h[j])+" "+str(sm[j])+" "+str(V[i])+" "+str(I[i])+"\n" )

        if ( (N[j] == n[i]) and ( math.log10(p[j])>= 1.8 ) ):
            dataC.write(N[j]+" "+str(p[j])+" "+str(h[j])+" "+str(sm[j])+" "+str(V[i])+" "+str(I[i])+"\n" )

        if ( (N[j] == n[i]) and ( math.log10(p[j])< 1.8 ) ):
            dataAB.write(N[j]+" "+str(p[j])+" "+str(h[j])+" "+str(sm[j])+" "+str(V[i])+" "+str(I[i])+"\n" )

        if ( N[j] == n[i] ):
            dataABC.write(N[j]+" "+str(p[j])+" "+str(h[j])+" "+str(sm[j])+" "+str(V[i])+" "+str(I[i])+"\n" )


dataA.close()
           
dataB.close()

dataAB.close()

dataC.close()

dataABC.close()
 
exit()
