import numpy as np
import os
import matplotlib.pyplot as plt


prefix='ParRF10-25'

myindexes=range(1,46)+[76]+range(107,182)
nrec=len(myindexes)

for id1 in range(0,nrec):
   id=myindexes[id1]
   mytemplate='%s/%s-faultreceiver-00%03d*' %(prefix, prefix,id)
   print(mytemplate)
   f = os.popen("ls "+mytemplate)
   now = f.read()
   myfile=now.strip()
   print(myfile)
   test = np.loadtxt(myfile,  skiprows=8)
   if id1==0:
      ti=test[:,0]
      nti=np.shape(ti)[0]
      aSR = np.zeros((nrec,nti))
   aSR[id1,:]=test[:,1]
   print(np.max(aSR[id1,:]))
   
xi = np.zeros((nrec+1))
for id in range(0,nrec+1):
   xi[id]=-20000. + 333.*id



# Generate a regular grid to interpolate the data.
Xi, Ti = np.meshgrid(xi, ti)

print(xi)
print(Xi)
# Plot the results
plt.figure()

plt.pcolormesh(Xi,Ti,np.transpose(aSR))
print(aSR)
plt.colorbar()
plt.clim(0,2)
#plt.axis('equal')
plt.xlim(-20e3, 20e3)
plt.xlabel('along strike distance (m)')
plt.ylabel('time (s)')
plt.show()
