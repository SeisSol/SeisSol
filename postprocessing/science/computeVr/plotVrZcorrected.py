import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata
import scipy.ndimage as ndimage


#Compute a rupture Velocity map from the contatenated RF file
#the fault has to be parametrizable in function of x and z
#in fact t and y are mapped in function of the structured grid xi, zi 
#in order to apply the gradient function of numpy


#####Parameters###########

cplotfile="ParRF100-600-RF-concat.dat"
# Size of regular grid (in m)
dx, dz = 40., 40.
#GDmethod = 'cubic'
GDmethod = 'linear'
#use Gaussian Filter?
useGaussianFilter = True
#Show a few percentiles of the Vr distribution for helping
#setting p1 and p2
ShowTailVrDistribution = False
#min max percentiles of V for plotting (eliminate artefacts)
p1, p2=1, 99
#plot Vr=f(x) and Vr=f(y)
plotHist = False
#boundaries values for Vr plot
xm,xp = -20e3, 20e3
zm,zp = -20e3, 0
##########################

#Read Cplot file
xyt = np.loadtxt(cplotfile,  skiprows=1)
print('done reading %s' %cplotfile)
x = xyt[:,0]
y = xyt[:,1]
z = xyt[:,2]
t = xyt[:,3]


# Generate a regular grid to interpolate the data.
Xi = np.arange(min(x), max(x), dx)
Zi = np.arange(min(z), max(z), dz)
xi, zi = np.meshgrid(Xi, Zi)

# Interpolate using delaunay triangularization 
ti = griddata((x,z), t, (xi, zi), method=GDmethod, fill_value = 1e20)
yi = griddata((x,z), y, (xi, zi), method=GDmethod, fill_value = 0.)


if useGaussianFilter:
   # Increase the value of sigma to increase the amount of blurring.
   # order=0 means gaussian kernel
   ti = ndimage.gaussian_filter(ti, sigma=2.0, order=0)
   yi = ndimage.gaussian_filter(yi, sigma=2.0, order=0)

grad = np.gradient(ti)
gradY = np.gradient(yi)

dy1 = gradY[0]
dy2 = gradY[1]

dtdx=grad[0]/np.sqrt(pow(dx,2)+pow(dy1,2))
dtdz=grad[1]/np.sqrt(pow(dx,2)+pow(dy2,2))

slowness = np.sqrt(np.square(dtdx) + np.square(dtdz))
V=1./ slowness

#process data where NaN
#where_are_NaNs = np.isnan(V)
#V[where_are_NaNs] = 0
#where_are_null = np.where(V<1.)
#V[where_are_null] = np.nan

#Show a few percentiles for helping setting up p1 and p2
if ShowTailVrDistribution:
   for i in range(1,20):
     V1=np.percentile(V, i)
     print("percentile %d: %f" %(i,V1))

   print(" ")
   for i in range(99,70,-1):
     V1=np.percentile(V, i)
     print("percentile %d: %f" %(i,V1))

V1=np.percentile(V, p1)
V50=np.percentile(V, 50)
V2=np.percentile(V, p2)

print("percentiles %d, 50, %d: %d %d %d" %(p1,p2,V1,V50,V2))

#plt.hist(V, bins=[1000*i for i in range(0,10)])
#plt.show()

# Plot the results
plt.figure()

masked_array=np.ma.masked_where(V>1e10, V)
cmap = cm.jet
cmap.set_bad('w',1.)
plt.pcolormesh(xi,zi,V)
#Eliminate Vr artefacts (at rupture rim)
plt.pcolormesh(xi,zi,masked_array,cmap=cmap)
#plt.clim(V1,V2)
plt.clim(0,5400)
plt.colorbar()

CS = plt.contour(Xi, Zi, ti,range(1,21),colors='k')
plt.clabel(CS, fontsize=9, inline=1, fmt='%d')
plt.xlim(xm,xp)
plt.ylim(zm-2e3,zp+2e3)
plt.axis('equal')

plt.show()

if plotHist:
	n=np.shape(V)[0]
	lV1, lV50, lV2 = np.zeros(n), np.zeros(n), np.zeros(n)
	for i in range(n):
	    subV=V[i,:]
	    lV1[i]=np.percentile(subV, 33)
	    lV50[i]=np.percentile(subV, 50)
	    lV2[i]=np.percentile(subV, 67)
	plt.plot(lV1, Zi, label = '33%')
	plt.plot(lV50, Zi, label = '50%')
	plt.plot(lV2, Zi, label = '67%')
	plt.legend()
	plt.title("Vr=f(y)")
	plt.show()

	n=np.shape(V)[1]
	lV1, lV50, lV2 = np.zeros(n), np.zeros(n), np.zeros(n)
	for i in range(n):
	    subV=V[:,i]
	    lV1[i]=np.percentile(subV, 33)
	    lV50[i]=np.percentile(subV, 50)
	    lV2[i]=np.percentile(subV, 67)
	plt.plot(Xi, lV1, label = '33%')
	plt.plot(Xi, lV50, label = '50%')
	plt.plot(Xi, lV2, label = '67%')
	plt.legend()
	plt.title("Vr=f(x)")
	plt.show()

