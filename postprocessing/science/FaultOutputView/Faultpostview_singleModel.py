#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is written for SeisSol, to view the fault output for dynamic rupture models using SeisSol.
Currently the parameter_viewlist={'ASl' 'Vr' 'stressdrop' 'magnitude' 'sliprate' 'PSR'};
To run this, you need to change the model directory and the xdmfFilename. 
You may also need to adjust the view_azimuth to have a better view
@author: boli
"""
import sys
sys.modules[__name__].__dict__.clear()
from submodules.pythonXdmfReader.pythonXdmfReader import *
import glob
import numpy as np
import statistics
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

dir='/import/freenas-m-05-seissol/bo/CHEESE/PSHA/HFF/simple_model/h5test/smoothed_close/West_M7.301_R055_fd06_Us06_Ud01_Dc0206_SH145_Nr15_plastomu_co10/output_o4/'
xdmfFilename=(dir + 'HFFtest-fault.xdmf')
view_azimuth = -80
###moment rate and magnitue
print ('start plotting the moment rate & calculating the magnitude')
filelist=glob.glob(dir+'*-EnF_t*')
En1=np.loadtxt(filelist[1],skiprows=1)
dt=np.mean(np.diff(En1[:,0]))
Mo_rate=np.zeros([En1.shape[0],1])

for i,fname in enumerate(filelist):
    En=np.loadtxt(filelist[i],skiprows=1)
    Mo_rate[:,0]=Mo_rate[:,0]+En[:,1]
Mw=2/3*np.log10(np.sum(Mo_rate*dt))-6.07
print(Mw)

fig = plt.figure()
plt.plot(En1[:,0],Mo_rate[:,0])
plt.xlabel('Time (s)',fontsize=15)
plt.ylabel('Moment Rate',fontsize=15)
plt.title('Mw = '+ str(round(Mw*1000)/1000),fontsize=15)
plt.savefig('Moment_rate.png')

connect = ReadConnect(xdmfFilename)
geo = ReadGeometry(xdmfFilename)
x=geo[:,0]
y=geo[:,1]
z=geo[:,2]

ele_x=x[np.ravel(connect)]  ### ele_x=x[connect.ravel()]
ele_y=y[np.ravel(connect)]
ele_z=z[np.ravel(connect)]

nElements = ReadNElements(xdmfFilename)
ndt=ReadNdt(xdmfFilename)  ##Here we are plotting the last time step

print ('start plotting the slip')
ASl = LoadData(xdmfFilename,'ASl',nElements, idt=ndt-1, oneDtMem=True)
ASl0=ASl[0]
ASl0[ASl0<0.1]=float("nan")
ASl1=ASl0[~np.isnan(ASl0)]
ASl_mean=statistics.mean(ASl1)
fig = plt.figure()
ax = fig.gca(projection='3d')
collec = ax.plot_trisurf(x, y, z, triangles=connect, cmap=plt.cm.jet, linewidth=0.2)
collec.set_array(ASl0)
collec.set_clim(0,np.percentile(ASl1,95))
plt.title('Total slip (m) \n '+'ASl mean = '+ str(round(ASl_mean*100)/100) + ' m',fontsize=15)
ax.view_init(10, view_azimuth)
cbar = fig.colorbar(collec)
fig
plt.savefig('ASl.png')

###plot the Vr
print('start plotting the rupture velocity')
Vr = LoadData(xdmfFilename,'Vr',nElements, idt=ndt-1, oneDtMem=True)
Vr0=Vr[0]
Vr0[ASl0<0.1]=float("nan")
Vr0[Vr0<100]=float("nan")
Vr0[Vr0>10000]=float("nan") ##exclude some abnormal high values
Vr1=Vr0[~np.isnan(Vr0)]
Vr_mean=statistics.mean(Vr1)

fig = plt.figure()
ax = fig.gca(projection='3d')
collec = ax.plot_trisurf(x, y, z, triangles=connect, cmap=plt.cm.jet, linewidth=0.2)
collec.set_array(Vr0)
collec.set_clim(0,np.percentile(Vr1,95))
plt.title('Rupture Velocity (m/s) \n '+'Vr mean = '+ str(round(Vr_mean)) + ' m/s',fontsize=15)
ax.view_init(10, view_azimuth)
cbar = fig.colorbar(collec)
fig
plt.savefig('Vr_faultview.png')

fig = plt.figure()
plt.hist(Vr1, bins = np.linspace(0, np.percentile(Vr1,99), endpoint=True, num=20))
plt.xlabel('Velocity (m/s)',fontsize=15)
plt.ylabel('Element counts',fontsize=15)
plt.title('Rupture Velocity Distribution Histogram \n '+'Vr mean = '+ str(round(Vr_mean)) + ' m/s',fontsize=15)
fig
plt.savefig('Vr_histogram.png')

print ('start plotting the stress drop')
T_s = LoadData(xdmfFilename,'T_s',nElements, idt=ndt-1, oneDtMem=True)
T_s0=T_s[0]
T_d = LoadData(xdmfFilename,'T_d',nElements, idt=ndt-1, oneDtMem=True)
T_d0=T_d[0]
stress_drop=np.sqrt(T_s0.ravel() ** 2 + T_d0.ravel() ** 2) / (1000000.0)
stress_drop[ASl0<0.1] = float("nan")
stress_drop1=stress_drop[~np.isnan(stress_drop)]
stress_drop_mean=statistics.mean(stress_drop1)
fig = plt.figure()
ax = fig.gca(projection='3d')
collec = ax.plot_trisurf(x, y, z, triangles=connect, cmap=plt.cm.jet, linewidth=0.2)
collec.set_array(stress_drop)
collec.set_clim(0,np.percentile(stress_drop1,95))
plt.title('mean stressdrop (Mpa) = '+ str(round(stress_drop_mean*100)/100) + ' Mpa',fontsize=15)
ax.view_init(10, view_azimuth)
cbar = fig.colorbar(collec)
fig
plt.savefig('stress_drop.png')

print ('start plotting the peak slip rate')
PSR = LoadData(xdmfFilename,'PSR',nElements, idt=ndt-1, oneDtMem=True)
PSR0=PSR[0]
PSR0[ASl0<0.1] = float("nan")
PSR1=PSR0[~np.isnan(PSR0)]
PSR_mean=statistics.mean(PSR1)
fig = plt.figure()
ax = fig.gca(projection='3d')
collec = ax.plot_trisurf(x, y, z, triangles=connect, cmap=plt.cm.jet, linewidth=0.2)
collec.set_array(PSR0)
collec.set_clim(0,np.percentile(PSR1,95))
plt.title('Peak Slip Rate (m/s) \n '+'PSR mean = '+ str(round(PSR_mean*100)/100) + ' m/s',fontsize=15)
ax.view_init(10, view_azimuth)
cbar = fig.colorbar(collec)
fig
plt.savefig('PSR.png')

print ('start plotting the Slip Rate')
step_time=0.25  ##time step for fault output
SR_window=np.round(np.linspace(1,ndt,endpoint=True,num=10))
SR_window=SR_window.astype(int)
##use the 3 window to help set the color bar limit
SRd = LoadData(xdmfFilename,'SRd',nElements, idt=SR_window[3]-1, oneDtMem=True)
SRd0=SRd[0]
SRs = LoadData(xdmfFilename,'SRs',nElements, idt=SR_window[3]-1, oneDtMem=True)
SRs0=SRs[0]
SR_plot=np.sqrt(SRs0.ravel() ** 2 + SRd0.ravel() ** 2)
SR_plot[SR_plot<0.0001]=float("nan")
SR_plot1=SR_plot[~np.isnan(SR_plot)]
SR_colorlimit=np.percentile(SR_plot1,98);

###select the middle 2-7 time step to show slip rate, change to your own case
fig = plt.figure()
for i in range(len(SR_window)-4):  
        SRd = LoadData(xdmfFilename,'SRd',nElements, idt=SR_window[i+1]-1, oneDtMem=True)
        SRd0=SRd[0]
        SRs = LoadData(xdmfFilename,'SRs',nElements, idt=SR_window[i+1]-1, oneDtMem=True)
        SRs0=SRs[0]
        slip_rate=np.sqrt(SRs0.ravel() ** 2 + SRd0.ravel() ** 2)
        ax = fig.add_subplot(3,2,i+1,projection='3d')
        collec = ax.plot_trisurf(x, y, z, triangles=connect, cmap=plt.cm.jet, linewidth=0.2)
        collec.set_array(slip_rate)
        collec.set_clim(0,SR_colorlimit)
        plt.title('SR at '+ str(round(step_time*SR_window[i+1]*100)/100) + 's',fontsize=10)
        ax.view_init(10, view_azimuth)
        plt.axis('off')
        cbar = fig.colorbar(collec)
fig
plt.savefig('SR.png')

