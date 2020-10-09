#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is written for SeisSol, to view the fault output of dynamic rupture models.
Currently the parameter_viewlist={'ASl' 'Vr' 'stressdrop' 'magnitude' 'sliprate' 'PSR'};
To run this, you need to change the model directory and the xdmfFilename. 
You may also need to adjust the view_azimuth to have a better view
This will plot the whole faults, including the ones not break
If you only want to plot the fautls broke, e.g. stress_drop[N_belowthreshold]=0, just replace the "0" with float("nan")
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

dir='/hppfs/scratch/0A/ru94tiz3/CHEESE/HFF_simple_closed/West_M7.301_R055_fd06_Us06_Ud01_Dc0206_SH145_Nr15_plastomu_co10/output_o4/'
xdmfFilename=(dir + 'HFFtest-fault.xdmf')
view_azimuth=-80
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

len_x=np.max(x)-np.min(x)
len_y=np.max(y)-np.min(y)
len_z=np.max(z)-np.min(z)

nElements = ReadNElements(xdmfFilename)
ndt=ReadNdt(xdmfFilename)

##default figsize is (8,6)
#This would work if the view elemention is not far away from 0
fig_x=len_x*np.sin(view_azimuth/180*np.pi)
fig_y=len_y*np.cos(view_azimuth/180*np.pi)
fig_length=np.sqrt(fig_x ** 2 +fig_y ** 2)/np.min([len_x,len_y,len_z])*4
fig_width=len_z/np.min([len_x,len_y,len_z])*4
fontsize_title=1.4*np.sqrt(fig_length ** 2 + fig_width ** 2)
fontsize_colorbar=1*np.sqrt(fig_length ** 2 + fig_width ** 2)

print ('start plotting the slip')
ASl = LoadData(xdmfFilename,'ASl',nElements, idt=ndt-1, oneDtMem=True)
ASl0=ASl[0]
N_abovethreshold=np.argwhere(ASl0>=0.1).ravel()
N_belowthreshold=np.argwhere(ASl0<0.1).ravel()
ASl_mean=statistics.mean(ASl0[N_abovethreshold])
fig = plt.figure(figsize=(fig_length,fig_width))
ax = fig.gca(projection='3d')
collec = ax.plot_trisurf(x, y, z, triangles=connect, cmap=plt.cm.jet, linewidth=0.2)
collec.set_array(ASl0)
collec.set_clim(0,np.percentile(ASl0[N_abovethreshold],95))
plt.axis('off')
plt.title('Total slip (m) mean = '+ str(round(ASl_mean*100)/100) + ' m',fontsize= fontsize_title)
ax.view_init(10, view_azimuth)
cbar = fig.colorbar(collec,orientation="horizontal")
cbar.ax.tick_params(labelsize=fontsize_colorbar)
plt.savefig('ASl.png')

###plot the Vr
print('start plotting the rupture velocity')
Vr = LoadData(xdmfFilename,'Vr',nElements, idt=ndt-1, oneDtMem=True)
Vr0=Vr[0]
Vr0[N_belowthreshold]=0 ##replace "0" with float("nan") if you only want to view the fault ruptured
Vr_mean=statistics.mean(Vr0[N_abovethreshold])

fig = plt.figure(figsize=(fig_length,fig_width))
ax = fig.gca(projection='3d')
collec = ax.plot_trisurf(x, y, z, triangles=connect, cmap=plt.cm.jet, linewidth=0.2)
collec.set_array(Vr0)
collec.set_clim(0,np.percentile(Vr0[N_abovethreshold],95))
plt.axis('off')
plt.title('Rupture Velocity (m/s) mean = '+ str(round(Vr_mean)) + ' m/s',fontsize= fontsize_title)
ax.view_init(10, view_azimuth)
cbar = fig.colorbar(collec,orientation="horizontal")
cbar.ax.tick_params(labelsize=fontsize_colorbar)
plt.savefig('Vr_faultview.png')

fig = plt.figure()
plt.hist(Vr0[N_abovethreshold], bins = np.linspace(0, np.percentile(Vr0[N_abovethreshold],99), endpoint=True, num=20))
plt.xlabel('Velocity (m/s)',fontsize=15)
plt.ylabel('Element counts',fontsize=15)
plt.title('Rupture Velocity Distribution Histogram \n '+'Vr mean = '+ str(round(Vr_mean)) + ' m/s',fontsize= 15)
plt.savefig('Vr_histogram.png')

print ('start plotting the stress drop')
T_s = LoadData(xdmfFilename,'T_s',nElements, idt=ndt-1, oneDtMem=True)
T_s0=T_s[0]
T_d = LoadData(xdmfFilename,'T_d',nElements, idt=ndt-1, oneDtMem=True)
T_d0=T_d[0]
stress_drop=np.sqrt(T_s0.ravel() ** 2 + T_d0.ravel() ** 2) / (1000000.0)
stress_drop[N_belowthreshold]=0 ##replace "0" with float("nan") if you only want to view the fault ruptured
stress_drop_mean=statistics.mean(stress_drop[N_abovethreshold])
fig = plt.figure(figsize=(fig_length,fig_width))
ax = fig.gca(projection='3d')
collec = ax.plot_trisurf(x, y, z, triangles=connect, cmap=plt.cm.jet, linewidth=0.2)
collec.set_array(stress_drop)
collec.set_clim(0,np.percentile(stress_drop[N_abovethreshold],95))
plt.axis('off')
plt.title('stressdrop (Mpa) mean= '+ str(round(stress_drop_mean*100)/100) + ' Mpa',fontsize= fontsize_title)
ax.view_init(10, view_azimuth)
cbar = fig.colorbar(collec,orientation="horizontal")
cbar.ax.tick_params(labelsize=fontsize_colorbar)
plt.savefig('stress_drop.png')

print ('start plotting the peak slip rate')
PSR = LoadData(xdmfFilename,'PSR',nElements, idt=ndt-1, oneDtMem=True)
PSR0=PSR[0]
PSR0[N_belowthreshold]=0 ##replace "0" with float("nan") if you only want to view the fault ruptured
PSR_mean=statistics.mean(PSR0[N_abovethreshold])
fig = plt.figure(figsize=(fig_length,fig_width))
ax = fig.gca(projection='3d')
collec = ax.plot_trisurf(x, y, z, triangles=connect, cmap=plt.cm.jet, linewidth=0.2)
collec.set_array(PSR0)
collec.set_clim(0,np.percentile(PSR0[N_abovethreshold],95))
plt.axis('off')
plt.title('Peak Slip Rate (m/s) mean = '+ str(round(PSR_mean*100)/100) + ' m/s',fontsize= fontsize_title)
ax.view_init(10, view_azimuth)
cbar = fig.colorbar(collec,orientation="horizontal")
cbar.ax.tick_params(labelsize=fontsize_colorbar)
plt.savefig('PSR.png')

print ('start plotting the Slip Rate')
step_time=0.25  ##time step for fault output
SR_window=np.round(np.linspace(1,ndt,endpoint=True,num=10))
SR_window=SR_window.astype(int)
##use the 3 window to help set the color bar limit, and use it for all frames
SRd = LoadData(xdmfFilename,'SRd',nElements, idt=SR_window[3]-1, oneDtMem=True)
SRd0=SRd[0]
SRs = LoadData(xdmfFilename,'SRs',nElements, idt=SR_window[3]-1, oneDtMem=True)
SRs0=SRs[0]
SR_plot=np.sqrt(SRs0.ravel() ** 2 + SRd0.ravel() ** 2)
SR_colorlimit=np.percentile(SR_plot[N_abovethreshold],98);

###select the middle 2-7 time step to show slip rate, change to your own case
subplot_row=3
subplot_col=2
fig = plt.figure(figsize=(fig_length*subplot_col,fig_width*subplot_row))
for i in range(len(SR_window)-4):  
        SRd = LoadData(xdmfFilename,'SRd',nElements, idt=SR_window[i+1]-1, oneDtMem=True)
        SRd0=SRd[0]
        SRs = LoadData(xdmfFilename,'SRs',nElements, idt=SR_window[i+1]-1, oneDtMem=True)
        SRs0=SRs[0]
        slip_rate=np.sqrt(SRs0.ravel() ** 2 + SRd0.ravel() ** 2)
        ax = fig.add_subplot(subplot_row,subplot_col,i+1,projection='3d')
        collec = ax.plot_trisurf(x, y, z, triangles=connect, cmap=plt.cm.jet, linewidth=0.2)
        collec.set_array(slip_rate)
        collec.set_clim(0,SR_colorlimit)
        plt.title('SR at '+ str(round(step_time*SR_window[i+1]*100)/100) + 's',fontsize=fontsize_title)
        ax.view_init(10, view_azimuth)
        plt.axis('off')
        cbar = fig.colorbar(collec,orientation="horizontal")
        cbar.ax.tick_params(labelsize=fontsize_colorbar)
plt.savefig('SR.png')
