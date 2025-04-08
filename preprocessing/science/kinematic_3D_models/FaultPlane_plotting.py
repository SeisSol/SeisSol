import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import math
import os

cmap = cm.get_cmap('viridis')

def roundup(x):
    return int(np.sign(x)*math.ceil(np.abs(x) / 10.0)) * 10

def plot_slip(fault1, fault2, savefig=None):
    A1, A2 = fault1.slip1, fault1.slip2
    AX, AY = fault1.x, fault1.y
    B1, B2 = fault2.slip1, fault2.slip2
    BX, BY = fault2.x, fault2.y
    
    p1max = roundup(max(A1.max(), B1.max()))
    p1min = roundup(min(A1.min(), B1.min()))
    p2max = roundup(max(A2.max(), B2.max()))
    p2min = roundup(min(A2.min(), B2.min()))

    fig, ax = plt.subplots(nrows=2,ncols=2,figsize=(9,9))
    ax[0,0].scatter(AX, AY, s=0.3 ,c='grey',label='Upsampled grid')
    ax[0,0].legend(loc='lower right')
    ax[0,0].scatter(AX, AY, c=A1,cmap=cmap,vmin=p1min, vmax=p1max)
    ax[0,1].contourf(AX, AY, A1, cmap=cmap, vmin=p1min, vmax=p1max,levels=20)

    a = ax[1,0].scatter(BX, BY, s=1.5, c=B1,cmap=cmap,vmin=p1min, vmax=p1max)
    ax[1,1].contourf(BX, BY, B1,cmap=cmap, vmin=p1min, vmax=p1max, levels=20)
    fig.suptitle('Along strike slip')
    fig.tight_layout()
    fig.subplots_adjust(right=0.93)
    cbar_ax = fig.add_axes([0.95, 0.2, 0.02, 0.6])
    fig.colorbar(a, cax=cbar_ax,label='slip (m)')
    
    if savefig is not None:
        fig.savefig(savefig+'_alongstrike.png')

    fig, ax = plt.subplots(nrows=2,ncols=2,figsize=(9,9))
    ax[0,0].scatter(BX, BY, s=0.2 ,c='grey')
    ax[0,0].scatter(AX, AY, c=A2, cmap=cmap, vmin=0, vmax=p2max)
    ax[0,1].contourf(AX, AY, A2, cmap=cmap, vmin=0, vmax=p2max,levels=20)

    a = ax[1,0].scatter(BX, BY, s=1.5, c=B2, cmap=cmap,vmin=0, vmax=p2max)
    ax[1,1].contourf(BX, BY, B2, cmap=cmap,vmin=0, vmax=p2max, levels=20)
    fig.suptitle('Along dip slip')
    fig.tight_layout()
    fig.subplots_adjust(right=0.93)
    cbar_ax = fig.add_axes([0.95, 0.2, 0.02, 0.6])
    fig.colorbar(a, cax=cbar_ax,label='slip (m)')
    
    if savefig is not None:
        fig.savefig(savefig+'_alongdip.png')

def plot_slip2(fault1, A1, A2, fault2, savefig=None):
    """
    Compare two slip
    Input:
        fault1 - fault object
        A1, A2 - slip1, slip2,
        fault2 - fault object
    """
    AX, AY = fault1.x, fault1.y
    B1, B2 = fault2.slip1, fault2.slip2
    BX, BY = fault2.x, fault2.y
    
    p1max = roundup(max(A1.max(), B1.max()))
    p1min = roundup(min(A1.min(), B1.min()))
    p2max = roundup(max(A2.max(), B2.max()))
    p2min = roundup(min(A2.min(), B2.min()))

    fig, ax = plt.subplots(nrows=2,ncols=2,figsize=(9,9))
    ax[0,0].scatter(AX, AY, s=0.3 ,c='grey',label='Upsampled grid')
    ax[0,0].legend(loc='lower right')
    ax[0,0].scatter(AX, AY, c=A1,cmap=cmap,vmin=p1min, vmax=p1max)
    ax[0,1].contourf(AX, AY, A1, cmap=cmap, vmin=p1min, vmax=p1max,levels=20)

    a = ax[1,0].scatter(BX, BY, s=1.5, c=B1,cmap=cmap,vmin=p1min, vmax=p1max)
    ax[1,1].contourf(BX, BY, B1,cmap=cmap, vmin=p1min, vmax=p1max, levels=20)
    fig.suptitle('Along strike slip')
    fig.tight_layout()
    fig.subplots_adjust(right=0.93)
    cbar_ax = fig.add_axes([0.95, 0.2, 0.02, 0.6])
    fig.colorbar(a, cax=cbar_ax,label='slip (m)')
    if savefig is not None:
        fig.savefig(savefig+'_alongstrike.png')

    fig, ax = plt.subplots(nrows=2,ncols=2,figsize=(9,9))
    ax[0,0].scatter(BX, BY, s=0.2 ,c='grey')
    ax[0,0].scatter(AX, AY, c=A2, cmap=cmap, vmin=0, vmax=p2max)
    ax[0,1].contourf(AX, AY, A2, cmap=cmap, vmin=0, vmax=p2max,levels=20)

    a = ax[1,0].scatter(BX, BY, s=1.5, c=B2, cmap=cmap,vmin=0, vmax=p2max)
    ax[1,1].contourf(BX, BY, B2, cmap=cmap,vmin=0, vmax=p2max, levels=20)
    fig.suptitle('Along dip slip')
    fig.tight_layout()
    fig.subplots_adjust(right=0.93)
    cbar_ax = fig.add_axes([0.95, 0.2, 0.02, 0.6])
    fig.colorbar(a, cax=cbar_ax,label='slip (m)')
    
    if savefig is not None:
        fig.savefig(savefig+'_alongdip.png')
        
def plot_netcdf(fault, savefig=None):
    AX, AY = fault.Xg, fault.Yg
    A1, A2 = fault.gSls , fault.gSld
    cmap = cm.get_cmap('viridis')
    
    fig, ax = plt.subplots(nrows=2,ncols=2,figsize=(9,9))
    ax[0,0].scatter(AX, AY, s=0.3 ,c='grey',label='Upsampled grid')
    ax[0,0].legend(loc='lower right')
    ax[0,0].scatter(AX, AY, c=A1,cmap=cmap,vmin=A1.min(), vmax=A1.max())
    c1 =ax[0,1].contourf(AX, AY, A1, cmap=cmap, vmin=A1.min(), vmax=A1.max(),levels=12)
    ax[1,0].scatter(AX, AY, s=1.5, c=A2, cmap=cmap,vmin=0, vmax=A2.max())
    c2 = ax[1,1].contourf(AX, AY, A2, cmap=cmap,vmin=0, vmax=A2.max(), levels=12)

    fig.suptitle('NetCDF ')
    fig.tight_layout()
    fig.subplots_adjust(right=0.93)
    cbar_ax = fig.add_axes([0.95, 0.6, 0.02, 0.3])
    fig.colorbar(c1, cax=cbar_ax,label='slip (m)')
    cbar_ax = fig.add_axes([0.95, 0.1, 0.02, 0.3])
    fig.colorbar(c2, cax=cbar_ax,label='slip (m)')
        
    if savefig is not None:
        fig.savefig(savefig)
        
def plot_tricontourf(Data, cm_n=None):
    """tricontourf is for non-gridded data
    contourf is for gridded data (z has to be 2d Array)"""
    if cm_n is not None:
        cmap = cm.get_cmap(cm_n)
    else:
        cmap = cm.get_cmap('viridis')
        
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(5,5))
    c1 =ax.tricontourf(Data[:,0], Data[:,1], Data[:,2], cmap=cmap,levels=12)
    fig.tight_layout()
    fig.subplots_adjust(right=0.93)
    cbar_ax = fig.add_axes([0.95, 0.2, 0.02, 0.5])
    fig.colorbar(c1, cax=cbar_ax)
