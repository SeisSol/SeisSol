#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from seissolxdmf import seissolxdmf as sx

plt.rcParams['savefig.facecolor'] = 'w'


def main():
    description = '''Creates an animation of the velocity/displacement
    wavefield from SeisSol free-surface XDMF output file.'''
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        'fin', help='SeisSol surface (XDMF) output file.')
    parser.add_argument(
        'fout', help='Output animation file.', default='animation.gif')
    parser.add_argument(
        'dataid', help='Which data component to use', default='v1')
    parser.add_argument(
        '--cmap', help='Matplotlib colormap', default='bwr')
    parser.add_argument(
        '--dpi', default=200, type=int)
    parser.add_argument(
        '--bounds', nargs='+', help='xmin xmax ymin ymax', type=float)
    args = parser.parse_args()

    ds = sx(args.fin)
    data = ds.ReadData(args.dataid)
    geo = ds.ReadGeometry()
    connect = ds.ReadConnect()
    dt = ds.ReadTimeStep()
    nt = data.shape[0]
    tarr = np.arange(0.0, nt * dt, dt)

    fig = plt.figure(figsize=(8, 7))
    ax = plt.subplot(111)
    dmax = abs(data).max()
    im = ax.tripcolor(
        geo.T[0], geo.T[1], connect, data[0], cmap=args.cmap, vmin=-dmax,
        vmax=dmax)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_xlim(args.bounds[0], args.bounds[1])
    ax.set_ylim(args.bounds[2], args.bounds[3])
    fig.colorbar(im, label=args.dataid)

    def animate(i):
        im.set_array(data[i])
        ax.set_title('$t=$%.1f' % tarr[i])

    def progress_callback(i, n):
        return print(f'Saving frame {i} of {n}')

    anim = animation.FuncAnimation(
        fig, animate, frames=data.shape[0],
        interval=dt * 1000)
    anim.save(args.fout, dpi=args.dpi, progress_callback=progress_callback)
    plt.close('all')


if __name__ == '__main__':
    main()
