#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
from seissolxdmf import seissolxdmf as sx


def main():
    description = '''Plots a map view of the ground motion intensity.'''
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        'fin', help='XDMF input file name.')
    parser.add_argument(
        'fout', help='Output image file name.')
    parser.add_argument(
        'dataid', help='Which intensity metric to use', default='PGV')
    parser.add_argument(
        '--cmap', help='Matplotlib colormap', default='viridis')
    parser.add_argument(
        '--dpi', default=200)
    parser.add_argument(
        '--log', help='Plot log of intensity', action='store_true',
        default=False)
    parser.add_argument(
        '--bounds', nargs='+', help='xmin xmax ymin ymax', type=float)
    parser.add_argument(
        '--contour', action='store_true', default=False)
    args = parser.parse_args()

    ds = sx(args.fin)
    geo = ds.ReadGeometry()
    connect = ds.ReadConnect()
    coords = geo[connect].mean(axis=1)
    x, y = coords.T[0], coords.T[1]
    imt = ds.ReadData(args.dataid)
    clabel = args.dataid
    if args.log:
        imt = np.log(imt)
        clabel = 'log(%s)' % clabel

    if args.bounds:
        idx = ((x > args.bounds[0]) & (x < args.bounds[1]) &
               (y > args.bounds[2]) & (y < args.bounds[3]))
        vmin = imt[idx].min()
        vmax = imt[idx].max()
    else:
        vmin, vmax = None, None

    plt.figure(figsize=(6, 5))

    if args.contour:
        plt.tricontourf(x, y, imt, vmin=vmin, vmax=vmax)
    else:
        plt.tripcolor(geo.T[0], geo.T[1], connect, imt, vmin=vmin, vmax=vmax)
    if args.bounds:
        plt.xlim(args.bounds[0], args.bounds[1])
        plt.ylim(args.bounds[2], args.bounds[3])
    plt.colorbar(label=clabel)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.savefig(args.fout, dpi=int(args.dpi))
    print('Image file saved to %s' % args.fout)
    plt.close('all')


if __name__ == '__main__':
    main()
