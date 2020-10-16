import numpy as np
import pyproj
import matplotlib.pyplot as plt
import scipy.ndimage
from scipy import interpolate
from math import floor
import os


class FaultPlane:
    def __init__(self, nx, ny):
        self.nx = nx
        self.ny = ny
        self.ndt = 0
        self.lon = np.zeros((ny, nx))
        self.lat = np.zeros((ny, nx))
        self.x = np.zeros((ny, nx))
        self.y = np.zeros((ny, nx))
        self.depth = np.zeros((ny, nx))
        self.t0 = np.zeros((ny, nx))
        self.slip1 = np.zeros((ny, nx))
        self.strike = np.zeros((ny, nx))
        self.dip = np.zeros((ny, nx))
        self.rake = np.zeros((ny, nx))
        self.aSR = 0
        self.PSarea = 0
        self.dt = 0
        self.myt = 0

    def init_aSR(self):
        self.aSR = np.zeros((self.ny, self.nx, self.ndt))

    def compute_xy_from_latlon(self, proj_string):
        if args.proj != None:
            lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")
            myproj = pyproj.Proj(proj_string)
            self.x, self.y = pyproj.transform(lla, myproj, self.lon, self.lat)
        else:
            print("no proj string specified!")
            self.x, self.y = self.lon, self.lat

    def compute_latlon_from_xy(self, proj_string):
        if args.proj != None:
            lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")
            myproj = pyproj.Proj(proj_string)
            self.lon, self.lat = pyproj.transform(myproj, lla, self.x, self.y)
        else:
            self.lon, self.lat = self.x, self.y

    def compute_time(self):
        self.myt = np.linspace(0, (self.ndt - 1) * self.dt, self.ndt)

    def write_SRF(self, fname):
        with open(fname, "w") as fout:
            fout.write("1.0\n")
            fout.write("POINTS %d\n" % (self.nx * self.ny))
            for j in range(self.ny):
                for i in range(self.nx):
                    fout.write("%g %g %g %g %g %e %g %g\n" % (self.lon[j, i], self.lat[j, i], self.depth[j, i], self.strike[j, i], self.dip[j, i], self.PSarea, self.t0[j, i], self.dt))
                    fout.write("%g %g %d %f %d %f %d\n" % (self.rake[j, i], self.slip1[j, i], self.ndt, 0.0, 0, 0.0, 0))
                    np.savetxt(fout, self.aSR[j, i, :], fmt="%g", newline=" ")
                    fout.write("\n")
        print("done writing", fname)


def ReadSRF(fname):
    with open(fname) as fid:
        # version
        line = fid.readline()
        version = float(line)
        if abs(version - 1.0) > 1e-03:
            print("SRF version: %s not supported" % (line))
            raise
        line_el = fid.readline().split()
        if line_el[0] != "PLANE":
            print("no plane specified")
            raise
        line_el = fid.readline().split()
        nx = int(line_el[2])
        ny = int(line_el[3])
        line_el = fid.readline().split()
        if line_el[0] != "POINTS":
            print("only one plane supported")
            raise
        # check that the plane data are consistent with the number of points
        assert int(line_el[1]) == nx * ny
        p1 = FaultPlane(nx, ny)
        for j in range(ny):
            for i in range(nx):
                # first header line
                line = fid.readline()
                p1.lon[j, i], p1.lat[j, i], p1.depth[j, i], p1.strike[j, i], p1.dip[j, i], p1.PSarea, p1.t0[j, i], p1.dt = [float(v) for v in line.split()]
                # second header line
                line = fid.readline()
                p1.rake[j, i], p1.slip1[j, i], ndt1, slip2, ndt2, slip3, ndt3 = [float(v) for v in line.split()]
                ndt1 = int(ndt1)
                if max(i, j) == 0:
                    p1.ndt = ndt1
                    p1.init_aSR()
                lSTF = []
                if ndt1 == 0:
                    continue
                if ndt1 > p1.ndt:
                    print("this script assumes that ndt is not larger than the first ndt in the file", ndt1, p1.ndt)
                    raise
                while True:
                    line = fid.readline()
                    lSTF.extend(line.split())
                    if len(lSTF) == ndt1:
                        p1.aSR[j, i, 0:ndt1] = np.array([float(v) for v in lSTF])
                        break
    return p1


def upsampleFault(p, spatial_order, spatial_zoom, temporal_zoom):
    # time vector
    ndt2 = (p.ndt - 1) * temporal_zoom + 1
    ny2, nx2 = p.ny * spatial_zoom, p.nx * spatial_zoom
    # resampled source
    p2 = FaultPlane(nx2, ny2)
    p2.ndt = ndt2
    p2.init_aSR()

    p2.dt = p.dt / temporal_zoom
    p2.compute_time()

    # upsample spatially all these quantitites
    allarr = np.array([p.x, p.y, p.depth, p.t0, p.slip1, p.strike, p.dip, p.rake])
    nd = allarr.shape[0]
    allarr0 = np.zeros((nd, p2.ny, p2.nx))

    for k in range(nd):
        allarr0[k, :, :] = scipy.ndimage.zoom(allarr[k, :, :], spatial_zoom, order=spatial_order)
    p2.x, p2.y, p2.depth, p2.t0, p2.slip1, p2.strike, p2.dip, p2.rake = allarr0
    p2.compute_latlon_from_xy(args.proj)
    p2.PSarea = p.PSarea / spatial_zoom ** 2

    aSRa = np.zeros((p2.ny, p2.nx, p.ndt))
    for k in range(p.ndt):
        aSRa[:, :, k] = scipy.ndimage.zoom(p.aSR[:, :, k], spatial_zoom, order=spatial_order)

    # interpolate temporally the AST
    for j in range(p2.ny):
        for i in range(p2.nx):
            f = interpolate.interp1d(p.myt, aSRa[j, i, :], kind="quadratic")
            p2.aSR[j, i, :] = f(p2.myt)
            if p2.slip1[j, i] < 0:
                p2.aSR[j, i, :] = 0
                continue
            # should be the SR
            integral_STF = np.trapz(np.abs(p2.aSR[j, i, :]), dx=p2.dt)
            if abs(integral_STF) > 0:
                p2.aSR[j, i, :] = p2.slip1[j, i] * p2.aSR[j, i, :] / integral_STF
    return p2


import argparse

parser = argparse.ArgumentParser(description="upsample temporally and spatially a kinematic model (should be a planar model) in the standard rupture format")
parser.add_argument("filename", help="filename of the srf file")
parser.add_argument("--proj", help="proj4 string (might be better to upsample the geometry in the local coordinate system)")
parser.add_argument("--spatial_order", nargs=1, metavar=("spatial_order"), default=([3]), help="spatial order of the interpolation", type=int)
parser.add_argument("--spatial_zoom", nargs=1, metavar=("spatial_zoom"), default=([3]), help="level of spatial upsampling", type=int)
parser.add_argument("--temporal_zoom", nargs=1, metavar=("temporal_zoom"), default=([10]), help="level of temporal upsampling", type=int)
args = parser.parse_args()


p1 = ReadSRF(args.filename)
p1.compute_xy_from_latlon(args.proj)
p1.compute_time()
p2 = upsampleFault(p1, spatial_order=args.spatial_order[0], spatial_zoom=args.spatial_zoom[0], temporal_zoom=args.temporal_zoom[0])
prefix, ext = os.path.splitext(args.filename)
fnout = prefix + "_resampled" + ".srf"
p2.write_SRF(fnout)
