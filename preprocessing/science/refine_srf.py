import numpy as np
import pyproj
import matplotlib.pyplot as plt
import scipy.ndimage
from scipy import interpolate
from math import floor
import os
import argparse


class FaultPlane:
    def __init__(self):
        self.nx = 0
        self.ny = 0
        self.ndt = 0
        self.PSarea = 0
        self.dt = 0
        # array member initialized to dummy value
        self.lon = 0
        self.lat = 0
        self.x = 0
        self.y = 0
        self.depth = 0
        self.t0 = 0
        self.slip1 = 0
        self.strike = 0
        self.dip = 0
        self.rake = 0
        self.aSR = 0
        self.myt = 0

    def init_spatial_arrays(self, nx, ny):
        self.nx = nx
        self.ny = ny
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

    def init_aSR(self):
        self.aSR = np.zeros((self.ny, self.nx, self.ndt))

    def extend_aSR(self, ndt_old, ndt_new):
        tmpSR = np.copy(self.aSR)
        self.ndt = ndt_new
        self.aSR = np.zeros((self.ny, self.nx, self.ndt))
        self.aSR[:, :, 0:ndt_old] = tmpSR[:, :, :]

    def compute_xy_from_latlon(self, proj_string):
        if args.proj:
            lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")
            myproj = pyproj.Proj(proj_string)
            self.x, self.y = pyproj.transform(lla, myproj, self.lon, self.lat)
        else:
            print("no proj string specified!")
            self.x, self.y = self.lon, self.lat

    def compute_latlon_from_xy(self, proj_string):
        if args.proj:
            lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")
            myproj = pyproj.Proj(proj_string)
            self.lon, self.lat = pyproj.transform(myproj, lla, self.x, self.y)
        else:
            self.lon, self.lat = self.x, self.y

    def compute_time(self):
        self.myt = np.linspace(0, (self.ndt - 1) * self.dt, self.ndt)

    def write_srf(self, fname):
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

    def init_from_srf(self, fname):
        with open(fname) as fid:
            # version
            line = fid.readline()
            version = float(line)
            if not (abs(version - 1.0) < 1e-03 or abs(version - 2.0) < 1e-03):
                print("srf version: %s not supported" % (line))
                raise
            # skip comments
            while True:
                line = fid.readline()
                if not line.startswith("#"):
                    break
            line_el = line.split()
            if line_el[0] != "PLANE":
                print("no plane specified")
                raise
            if line_el[1] != "1":
                print("only one plane supported")
                raise
            line_el = fid.readline().split()
            nx = int(line_el[2])
            ny = int(line_el[3])
            line_el = fid.readline().split()
            # check that the plane data are consistent with the number of points
            assert int(line_el[1]) == nx * ny
            self.init_spatial_arrays(nx, ny)
            for j in range(ny):
                for i in range(nx):
                    # first header line
                    line = fid.readline()
                    if abs(version - 1.0) < 1e-03:
                        self.lon[j, i], self.lat[j, i], self.depth[j, i], self.strike[j, i], self.dip[j, i], self.PSarea, self.t0[j, i], dt = [float(v) for v in line.split()]
                    else:
                        # srf version 2 has also density and Vs in that line
                        self.lon[j, i], self.lat[j, i], self.depth[j, i], self.strike[j, i], self.dip[j, i], self.PSarea, self.t0[j, i], dt = [float(v) for v in line.split()[0:-2]]
                    # second header line
                    line = fid.readline()
                    self.rake[j, i], self.slip1[j, i], ndt1, slip2, ndt2, slip3, ndt3 = [float(v) for v in line.split()]
                    if max(slip2, slip3) > 0.0:
                        print("this script assumes slip2 and slip3 are zero", slip2, slip3)
                        raise
                    ndt1 = int(ndt1)
                    if max(i, j) == 0:
                        self.ndt = ndt1
                        self.dt = dt
                        self.init_aSR()
                    lSTF = []
                    if ndt1 == 0:
                        continue
                    if ndt1 > self.ndt:
                        print("a larger ndt (%d> %d)was found for point source (i,j) = (%d, %d) extending aSR array..." % (ndt1, self.ndt, i, j))
                        self.extend_aSR(self.ndt, ndt1)
                    if abs(dt - self.dt) > 1e-6:
                        print("this script assumes that dt is the same for all sources", dt, self.dt)
                        raise
                    while True:
                        line = fid.readline()
                        lSTF.extend(line.split())
                        if len(lSTF) == ndt1:
                            self.aSR[j, i, 0:ndt1] = np.array([float(v) for v in lSTF])
                            break

    def upsample_fault(self, spatial_order, spatial_zoom, temporal_zoom):
        # time vector
        ndt2 = (self.ndt - 1) * temporal_zoom + 1
        ny2, nx2 = self.ny * spatial_zoom, self.nx * spatial_zoom
        # resampled source
        pf = FaultPlane()
        pf.init_spatial_arrays(nx2, ny2)
        pf.ndt = ndt2
        pf.init_aSR()

        pf.dt = self.dt / temporal_zoom
        pf.compute_time()

        # upsample spatially all these quantitites
        allarr = np.array([self.x, self.y, self.depth, self.t0, self.slip1, self.strike, self.dip, self.rake])
        nd = allarr.shape[0]
        allarr0 = np.zeros((nd, pf.ny, pf.nx))

        for k in range(nd):
            allarr0[k, :, :] = scipy.ndimage.zoom(allarr[k, :, :], spatial_zoom, order=spatial_order)
        pf.x, pf.y, pf.depth, pf.t0, pf.slip1, pf.strike, pf.dip, pf.rake = allarr0
        pf.compute_latlon_from_xy(args.proj)
        pf.PSarea = self.PSarea / spatial_zoom ** 2

        aSRa = np.zeros((pf.ny, pf.nx, self.ndt))
        for k in range(self.ndt):
            aSRa[:, :, k] = scipy.ndimage.zoom(self.aSR[:, :, k], spatial_zoom, order=spatial_order)

        # interpolate temporally the AST
        for j in range(pf.ny):
            for i in range(pf.nx):
                f = interpolate.interp1d(self.myt, aSRa[j, i, :], kind="quadratic")
                pf.aSR[j, i, :] = f(pf.myt)
                # With a cubic interpolation, the interpolated slip1 may be negative which does not make sense.
                if pf.slip1[j, i] < 0:
                    pf.aSR[j, i, :] = 0
                    continue
                # should be the SR
                integral_STF = np.trapz(np.abs(pf.aSR[j, i, :]), dx=pf.dt)
                if abs(integral_STF) > 0:
                    pf.aSR[j, i, :] = pf.slip1[j, i] * pf.aSR[j, i, :] / integral_STF
        return pf


parser = argparse.ArgumentParser(description="upsample temporally and spatially a kinematic model (should be a planar model) in the standard rupture format")
parser.add_argument("filename", help="filename of the srf file")
parser.add_argument("--proj", help="proj4 string (might be better to upsample the geometry in the local coordinate system)")
parser.add_argument("--spatial_order", nargs=1, metavar=("spatial_order"), default=([3]), help="spatial order of the interpolation", type=int)
parser.add_argument("--spatial_zoom", nargs=1, metavar=("spatial_zoom"), required=True, help="level of spatial upsampling", type=int)
parser.add_argument("--temporal_zoom", nargs=1, metavar=("temporal_zoom"), required=True, help="level of temporal upsampling", type=int)
args = parser.parse_args()

p1 = FaultPlane()
p1.init_from_srf(args.filename)
p1.compute_xy_from_latlon(args.proj)
p1.compute_time()

p2 = p1.upsample_fault(spatial_order=args.spatial_order[0], spatial_zoom=args.spatial_zoom[0], temporal_zoom=args.temporal_zoom[0])
prefix, ext = os.path.splitext(args.filename)
fnout = prefix + "_resampled" + ".srf"
p2.write_srf(fnout)
