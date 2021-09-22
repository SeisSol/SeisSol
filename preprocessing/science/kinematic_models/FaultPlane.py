import numpy as np
import pyproj
import scipy.ndimage
from scipy import interpolate
from netCDF4 import Dataset
from Yoffe import regularizedYoffe

def writeNetcdf(sname, lDimVar, lName, lData, paraview_readable=False):
    """create a netcdf file either readable by ASAGI
    or by paraview (paraview_readable=True)

    Parameters
    ----------
    sname: str prefix name of output file
    lDimVar: list of 1d numpy array containing the dimension variables
    lName: list if str containing the name of the nd variables
    lData: list if n-d numpy array containing the data
    paraview_readable: bool
    """
    fname = f"{sname}.nc"
    print("writing " + fname)

    with Dataset(fname, "w", format="NETCDF4") as rootgrp:
        #Create dimension and 1d variables
        sdimVarNames='uvwxyz'
        dims = []
        for i, xi in enumerate(lDimVar):
            nxi = xi.shape[0]
            dimName=sdimVarNames[i]
            dims.append(dimName)
            rootgrp.createDimension(dimName, nxi)
            vx = rootgrp.createVariable(dimName, "f4", (dimName,))
            vx[:] = xi
        dims.reverse()
        dims = tuple(dims)

        if paraview_readable:
            for i in range(len(lName)):
                vTd = rootgrp.createVariable(lName[i], "f4", dims)
                vTd[:] = lData[i]
        else:
            ldata4 = [(name, "f4") for name in lName]
            ldata8 = [(name, "f8") for name in lName]
            mattype4 = np.dtype(ldata4)
            mattype8 = np.dtype(ldata8)
            mat_t = rootgrp.createCompoundType(mattype4, "material")

            # this transform the nD array into an array of tuples
            arr = np.stack([lData[i] for i in range(len(lName))], axis=len(dims))
            newarr = arr.view(dtype=mattype8)
            newarr = newarr.reshape(newarr.shape[:-1])
            mat = rootgrp.createVariable("data", mat_t, dims)
            mat[:] = newarr


def interpolate_nan_from_neighbors(array):
    """rise_time and tacc may not be defined where there is no slip (no SR function).
    in this case, we interpolate from neighbors
    source: https://stackoverflow.com/questions/37662180/interpolate-missing-values-2d-python
    """
    x = np.arange(0, array.shape[1])
    y = np.arange(0, array.shape[0])
    # mask invalid values
    array = np.ma.masked_invalid(array)
    xx, yy = np.meshgrid(x, y)
    # get only the valid values
    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]
    return interpolate.griddata((x1, y1), newarr.ravel(), (xx, yy), method="linear", fill_value=np.average(array))


def upsample_quantities(allarr, spatial_order, spatial_zoom, padding="constant", extra_padding_layer=False):
    """1. pad
    2. upsample, adding spatial_zoom per node
    """
    nd = allarr.shape[0]
    ny, nx = [val * spatial_zoom for val in allarr[0].shape]
    if extra_padding_layer:
        # required for vertex aligned netcdf format
        nx = nx + 2
        ny = ny + 2
    allarr0 = np.zeros((nd, ny, nx))
    for k in range(nd):
        if padding == "extrapolate":
            my_array = np.pad(allarr[k, :, :], ((1, 1), (1, 1)), "reflect", reflect_type="odd")
        else:
            my_array = np.pad(allarr[k, :, :], ((1, 1), (1, 1)), padding)
        if extra_padding_layer:
            ncrop = spatial_zoom - 1
        else:
            ncrop = spatial_zoom
        my_array = scipy.ndimage.zoom(my_array, spatial_zoom, order=spatial_order, mode="grid-constant", grid_mode=True)
        if ncrop > 0:
            allarr0[k, :, :] = my_array[ncrop:-ncrop, ncrop:-ncrop]

    return allarr0


class FaultPlane:
    def __init__(self):
        self.nx = 0
        self.ny = 0
        self.ndt = 0
        self.PSarea_cm2 = 0
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
        "extend aSR array to more time samplings"
        tmpSR = np.copy(self.aSR)
        self.ndt = ndt_new
        self.aSR = np.zeros((self.ny, self.nx, self.ndt))
        self.aSR[:, :, 0:ndt_old] = tmpSR[:, :, :]

    def compute_xy_from_latlon(self, proj):
        if proj:
            from pyproj import Transformer

            transformer = Transformer.from_crs("epsg:4326", proj[0], always_xy=True)
            self.x, self.y = transformer.transform(self.lon, self.lat)
        else:
            print("no proj string specified!")
            self.x, self.y = self.lon, self.lat

    def compute_latlon_from_xy(self, proj):
        if proj:
            from pyproj import Transformer

            transformer = Transformer.from_crs(proj[0], "epsg:4326", always_xy=True)
            self.lon, self.lat = transformer.transform(self.x, self.y)
        else:
            self.lon, self.lat = self.x, self.y

    def compute_time_array(self):
        self.myt = np.linspace(0, (self.ndt - 1) * self.dt, self.ndt)

    def write_srf(self, fname):
        "write kinematic model to a srf file (standard rutpure format)"
        with open(fname, "w") as fout:
            fout.write("1.0\n")
            fout.write("POINTS %d\n" % (self.nx * self.ny))
            for j in range(self.ny):
                for i in range(self.nx):
                    fout.write("%g %g %g %g %g %e %g %g\n" % (self.lon[j, i], self.lat[j, i], self.depth[j, i], self.strike[j, i], self.dip[j, i], self.PSarea_cm2, self.t0[j, i], self.dt))
                    fout.write("%g %g %d %f %d %f %d\n" % (self.rake[j, i], self.slip1[j, i], self.ndt, 0.0, 0, 0.0, 0))
                    np.savetxt(fout, self.aSR[j, i, :], fmt="%g", newline=" ")
                    fout.write("\n")
        print("done writing", fname)

    def init_from_srf(self, fname):
        "init object by reading a srf file (standard rutpure format)"
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
                raise NotImplementedError
            line_el = fid.readline().split()
            nx, ny = [int(val) for val in line_el[2:4]]
            line_el = fid.readline().split()
            # check that the plane data are consistent with the number of points
            assert int(line_el[1]) == nx * ny
            self.init_spatial_arrays(nx, ny)
            for j in range(ny):
                for i in range(nx):
                    # first header line
                    line = fid.readline()
                    # rho_vs are only present for srf version 2
                    self.lon[j, i], self.lat[j, i], self.depth[j, i], self.strike[j, i], self.dip[j, i], self.PSarea_cm2, self.t0[j, i], dt, *rho_vs = [float(v) for v in line.split()]
                    # second header line
                    line = fid.readline()
                    self.rake[j, i], self.slip1[j, i], ndt1, slip2, ndt2, slip3, ndt3 = [float(v) for v in line.split()]
                    if max(slip2, slip3) > 0.0:
                        print("this script assumes slip2 and slip3 are zero", slip2, slip3)
                        raise NotImplementedError
                    ndt1 = int(ndt1)
                    if max(i, j) == 0:
                        self.ndt = ndt1
                        self.dt = dt
                        self.init_aSR()
                    lSTF = []
                    if ndt1 == 0:
                        continue
                    if ndt1 > self.ndt:
                        print(f"a larger ndt ({ndt1}> {self.ndt}) was found for point source (i,j) = ({i}, {j}) extending aSR array...")
                        self.extend_aSR(self.ndt, ndt1)
                    if abs(dt - self.dt) > 1e-6:
                        print("this script assumes that dt is the same for all sources", dt, self.dt)
                        raise NotImplementedError
                    while True:
                        line = fid.readline()
                        lSTF.extend(line.split())
                        if len(lSTF) == ndt1:
                            self.aSR[j, i, 0:ndt1] = np.array([float(v) for v in lSTF])
                            break

    def assess_Yoffe_parameters(self):
        "compute rise_time (slip duration) and t_acc (peak SR) from SR time histories"
        self.rise_time = np.zeros((self.ny, self.nx))
        self.tacc = np.zeros((self.ny, self.nx))
        for j in range(self.ny):
            for i in range(self.nx):
                if not self.slip1[j, i]:
                    self.rise_time[j, i] = np.nan
                    self.tacc[j, i] = np.nan
                else:
                    first_non_zero = np.amin(np.where(self.aSR[j, i, :])[0])
                    last_non_zero = np.amax(np.where(self.aSR[j, i, :])[0])
                    id_max = np.where(self.aSR[j, i, :] == np.amax(self.aSR[j, i, :]))[0]
                    self.rise_time[j, i] = (last_non_zero - first_non_zero + 1) * self.dt
                    self.tacc[j, i] = (id_max - first_non_zero + 1) * self.dt
                    self.t0[j, i] += first_non_zero * self.dt
        self.rise_time = interpolate_nan_from_neighbors(self.rise_time)
        self.tacc = interpolate_nan_from_neighbors(self.tacc)

        print("slip rise_time (min, 50%, max)", np.amin(self.rise_time), np.median(self.rise_time), np.amax(self.rise_time))
        print("tacc (min, 50%, max)", np.amin(self.tacc), np.median(self.tacc), np.amax(self.tacc))

    def upsample_fault(self, spatial_order, spatial_zoom, temporal_zoom, proj, use_Yoffe=False):
        "increase spatial and temporal resolution of kinematic model by interpolation"
        # time vector
        ndt2 = (self.ndt - 1) * temporal_zoom + 1
        ny2, nx2 = self.ny * spatial_zoom, self.nx * spatial_zoom
        # resampled source
        pf = FaultPlane()
        pf.init_spatial_arrays(nx2, ny2)
        pf.ndt = ndt2
        pf.init_aSR()

        pf.dt = self.dt / temporal_zoom
        pf.compute_time_array()

        # upsample spatially geometry (bilinear interpolation)
        allarr = np.array([self.x, self.y, self.depth])
        pf.x, pf.y, pf.depth = upsample_quantities(allarr, spatial_order=1, spatial_zoom=spatial_zoom, padding="extrapolate")

        # upsample other quantities
        allarr = np.array([self.t0, self.strike, self.dip, self.rake])
        pf.t0, pf.strike, pf.dip, pf.rake = upsample_quantities(allarr, spatial_order, spatial_zoom, padding="edge")
        # the interpolation may generate some acausality that we here prevent
        pf.t0 = np.maximum(pf.t0, np.amin(self.t0))

        allarr = np.array([self.slip1])
        (pf.slip1,) = upsample_quantities(allarr, spatial_order, spatial_zoom, padding="constant")
        pf.compute_latlon_from_xy(proj)
        pf.PSarea_cm2 = self.PSarea_cm2 / spatial_zoom ** 2
        ratio_potency = np.sum(pf.slip1) * pf.PSarea_cm2 / (np.sum(self.slip1) * self.PSarea_cm2)
        print(f"seismic potency ratio (upscaled over initial): {ratio_potency}")

        if use_Yoffe:
            self.assess_Yoffe_parameters()
            allarr = np.array([self.rise_time, self.tacc])
            pf.rise_time, pf.tacc = upsample_quantities(allarr, spatial_order, spatial_zoom, padding="edge")
            pf.rise_time = np.maximum(pf.rise_time, np.amin(self.rise_time))
            pf.tacc = np.maximum(pf.tacc, np.amin(self.tacc))
            print('using ts = tacc / 1.27 to compute the regularized Yoffe')
            ts = pf.tacc / 1.27
            tr = pf.rise_time - 2.0 * ts
            for j in range(pf.ny):
                for i in range(pf.nx):
                    for k, tk in enumerate(pf.myt):
                        pf.aSR[j, i, k] = pf.slip1[j, i] * regularizedYoffe(tk, ts[j, i], tr[j, i])
        else:
            aSRa = np.zeros((pf.ny, pf.nx, self.ndt))
            for k in range(self.ndt):
                aSRa[:, :, k] = upsample_quantities(np.array([self.aSR[:, :, k]]), spatial_order, spatial_zoom, padding="constant")

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

    def compute_corrected_slip_for_differing_area(self, proj):
        """
        self.PSarea_cm2 may slightly differ from the patch area from the fault geometry
        (e.g. due to the projection)
        Therefore, we need to update slip to keep seismic potency (area*slip) unchanged
        """
        cm2m = 0.01
        km2m = 1e3
        PSarea_m2 = self.PSarea_cm2 * cm2m * cm2m
        self.compute_xy_from_latlon(proj)
        nx, ny = self.nx, self.ny
        # Compute actual dx and dy from coordinates
        dy = np.zeros((ny, nx))
        dx = np.zeros((ny, nx))
        # central difference for the inside
        coords = np.array((self.x, self.y, -km2m * self.depth))
        for i in range(0, nx):
            p0 = coords[:, 0 : ny - 2, i] - coords[:, 2:ny, i]
            dy[1 : ny - 1, i] = 0.5 * np.linalg.norm(p0, axis=0)
        # special case of 0 and ny-1
        p0 = coords[:, 1, :] - coords[:, 0, :]
        dy[0, :] = np.linalg.norm(p0, axis=0)
        p0 = coords[:, ny - 1, :] - coords[:, ny - 2, :]
        dy[ny - 1, :] = np.linalg.norm(p0, axis=0)
        # dx for coordinates
        for j in range(0, ny):
            p0 = coords[:, j, 0 : nx - 2] - coords[:, j, 2:nx]
            dx[j, 1 : nx - 1] = 0.5 * np.linalg.norm(p0, axis=0)
        p0 = coords[:, :, 1] - coords[:, :, 0]
        dx[:, 0] = np.linalg.norm(p0, axis=0)
        p0 = coords[:, :, nx - 1] - coords[:, :, nx - 2]
        dx[:, nx - 1] = np.linalg.norm(p0, axis=0)
        factor_area = dx[:, :] * dy[:, :] / PSarea_m2
        slip1 = self.slip1 * factor_area
        print(
            f"done correcting slip for area. \
The correcting factor ranges between {np.amin(factor_area)} and {np.amax(factor_area)}"
        )
        return slip1

    def generate_netcdf_fl33(self, prefix, spatial_order, spatial_zoom, proj, write_paraview):
        "generate netcdf files to be used with SeisSol friction law 33"

        cm2m = 0.01
        km2m = 1e3
        # a kinematic model defines the fault quantities at the subfault center
        # a netcdf file defines the quantities at the nodes
        # therefore the extra_padding_layer=True, and the added di below
        cslip = self.compute_corrected_slip_for_differing_area(proj)
        (slip,) = upsample_quantities(np.array([cslip]), spatial_order, spatial_zoom, padding="constant", extra_padding_layer=True)
        allarr = np.array([self.t0, self.rake, self.rise_time, self.tacc])
        rupttime, rake, rise_time, tacc = upsample_quantities(allarr, spatial_order, spatial_zoom, padding="edge", extra_padding_layer=True)
        # upsampled duration, rise_time and acc_time may not be smaller than initial values
        # at least rise_time could lead to a non-causal kinematic model
        rupttime = np.maximum(rupttime, np.amin(self.t0))
        rise_time = np.maximum(rise_time, np.amin(self.rise_time))
        tacc = np.maximum(tacc, np.amin(self.tacc))

        rake_rad = np.radians(rake)
        strike_slip = slip * np.cos(rake_rad) * cm2m
        dip_slip = slip * np.sin(rake_rad) * cm2m

        ny, nx = slip.shape
        dx = np.sqrt(self.PSarea_cm2 * cm2m * cm2m)
        ldataName = ["strike_slip", "dip_slip", "rupture_onset", "effective_rise_time", "acc_time"]
        lgridded_myData = [strike_slip, dip_slip, rupttime, rise_time, tacc]
        # we could do directly
        # di = 1.0 / (2 * spatial_zoom)
        # xb = np.linspace(-di, self.nx + di, nx) * dx
        # yb = np.linspace(-di, self.ny + di, ny) * dx
        # But we compute xb and yb based on the upsampled coordinates, to account for possible warping due to the projection
        allarr = np.array([self.x, self.y, -km2m * self.depth])
        coords = upsample_quantities(allarr, spatial_order=1, spatial_zoom=spatial_zoom, padding="extrapolate", extra_padding_layer=True)

        p0 = coords[:, (ny - 1) // 2, :] - coords[:, (ny - 1) // 2 - 1, :]
        dx1 = np.linalg.norm(p0, axis=0)
        xb = np.cumsum(dx1) - 1.5 * dx1[0]

        p0 = coords[:, :, (nx - 1) // 2] - coords[:, :, (nx - 1) // 2 - 1]
        dy1 = np.linalg.norm(p0, axis=0)
        yb = np.cumsum(dy1) - 1.5 * dy1[0]

        prefix2 = f"{prefix}_{spatial_zoom}_o{spatial_order}"
        if write_paraview:
            # see comment above
            for i, sdata in enumerate(ldataName):
                writeNetcdf(prefix2 + sdata, [xb, yb], [sdata], [lgridded_myData[i]], paraview_readable=True)
        writeNetcdf(prefix2, [xb, yb], ldataName, lgridded_myData)

    def generate_fault_yaml_fl33(self, prefix, spatial_order, spatial_zoom, proj):
        cm2m = 0.01
        km2m = 1e3
        self.compute_xy_from_latlon(proj)
        nx, ny = self.nx, self.ny
        p0 = np.array([self.x[0, 0], self.y[0, 0], -km2m * self.depth[0, 0]])
        p1 = np.array([self.x[ny - 1, 0], self.y[ny - 1, 0], -km2m * self.depth[ny - 1, 0]])
        p2 = np.array([self.x[0, nx - 1], self.y[0, nx - 1], -km2m * self.depth[0, nx - 1]])
        hw = p1 - p0
        dx1 = np.linalg.norm(hw) / (ny - 1)
        hw = hw / np.linalg.norm(hw)
        hh = p2 - p0
        dx2 = np.linalg.norm(hh) / (nx - 1)
        hh = hh / np.linalg.norm(hh)
        dx = np.sqrt(self.PSarea_cm2 * cm2m * cm2m)
        # a kinematic model defines the fault quantities at the subfault center
        # a netcdf file defines the quantities at the nodes
        # therefore the dx/2
        # the term dxi/np.sqrt(dx1*dx2) allow accounting for non-square patches
        t1 = -np.dot(p0, hh) + dx * 0.5 * dx1 / np.sqrt(dx1 * dx2)
        t2 = -np.dot(p0, hw) + dx * 0.5 * dx2 / np.sqrt(dx1 * dx2)

        template_yaml = f"""!Switch
[strike_slip, dip_slip, rupture_onset, tau_S, tau_R]: !EvalModel
    parameters: [strike_slip, dip_slip, rupture_onset, effective_rise_time, acc_time]
    model: !Switch
        [strike_slip, dip_slip, rupture_onset, effective_rise_time, acc_time]: !AffineMap
              matrix:
                ua: [{hh[0]}, {hh[1]}, {hh[2]}]
                ub: [{hw[0]}, {hw[1]}, {hw[2]}]
              translation:
                ua: {t1}
                ub: {t2}
              components: !Any
                - !ASAGI
                    file: {prefix}_{spatial_zoom}_o{spatial_order}_ASAGI.nc
                    parameters: [strike_slip, dip_slip, rupture_onset, effective_rise_time, acc_time]
                    var: data
                    interpolation: linear
                - !ConstantMap
                  map:
                    strike_slip: 0.0
                    dip_slip:    0.0
                    rupture_onset:    0.0
                    acc_time:  1e100
                    effective_rise_time:  1e100
    components: !FunctionMap
       map:
          #Note the minus on strike_slip to acknowledge the different convention of SeisSol (T_s>0 means right-lateral)
          strike_slip: return -strike_slip;
          dip_slip: return dip_slip;
          rupture_onset: return rupture_onset;
          tau_S: return acc_time/1.27;
          tau_R: return effective_rise_time - 2.*acc_time/1.27;
        """
        fname = f"{prefix}_fault.yaml"
        with open(fname, "w") as fid:
            fid.write(template_yaml)
        print(f"done writing {fname}")
