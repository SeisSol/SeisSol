import numpy as np
import pyproj
import scipy.ndimage
from scipy import interpolate
from netCDF4 import Dataset
from Yoffe import regularizedYoffe
from scipy import ndimage
from GaussianSTF import GaussianSTF
from scipy.interpolate import RegularGridInterpolator
import xarray as xr


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
        # Create dimension and 1d variables
        sdimVarNames = "uvwxyz"
        dims = []
        for i, xi in enumerate(lDimVar):
            nxi = xi.shape[0]
            dimName = sdimVarNames[i]
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
            arr = np.ascontiguousarray(arr)
            newarr = arr.view(dtype=mattype8)
            newarr = newarr.reshape(newarr.shape[:-1])
            mat = rootgrp.createVariable("data", mat_t, dims)
            mat[:] = newarr


def cosine_taper(npts, p=0.1, freqs=None, flimit=None, halfcosine=True, sactaper=False):
    """
    Cosine Taper. (copied from obspy:
    https://docs.obspy.org/master/_modules/obspy/signal/invsim.html#cosine_taper

    :type npts: int
    :param npts: Number of points of cosine taper.
    :type p: float
    :param p: Decimal percentage of cosine taper (ranging from 0 to 1). Default
        is 0.1 (10%) which tapers 5% from the beginning and 5% form the end.
    :rtype: float NumPy :class:`~numpy.ndarray`
    :return: Cosine taper array/vector of length npts.
    :type freqs: NumPy :class:`~numpy.ndarray`
    :param freqs: Frequencies as, for example, returned by fftfreq
    :type flimit: list or tuple of floats
    :param flimit: The list or tuple defines the four corner frequencies
        (f1, f2, f3, f4) of the cosine taper which is one between f2 and f3 and
        tapers to zero for f1 < f < f2 and f3 < f < f4.
    :type halfcosine: bool
    :param halfcosine: If True the taper is a half cosine function. If False it
        is a quarter cosine function.
    :type sactaper: bool
    :param sactaper: If set to True the cosine taper already tapers at the
        corner frequency (SAC behavior). By default, the taper has a value
        of 1.0 at the corner frequencies.

    .. rubric:: Example

    >>> tap = cosine_taper(100, 1.0)
    >>> tap2 = 0.5 * (1 + np.cos(np.linspace(np.pi, 2 * np.pi, 50)))
    >>> np.allclose(tap[0:50], tap2)
    True
    >>> npts = 100
    >>> p = 0.1
    >>> tap3 = cosine_taper(npts, p)
    >>> (tap3[int(npts*p/2):int(npts*(1-p/2))]==np.ones(int(npts*(1-p)))).all()
    True
    """
    if p < 0 or p > 1:
        msg = "Decimal taper percentage must be between 0 and 1."
        raise ValueError(msg)
    if p == 0.0 or p == 1.0:
        frac = int(npts * p / 2.0)
    else:
        frac = int(npts * p / 2.0 + 0.5)

    if freqs is not None and flimit is not None:
        fl1, fl2, fl3, fl4 = flimit
        idx1 = np.argmin(abs(freqs - fl1))
        idx2 = np.argmin(abs(freqs - fl2))
        idx3 = np.argmin(abs(freqs - fl3))
        idx4 = np.argmin(abs(freqs - fl4))
    else:
        idx1 = 0
        idx2 = frac - 1
        idx3 = npts - frac
        idx4 = npts - 1
    if sactaper:
        # in SAC the second and third
        # index are already tapered
        idx2 += 1
        idx3 -= 1

    # Very small data lengths or small decimal taper percentages can result in
    # idx1 == idx2 and idx3 == idx4. This breaks the following calculations.
    if idx1 == idx2:
        idx2 += 1
    if idx3 == idx4:
        idx3 -= 1

    # the taper at idx1 and idx4 equals zero and
    # at idx2 and idx3 equals one
    cos_win = np.zeros(npts)
    if halfcosine:
        # cos_win[idx1:idx2+1] =  0.5 * (1.0 + np.cos((np.pi * \
        #    (idx2 - np.arange(idx1, idx2+1)) / (idx2 - idx1))))
        cos_win[idx1 : idx2 + 1] = 0.5 * (1.0 - np.cos((np.pi * (np.arange(idx1, idx2 + 1) - float(idx1)) / (idx2 - idx1))))
        cos_win[idx2 + 1 : idx3] = 1.0
        cos_win[idx3 : idx4 + 1] = 0.5 * (1.0 + np.cos((np.pi * (float(idx3) - np.arange(idx3, idx4 + 1)) / (idx4 - idx3))))
    else:
        cos_win[idx1 : idx2 + 1] = np.cos(-(np.pi / 2.0 * (float(idx2) - np.arange(idx1, idx2 + 1)) / (idx2 - idx1)))
        cos_win[idx2 + 1 : idx3] = 1.0
        cos_win[idx3 : idx4 + 1] = np.cos((np.pi / 2.0 * (float(idx3) - np.arange(idx3, idx4 + 1)) / (idx4 - idx3)))

    # if indices are identical division by zero
    # causes NaN values in cos_win
    if idx1 == idx2:
        cos_win[idx1] = 0.0
    if idx3 == idx4:
        cos_win[idx3] = 0.0
    return cos_win


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


def compute_block_mean(ar, fact):
    a = xr.DataArray(ar, dims=["x", "y"])
    return a.coarsen(x=fact, y=fact).mean().to_numpy()


def upsample_quantities(allarr, spatial_order, spatial_zoom, padding="constant", extra_padding_layer=False, minimize_block_average_variations=False):
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
            my_array0 = np.pad(allarr[k, :, :], ((1, 1), (1, 1)), "reflect", reflect_type="odd")
        else:
            my_array0 = np.pad(allarr[k, :, :], ((1, 1), (1, 1)), padding)
        if extra_padding_layer:
            ncrop = spatial_zoom - 1
        else:
            ncrop = spatial_zoom
        my_array = scipy.ndimage.zoom(my_array0, spatial_zoom, order=spatial_order, mode="grid-constant", grid_mode=True)
        if minimize_block_average_variations:
            # inspired by Tinti et al. (2005) (Appendix A)
            # This is for the specific case of fault slip.
            # We want to preserve the seismic moment of each subfault after interpolation
            # the rock rigidity is not know by this script (would require some python binding of easi).
            # the subfault area is typically constant over the kinematic model
            # So we just want to perserve subfault average.
            print("trying to perserve subfault average...")
            my_array = np.maximum(0, my_array)
            best_misfit = float("inf")
            # The algorithm does not seem to converge, but produces better model
            # (given the misfit) that inital after 2-3 iterations
            niter = 30
            for i in range(niter):
                block_average = compute_block_mean(my_array, spatial_zoom)
                correction = my_array0 / block_average
                # having a misfit as misfit = np.linalg.norm(correction) does not makes sense as for almost 0 slip, correction can be large
                misfit = np.linalg.norm(my_array0 - block_average) / len(my_array0)
                if best_misfit > misfit:
                    if i == 0:
                        print(f"misfit at iter {i}: {misfit}")
                    else:
                        print(f"misfit improved at iter {i}: {misfit}")
                    best_misfit = misfit
                    best = np.copy(my_array)
                my_array = scipy.ndimage.zoom(correction * my_array0, spatial_zoom, order=spatial_order, mode="grid-constant", grid_mode=True)
                my_array = np.maximum(0, my_array)
            my_array = best
        if ncrop > 0:
            allarr0[k, :, :] = my_array[ncrop:-ncrop, ncrop:-ncrop]
        else:
            allarr0[k, :, :] = my_array

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

    def assess_STF_parameters(self, threshold):
        "compute rise_time (slip duration) and t_acc (peak SR) from SR time histories"
        assert threshold >= 0.0 and threshold < 1
        self.rise_time = np.zeros((self.ny, self.nx))
        self.tacc = np.zeros((self.ny, self.nx))
        misfits_Yoffe = []
        misfits_Gaussian = []
        for j in range(self.ny):
            for i in range(self.nx):
                if not self.slip1[j, i]:
                    self.rise_time[j, i] = np.nan
                    self.tacc[j, i] = np.nan
                else:
                    peakSR = np.amax(self.aSR[j, i, :])
                    id_max = np.where(self.aSR[j, i, :] == peakSR)[0][0]
                    ids_greater_than_threshold = np.where(self.aSR[j, i, :] > threshold * peakSR)[0]
                    first_non_zero = np.amin(ids_greater_than_threshold)
                    last_non_zero = np.amax(ids_greater_than_threshold)
                    self.rise_time[j, i] = (last_non_zero - first_non_zero + 1) * self.dt
                    self.tacc[j, i] = (id_max - first_non_zero + 1) * self.dt
                    t0_increment = first_non_zero * self.dt
                    self.t0[j, i] += t0_increment
                    # 2 dims: 0: Yoffe 1: Gaussian
                    newSR = np.zeros((self.ndt, 2))
                    # Ts and Td parameters of the Yoffe function have no direct physical meaning
                    # Tinti et al. (2005) suggest that Ts can nevertheless be associated with the acceleration time tacc
                    # Empirically, they find that Ts and Tacc are for most (Ts,Td) parameter sets linearly related
                    # with the 'magic' number 1.27
                    ts = self.tacc[j, i] / 1.27
                    tr = self.rise_time[j, i] - 2.0 * ts
                    tr = max(tr, ts)
                    for k, tk in enumerate(self.myt):
                        newSR[k, 0] = regularizedYoffe(tk - t0_increment, ts, tr)
                        newSR[k, 1] = GaussianSTF(tk - t0_increment, self.rise_time[j, i], self.dt)
                    integral_aSTF = np.trapz(np.abs(self.aSR[j, i, :]), dx=self.dt)
                    integral_Yoffe = np.trapz(np.abs(newSR[:, 0]), dx=self.dt)
                    integral_Gaussian = np.trapz(np.abs(newSR[:, 1]), dx=self.dt)
                    if integral_aSTF > 0:
                        misfits_Yoffe.append(np.linalg.norm(self.aSR[j, i, :] / integral_aSTF - newSR[:, 0] / integral_Yoffe))
                        misfits_Gaussian.append(np.linalg.norm(self.aSR[j, i, :] / integral_aSTF - newSR[:, 1] / integral_Gaussian))
        misfits_Yoffe = np.array(misfits_Yoffe)
        misfits_Gaussian = np.array(misfits_Gaussian)
        print(f"misfit Yoffe (10-50-90%): {np.percentile(misfits_Yoffe,10):.2f} {np.percentile(misfits_Yoffe,50):.2f} {np.percentile(misfits_Yoffe,90):.2f}")
        print(f"misfit Gaussian (10-50-90%): {np.percentile(misfits_Gaussian,10):.2f} {np.percentile(misfits_Gaussian,50):.2f} {np.percentile(misfits_Gaussian,90):.2f}")

        self.rise_time = interpolate_nan_from_neighbors(self.rise_time)
        self.tacc = interpolate_nan_from_neighbors(self.tacc)

        print("slip rise_time (min, 50%, max)", np.amin(self.rise_time), np.median(self.rise_time), np.amax(self.rise_time))
        print("tacc (min, 50%, max)", np.amin(self.tacc), np.median(self.tacc), np.amax(self.tacc))

    def upsample_fault(self, spatial_order, spatial_zoom, temporal_zoom, proj, use_Yoffe=False, time_smoothing_kernel_as_dt_fraction=0.5):
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
        (pf.slip1,) = upsample_quantities(allarr, spatial_order, spatial_zoom, padding="constant", minimize_block_average_variations=True)
        pf.compute_latlon_from_xy(proj)
        pf.PSarea_cm2 = self.PSarea_cm2 / spatial_zoom**2
        ratio_potency = np.sum(pf.slip1) * pf.PSarea_cm2 / (np.sum(self.slip1) * self.PSarea_cm2)
        print(f"seismic potency ratio (upscaled over initial): {ratio_potency}")

        if use_Yoffe:
            allarr = np.array([self.rise_time, self.tacc])
            pf.rise_time, pf.tacc = upsample_quantities(allarr, spatial_order, spatial_zoom, padding="edge")
            pf.rise_time = np.maximum(pf.rise_time, np.amin(self.rise_time))
            pf.tacc = np.maximum(pf.tacc, np.amin(self.tacc))
            # see comment above explaining why the 1.27 factor
            print("using ts = tacc / 1.27 to compute the regularized Yoffe")
            ts = pf.tacc / 1.27
            tr = pf.rise_time - 2.0 * ts
            tr = np.maximum(tr, ts)
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
                    # 1. upsample with linear interpolation
                    # 2. apply a gauss kernel to smooth out sharp edges
                    # 3. tapper the signal smoothly to 0 at both time ends
                    # 4. rescale SR to ensure integral (SR) = slip
                    f = interpolate.interp1d(self.myt, aSRa[j, i, :], kind="linear")
                    pf.aSR[j, i, :] = f(pf.myt)
                    tapper = cosine_taper(pf.ndt, self.dt / (pf.ndt * pf.dt))
                    pf.aSR[j, i, :] = tapper * ndimage.gaussian_filter1d(pf.aSR[j, i, :], time_smoothing_kernel_as_dt_fraction * self.dt / pf.dt, mode="constant")
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

    def compute_1d_dimension_arrays(self, spatial_zoom):
        self.spatial_zoom = spatial_zoom
        # Compute dimension arrays
        km2m = 1e3
        coords = np.array([self.x, self.y, -km2m * self.depth])
        ny, nx = coords.shape[1:3]
        center_row = coords[:, (ny - 1) // 2, :] - coords[:, (ny - 1) // 2 - 1, :]
        dx1 = np.linalg.norm(center_row, axis=0)
        # with this convention the first data point is in local coordinate (0,0)
        xb = np.cumsum(dx1) - dx1[0]

        center_col = coords[:, :, (nx - 1) // 2] - coords[:, :, (nx - 1) // 2 - 1]
        dy1 = np.linalg.norm(center_col, axis=0)
        yb = np.cumsum(dy1) - dy1[0]

        self.xb = np.pad(xb, ((1), (1)), "reflect", reflect_type="odd")
        self.yb = np.pad(yb, ((1), (1)), "reflect", reflect_type="odd")

        # we want to cover all the fault, that is up to -dx/2.
        # With this strategy We will cover a bit more than that, but it is probably not a big deal

        ncrop = spatial_zoom - 1
        self.x_up = scipy.ndimage.zoom(self.xb, spatial_zoom, order=1, mode="grid-constant", grid_mode=True)[ncrop:-ncrop]
        self.y_up = scipy.ndimage.zoom(self.yb, spatial_zoom, order=1, mode="grid-constant", grid_mode=True)[ncrop:-ncrop]

        # used for the interpolation
        yg, xg = np.meshgrid(self.y_up, self.x_up)
        self.yx = np.array([yg.ravel(), xg.ravel()]).T

    def upsample_quantity_RGInterpolator_core(self, arr, method, is_slip=False):
        if is_slip:
            # tapper to 0 slip except at the top
            print("tapper slip to 0, except at the top (hardcoded)")
            padded_arr = np.pad(arr, ((1, 0), (1, 1)), "constant")
            padded_arr = np.pad(padded_arr, ((0, 1), (0, 0)), "edge")
        else:
            padded_arr = np.pad(arr, ((1, 1), (1, 1)), "edge")
        interp = RegularGridInterpolator([self.yb, self.xb], padded_arr)
        return interp(self.yx, method=method).reshape(self.x_up.shape[0], self.y_up.shape[0]).T

    def upsample_quantity_RGInterpolator(self, arr, method, is_slip=False):
        my_array = self.upsample_quantity_RGInterpolator_core(arr, method, is_slip)
        minimize_block_average_variations = is_slip
        if minimize_block_average_variations:
            # inspired by Tinti et al. (2005) (Appendix A)
            # This is for the specific case of fault slip.
            # We want to preserve the seismic moment of each subfault after interpolation
            # the rock rigidity is not know by this script (would require some python binding of easi).
            # the subfault area is typically constant over the kinematic model
            # So we just want to perserve subfault average.
            print("trying to perserve subfault average...")
            my_array = np.maximum(0, my_array)
            best_misfit = float("inf")
            # The algorithm does not seem to converge, but produces better model
            # (given the misfit) that inital after 2-3 iterations
            niter = 3
            for i in range(niter):
                block_average = compute_block_mean(my_array[1:-1, 1:-1], self.spatial_zoom)
                print(arr.shape, block_average.shape)
                correction = arr / block_average
                # having a misfit as misfit = np.linalg.norm(correction) does not makes sense as for almost 0 slip, correction can be large
                misfit = np.linalg.norm(arr - block_average) / len(arr)
                if best_misfit > misfit:
                    if i == 0:
                        print(f"misfit at iter {i}: {misfit}")
                    else:
                        print(f"misfit improved at iter {i}: {misfit}")
                    best_misfit = misfit
                    best = np.copy(my_array)
                my_array = self.upsample_quantity_RGInterpolator_core(correction * arr, method, is_slip)
                my_array = np.maximum(0, my_array)
            my_array = best
        return my_array

    def generate_netcdf_fl33(self, prefix, method, spatial_zoom, proj, write_paraview):
        "generate netcdf files to be used with SeisSol friction law 33"

        cm2m = 0.01
        km2m = 1e3
        # a kinematic model defines the fault quantities at the subfault center
        # a netcdf file defines the quantities at the nodes
        # therefore the extra_padding_layer=True, and the added di below
        cslip = self.compute_corrected_slip_for_differing_area(proj)
        self.compute_1d_dimension_arrays(spatial_zoom)

        upsampled_arrays = []

        slip = self.upsample_quantity_RGInterpolator(cslip, method, is_slip=True)
        for arr in [cslip, self.t0, self.rake, self.rise_time, self.tacc]:
            upsampled_arrays.append(self.upsample_quantity_RGInterpolator(arr, method))

        slip, rupttime, rake, rise_time, tacc = upsampled_arrays

        # upsampled duration, rise_time and acc_time may not be smaller than initial values
        # at least rise_time could lead to a non-causal kinematic model
        rupttime = np.maximum(rupttime, np.amin(self.t0))
        rise_time = np.maximum(rise_time, np.amin(self.rise_time))
        tacc = np.maximum(tacc, np.amin(self.tacc))

        rake_rad = np.radians(rake)
        strike_slip = slip * np.cos(rake_rad) * cm2m
        dip_slip = slip * np.sin(rake_rad) * cm2m

        dx = np.sqrt(self.PSarea_cm2 * cm2m * cm2m)
        ldataName = ["strike_slip", "dip_slip", "rupture_onset", "effective_rise_time", "acc_time"]
        lgridded_myData = [strike_slip, dip_slip, rupttime, rise_time, tacc]

        prefix2 = f"{prefix}_{spatial_zoom}_{method}"
        if write_paraview:
            # see comment above
            for i, sdata in enumerate(ldataName):
                writeNetcdf(
                    f"{prefix2}_{sdata}",
                    [self.x_up, self.y_up],
                    [sdata],
                    [lgridded_myData[i]],
                    paraview_readable=True,
                )
            writeNetcdf(f"{prefix2}_slip", [self.x_up, self.y_up], ["slip"], [slip], paraview_readable=True)
        writeNetcdf(prefix2, [self.x_up, self.y_up], ldataName, lgridded_myData)

    def generate_fault_ts_yaml_fl33(self, prefix, method, spatial_zoom, proj):
        """Generate yaml file initializing FL33 arrays and ts file describing the planar fault geometry."""
        # Generate yaml file loading ASAGI file
        cm2m = 0.01
        km2m = 1e3
        self.compute_xy_from_latlon(proj)
        nx, ny = self.nx, self.ny
        p0 = np.array([self.x[0, 0], self.y[0, 0], -km2m * self.depth[0, 0]])
        p1 = np.array([self.x[ny - 1, 0], self.y[ny - 1, 0], -km2m * self.depth[ny - 1, 0]])
        p2 = np.array([self.x[0, nx - 1], self.y[0, nx - 1], -km2m * self.depth[0, nx - 1]])
        p3 = np.array([self.x[ny - 1, nx - 1], self.y[ny - 1, nx - 1], -km2m * self.depth[ny - 1, nx - 1]])

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
        # the term dxi/np.sqrt(dx1*dx2) allows accounting for non-square patches
        non_square_factor = dx / np.sqrt(dx1 * dx2)
        t1 = -np.dot(p0, hh)
        t2 = -np.dot(p0, hw)

        template_yaml = f"""!Switch
[strike_slip, dip_slip, rupture_onset, tau_S, tau_R, rupture_rise_time]: !EvalModel
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
                    file: {prefix}_{spatial_zoom}_{method}.nc
                    parameters: [strike_slip, dip_slip, rupture_onset, effective_rise_time, acc_time]
                    var: data
                    interpolation: linear
                - !ConstantMap
                  map:
                    strike_slip: 0.0
                    dip_slip:    0.0
                    rupture_onset:    0.0
                    acc_time:  1e100
                    effective_rise_time:  2e100
    components: !FunctionMap
       map:
          #Note the minus on strike_slip to acknowledge the different convention of SeisSol (T_s>0 means right-lateral)
          strike_slip: return -strike_slip;
          dip_slip: return dip_slip;
          rupture_onset: return rupture_onset;
          tau_S: return acc_time/1.27;
          tau_R: return effective_rise_time - 2.*acc_time/1.27;
          rupture_rise_time: return effective_rise_time;
        """
        fname = f"{prefix}_fault.yaml"
        with open(fname, "w") as fid:
            fid.write(template_yaml)
        print(f"done writing {fname}")

        # Generate ts file containing mesh geometry
        vertex = np.zeros((4, 3))
        vertex[0, :] = p0 + 0.5 * (-hh * dx1 - hw * dx2) * non_square_factor
        vertex[1, :] = p2 + 0.5 * (hh * dx1 - hw * dx2) * non_square_factor
        vertex[2, :] = p3 + 0.5 * (hh * dx1 + hw * dx2) * non_square_factor
        vertex[3, :] = p1 + 0.5 * (-hh * dx1 + hw * dx2) * non_square_factor

        connect = np.zeros((2, 3), dtype=int)
        connect[0, :] = [1, 2, 3]
        connect[1, :] = [1, 3, 4]
        fname = f"{prefix}_fault.ts"
        with open(fname, "w") as fout:
            fout.write("GOCAD TSURF 1\nHEADER {\nname:%s\nborder: true\nmesh: false\n*border*bstone: true\n}\nTFACE\n" % (fname))
            for ivx in range(1, 5):
                fout.write("VRTX %s %s %s %s\n" % (ivx, vertex[ivx - 1, 0], vertex[ivx - 1, 1], vertex[ivx - 1, 2]))

            for i in range(2):
                fout.write("TRGL %d %d %d\n" % (connect[i, 0], connect[i, 1], connect[i, 2]))
            fout.write("END\n")
        print(f"done writing {fname}")
