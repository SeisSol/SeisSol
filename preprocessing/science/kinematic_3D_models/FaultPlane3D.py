import numpy as np
import pyproj
import h5py
import copy
import scipy
from netCDF4 import Dataset
from scipy import interpolate, ndimage
from scipy.interpolate import RegularGridInterpolator, interp2d
from scipy.spatial import cKDTree as KDTree
import os

""" 
Generalize FaultPlane.py by Thomas to non-planar SRF and varying rake angle

- Functions:
    WriteNetcdf
    upsample_quantities       - allow negative slip in seismic potency conservation
    compute_block_mean
    generate_netcdf_fl33
    compute_grid              - for ASAGI strutural mesh
    compute_distance_to_cloud - for ASAGI strutural mesh
    interpolate_3D            - for ASAGI strutural mesh
    
    
- Object - Fault:
    - attributes
        - Basic representation of SRF model
    - Functions:
        - File i/o   - updated to point-wise SRF file
        - Project to and from local coordinate 
        - Upscale
        - Area correction - updated with considering varying dip angle
        - Write FL33 required files

Usage: refer to Tohoku_gen_FL33_notebook jupter notebook \
        or Gerneate_non_planar_FL33.py
"""



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
            newarr = arr.view(dtype=mattype8)
            newarr = newarr.reshape(newarr.shape[:-1])
            mat = rootgrp.createVariable("data", mat_t, dims)
            mat[:] = newarr



def compute_block_mean(ar, fact):
    """
    dowsample array ar by factor fact
    https://stackoverflow.com/questions/18666014/downsample-array-in-python
    """
    assert isinstance(fact, int), type(fact)
    sx, sy = ar.shape
    X, Y = np.ogrid[0:sx, 0:sy]
    regions = sy // fact * (X // fact) + Y // fact
    res = ndimage.mean(ar, labels=regions, index=np.arange(regions.max() + 1)) 
    res.shape = (sx // fact, sy // fact)
    return res 



def upsample_quantities(allarr, spatial_order, spatial_zoom, padding="constant", extra_padding_layer=False, \
                        minimize_block_average_variations=False, \
                        minimize_block_average_variations_non_pos=False, verbose=False):
    """
    1. pad
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
        
        # Main interpolation process
        my_array = ndimage.zoom(my_array0, spatial_zoom, order=spatial_order, \
                                mode="grid-constant", grid_mode=True)
        
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
                elif verbose:
                    print(f"Misfit did not improved at iter{i}: {misfit}")
                my_array = ndimage.zoom(correction * my_array0, spatial_zoom, \
                                        order=spatial_order, mode="grid-constant", grid_mode=True)
                my_array = np.maximum(0, my_array)
            my_array = best
        if minimize_block_average_variations_non_pos:
            # this special case is allowing negative slip by computing the correction in
            # absolute slip
            print("trying to perserve subfault average...")
            #my_array = np.maximum(0, my_array)
            best_misfit = float("inf")
            
            # The algorithm does not seem to converge, but produces better model
            # (given the misfit) that inital after 2-3 iterations
            niter = 30
            for i in range(niter):
                block_average = compute_block_mean(np.abs(my_array), spatial_zoom)
                correction = np.abs(my_array0) / block_average
                # having a misfit as misfit = np.linalg.norm(correction) does not makes sense as for almost 0 slip, correction can be large
                misfit = np.linalg.norm(np.abs(my_array0) - block_average) / len(my_array0)
                if best_misfit > misfit:
                    if i == 0:
                        print(f"misfit at iter {i}: {misfit}")
                    else:
                        print(f"misfit improved at iter {i}: {misfit}")
                    best_misfit = misfit
                    best = np.copy(my_array)
                elif verbose:
                    print(f"Misfit did not improved at iter{i}: {misfit}")
                my_array = ndimage.zoom(correction * my_array0, spatial_zoom, \
                                        order=spatial_order, mode="grid-constant", grid_mode=True)
                #my_array = np.maximum(0, my_array)
            my_array = best
        if ncrop > 0:
            allarr0[k, :, :] = my_array[ncrop:-ncrop, ncrop:-ncrop]

    return allarr0


### used for computing Slip rate function (Tinti, 2005) (not used in Tohoku case)
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
    return interpolate.griddata((x1, y1), newarr.ravel(), (xx, yy), \
                                method="linear", fill_value=np.average(array))



### section for upscaled fault to FL33 NetCDF
def compute_grid(xyzc, nx, ny):
    x0, x1 = np.amin(xyzc[:, 0]), np.amax(xyzc[:, 0])
    y0, y1 = np.amin(xyzc[:, 1]), np.amax(xyzc[:, 1])
    X = np.linspace(x0, x1, nx)
    Y = np.linspace(y0, y1, ny)
    Xg, Yg = np.meshgrid(X, Y)
    return X,Y,Xg,Yg


def compute_distance_to_cloud(Xg, Yg, xyzc):
    tree = KDTree(xyzc[:, 0:2])
    dist, _ = tree.query(np.c_[Xg.ravel(), Yg.ravel()], k=1)
    dist = dist.reshape(Xg.shape)
    return dist


def interpolate_3D(xyzc, Xg, Yg, dist_to_cloud, inputVar):
    gOutput = scipy.interpolate.griddata(
        xyzc[:, 0:2], inputVar[:], (Xg, Yg), method="linear"
    )

    # remove points too far away from the point Cloud
    gOutput[dist_to_cloud > 15e3] = 0.
    gOutput[np.isnan(gOutput)]  = 0.
    return gOutput.reshape(Xg.shape)


class FaultPlane:
    def __init__(self):
        self.nx = 0
        self.ny = 0
        self.nz = 0
        self.ndt = 0
        self.dt = 0
        
        self.lat0 = 0
        self.lon0 = 0
        self.lat = 0
        self.lon = 0
        self.depth = 0
        self.t0 = 0
        self.slip1 = 0
        self.slip2 = 0
        self.strike = 0
        self.dip = 0
        self.rake = 0
        
        # two direction of slip is used here
        self.aSR1 = 0
        self.aSR2 = 0
        self.myt = 0

    def init_spatial_arrays(self, nx, ny, nz=None):
        self.nx = nx
        self.ny = ny
        self.lon = np.zeros((ny,nx))
        self.lat = np.zeros((ny,nx))
        self.x = np.zeros((ny,nx))
        self.y = np.zeros((ny,nx))
        self.depth = np.zeros((ny,nx))
        self.t0 = np.zeros((ny,nx))
        self.slip1 = np.zeros((ny,nx))
        self.slip2 = np.zeros((ny,nx))
        self.strike = np.zeros((ny,nx))
        self.dip = np.zeros((ny,nx))
        self.rake = np.zeros((ny,nx))

    def init_aSR(self):
        # aSR1 along strike, aSR2 along dip direction
        self.aSR1 = np.zeros((self.ny,self.nx,self.ndt))
        self.aSR2 = np.zeros((self.ny,self.nx,self.ndt))

    def extend_aSR(self, ndt_old, ntd_new):
        "extend aSR array to more time stepping"
        self.ndt = ndt_new
        tempSR1 = self.aSR1
        tempSR2 = self.aSR2
        self.aSR1 = np.zeros((self.ny, self.nx, self.ndt))
        self.aSR2 = np.zeros((self.ny, self.nx, self.ndt))
        self.aSR1[:, :, 0:ndt_old] = tmpSR1[:, :, :]
        self.aSR2[:, :, 0:ndt_old] = tmpSR2[:, :, :]

    def compute_xy_from_latlon(self, proj):
        if proj:
            from pyproj import Transformer
            transformer = Transformer.from_crs('epsg:4326', \
                                                proj[0], \
                                                always_xy=True)
            self.x, self.y = transformer.transform(self.lon, self.lat)
        else:
            print('No proj string specified')
            self.x, self.y = self.lon, self.lat

    def compute_latlon_from_xy(self, proj):
        if proj:
            from pyproj import Transformer
            transformer = Transformer.from_crs(proj[0], \
                                                'epsg:4326',\
                                                always_xy=True)
            self.lon, self.lat = transformer.transform(self.x, self.y)
        else:
            print('No proj string specified')
            self.lon, self.lat = self.x, self.y

    def compute_time_array(self):
        self.myt = np.linspace(0, (self.ndt-1)*self.dt, self.ndt)


    def write_srf(self, fname):
        """
        Write kinematic model to a srf file (standard rupture format)
        and include .init_spatial_array
        """
        with open(fname,"w") as fout:
            fout.write("1.0\n")
            fout.write("POINTS %d\n" % (self.nx * self.ny))
            for j in range(self.ny):
                for i in range(self.nx):
                    fout.write("%g %g %g %g %g %e %g %g\n" % \
                            (self.lon[j, i], self.lat[j, i], self.depth[j, i], 
                             self.strike[j, i], self.dip[j, i], self.PSarea_cm2,\
                             self.t0[j, i], self.dt))
                    fout.write("%g %g %d %f %d %f %d\n" % (self.rake[j, i], \
                                                           self.slip1[j, i], \
                                                           self.ndt, \
                                                           self.slip2[j,i], \
                                                           self.ndt, 0.0, 0))
                    np.savetxt(fout, self.aSR1[j, i, :], fmt="%g", newline=" ")
                    fout.write("\n")
                    np.savetxt(fout, self.aSR2[j, i, :], fmt="%g", newline=" ")
                    fout.write("\n")
        print("done writing", fname)

    def init_from_srf(self, fname, nx, ny):
        """Init object by reading a srf file (standard rupture format)
        nx ny requried for point-wise SRF file
        """
        self.nx, self.ny = int(nx), int(ny)
        with open(fname, "r") as fid:
            # version
            line = fid.readline()
            version = float(line)
            if not (abs(version - 1.0) < 1e-3 or abs(version - 2.0) < 1e-3):
                print("srf version: %s not supported" % (line))
                raise
            # skip comments
            while True:
                line = fid.readline()
                if not line.startswith("#"):
                    break
            # PLANE
            line_el = line.split()
            if line_el[0]=="PLANE":
               print("Plane is skiped")
               skip_line = int(line_el[1])
               for i in range(skip_line*2):
                   _ = fid.readline()
            # POINTS
            line_el = line.split()
            assert line_el[0]=="POINTS", "no points specified"
            assert int(line_el[1])==nx*ny, "no. of points does not match nx*ny"
            # init_spatial_arrays
            self.init_spatial_arrays(nx, ny)
            
            for j in range(ny):
                for i in range(nx):
                    # first header line
                    line = fid.readline()
                    self.lon[j, i], self.lat[j, i], self.depth[j, i], self.strike[j, i], self.dip[j, i], self.PSarea_cm2, self.t0[j, i], dt, *rho_vs = [float(v) for v in line.split()]
                    # second header line
                    line = fid.readline()
                    self.rake[j, i], self.slip1[j, i], ndt1, self.slip2[j, i], ndt2, slip3, ndt3 = [float(v) for v in line.split()]
                    assert slip3 == 0.0, 'this scripts assumes no fault opening (slip3 is zero)'

                    ndt1 = int(ndt1)
                    if max(i,j)==0:
                        self.ndt = ndt1
                        self.dt = dt
                        self.init_aSR()
                    STF1, STF2 = [], []
                    if ndt1==0:
                        continue
                    if ndt1 > self.ndt:
                        print(f"a larger ndt ({ndt1}> {self.ndt}) was found for point source (i,j) = ({i}, {j}) extending aSR array...")
                        self.extend_aSR(self.ndt, ndt1)
                    assert abs(dt-self.dt) < 1e-6, "this script assumes dt are the same for all sources"
                    while True:
                        line = fid.readline()
                        STF1.extend(line.split())
                        line = fid.readline()
                        STF2.extend(line.split())
                        if len(STF1) == ndt1:
                            self.aSR1[j, i, 0:ndt1] = np.array([float(v) for v in STF1])
                            self.aSR2[j, i, 0:ndt1] = np.array([float(v) for v in STF2])
                            break

    def assess_STF_parameters(self):
        """
        Depreciated for 3D fault model
        compute rise_time (slip duration) and t_acc (peak SR) from SR time histories
        """
        self.rise_time = np.zeros((self.ny, self.nx))
        self.tacc = np.zeros((self.ny, self.nx))
        for j in range(self.ny):
            for i in range(self.nx):
                if not self.slip1[j,i]:
                    self.rise_time[j,i] = np.nan
                    self.taac[j,i] = np.nan
                else:
                    first_non_zero = np.amin(np.where(self.aSR1[j,i,:])[0])                   
                    last_non_zero = np.amax(np.where(self.aSR1[j, i, :])[0])
                    id_max = np.where(self.aSR[j, i, :] == np.amax(self.aSR[j, i, :]))[0]
                    self.rise_time[j, i] = (last_non_zero - first_non_zero + 1) * self.dt
                    self.tacc[j, i] = (id_max - first_non_zero + 1) * self.dt
                    self.t0[j, i] += first_non_zero * self.dt
        self.rise_time = interpolate_nan_from_neighbors(self.rise_time)
        self.tacc = interpolate_nan_from_neighbors(self.tacc)

        print("slip rise_time (min, 50%, max)", np.amin(self.rise_time), np.median(self.rise_time), np.amax(self.rise_time))
        print("tacc (min, 50%, max)", np.amin(self.tacc), np.median(self.tacc), np.amax(self.tacc))


    def upsample_fault_slip(self, spatial_order, spatial_zoom, proj, verbose=False):
        """ a simplified version of upsample_fault 
        - Upsample fault slip (slip1 , slip2) without the slip rate funcion
        - Upsample dip, strike to compute the slip vector
        - To preserve rupture reaching boundary, we updated paddint to nearest
        
        output: upsampled fault object
        
        Note: minimize_block_average_variations from Thomas would limit positive slip
        """
        # time vector
        ny2, nx2 = self.ny * spatial_zoom, self.nx * spatial_zoom
        
        # resample source
        pf = FaultPlane()
        pf.init_spatial_arrays(nx2, ny2)
        pf.ndt = self.ndt
        pf.init_aSR()
        
        
        # upsample spatially geometry (bilinear interpolation)
        allarr = np.array([self.x, self.y, self.depth])
        pf.x, pf.y, pf.depth = upsample_quantities(allarr, spatial_order=1, spatial_zoom=spatial_zoom, padding="extrapolate")
        
        # upsample other quantities # edge  (bilinear interpolation)
        allarr = np.array([self.t0, self.strike, self.dip, self.rake])
        pf.t0, pf.strike, pf.dip, pf.rake = upsample_quantities(allarr, spatial_order=1, spatial_zoom=spatial_zoom, padding="edge")
        # the interpolation may generate some acausality that we here prevent
        pf.t0 = np.maximum(pf.t0, np.amin(self.t0))
        
        # To keep rupture reaching the boundary, we use padding='edge'
        allarr = np.array([self.slip1])
        (pf.slip1,) = upsample_quantities(allarr, spatial_order, spatial_zoom, padding="edge", \
                                          minimize_block_average_variations_non_pos=True, verbose=False)
        
        allarr = np.array([self.slip2])
        (pf.slip2,) = upsample_quantities(allarr, spatial_order, spatial_zoom, padding="edge", 
                                          minimize_block_average_variations_non_pos=True, verbose=False)
        
        pf.compute_latlon_from_xy(proj)
        pf.PSarea_cm2 = self.PSarea_cm2 / spatial_zoom**2
        ratio_potency1 = np.sum(np.abs(pf.slip1)) * pf.PSarea_cm2/ (np.sum(np.abs(self.slip1)) * self.PSarea_cm2)
        ratio_potency2 = np.sum(np.abs(pf.slip2)) * pf.PSarea_cm2/ (np.sum(np.abs(self.slip2)) * self.PSarea_cm2)
        print(f"seismic potency ratio - slip direction 1(upscaled over initial): {ratio_potency1}" )
        print(f"seismic potency ratio - slip direction 2(upscaled over initial): {ratio_potency2}" )
        
        return pf
        
        
    def upsample_fault(self, spatial_order, spatial_zoom, temporal_zoom, proj, \
                       use_Yoffe=False, time_smoothing_kernal_as_dt_fraction=0.5, verbose=False):
        """ Increase spatial and temproal resolution model by interpolation"""
        # time vector
        ndt2 = (self.ndt - 1)*temporal_zoom + 1
        ny2, nx2 = self.ny * spatial_zoom, self.nx * spatial_zoom
        
        # resample source
        pf = FaultPlane()
        pf.init_spatial_arrays(nx2, ny2)
        pf.ndt = ndt2
        pf.init_aSR()
        
        pf.dt = self.dt/ temporal_zoom
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
        (pf.slip1,) = upsample_quantities(allarr, spatial_order, spatial_zoom, padding="constant", \
                                          minimize_block_average_variations=True, verbose=False)
        
        allarr = np.array([self.slip2])
        (pf.slip2,) = upsample_quantities(allarr, spatial_order, spatial_zoom, padding="constant", \
                                          minimize_block_average_variations=True, verbose=False)
        
        pf.compute_latlon_from_xy(proj)
        pf.PSarea_cm2 = self.PSarea_cm2 / spatial_zoom**2
        ratio_potency1 = np.sum(pf.slip1) * pf.PSarea_cm2/ (np.sum(self.slip1) * self.PSarea_cm2)
        ratio_potency2 = np.sum(pf.slip2) * pf.PSarea_cm2/ (np.sum(self.slip2) * self.PSarea_cm2)
        print(f"seismic potency ratio - slip direction 1(upscaled over initial): {ratio_potency1}" )
        print(f"seismic potency ratio - slip direction 2(upscaled over initial): {ratio_potency2}" )
        
        if use_Yoffe:
            self.assess_STF_parameters()
            allarr = np.array([self.rise_time, self.tacc])
            pf.rise_time, pf.tacc = upsample_quantities(allarr, spatial_order, spatial_zoom, padding="edge")
            pf.rise_time = np.maximum(pf.ries_time, np.amin(self.rise_time))
            pf.tacc = np.maximum(pf.tacc, np.amin(self.tacc))
            print("using ts = tacc / 1.27 to compute the regularized Yoffe")
            ts = pf.tacc / 1.27
            tr = pf.rise_time - 2.0 * ts
            for j in range(pf.ny):
                for i in range(pf.nx):
                    for k, tk in enumerate(pf.myt):
                        pf.aSR1[j, i, k] = pf.slip1[j, i] * regularizedYoffe(tk, ts[j, i], tr[j, i])
                        pf.aSR2[j, i, k] = pf.slip2[j, i] * regularizedYoffe(tk, ts[j, i], tr[j, i])
        else:
            aSRa = np.zeros((pf.ny, pf.nx, self.ndt))
            for k in range(self.ndt):
                aSRa[:, :, k] = upsample_quantities(np.array([self.aSR[:, :, k]]), spatial_order, spatial_zoom, padding="constant")
            
            # interpolate temporally the AST
            for j in range(pf.ny):
                for i in range(pf.nx):
                    # 1. upsample with linear interpolation
                    # 2. apply a gauss kernal to smooth out sharp edges
                    # 3. tapper the signal smoothly to 0 at both ends
                    # 4. rescale SR ensure integral (SR) = slip
                    f = interpolate.interp1d(self.myt, aSRa[j,i,:], kind="linear")
                    pf.aSR[j, i, :] = f(pf.myt)
                    tapper = cosine_taper(pf.ndt, self. dt/ (pf.ndt * pf.dt))
                    pf.aSR[j ,i ,:] = tapper * ndimage.gaussian_filter1d(pf.aSR[j, i, :], time_smoothing_kernal_as_dt_fraction * self.dt / pf.dt, mode="constant")
                    
                    # with a cubic interpolation, the interpolated slip1 maybe neg. which does not make sense
                    if pf.slip1[j,i] < 0:
                        pf.aSR[j, i, :] = 0
                        continue
                    # should be the SR
                    integral_STF = np.trapz(np.abs(pf.aSR[j, i, :]), dx=pf.dt)
                    if abs(integral_STF) > 0:
                        pf.aSR[j,i,:] = pf.slip1[j,i] * aSR[j,i, :] / integrawl_STF
        return pf

    
    def area_correction(self, proj):
        """ 
        Kinematic models assume equal fault area for each nodes. However, 
        self.PSarea_cm2 may slightly differ from the patch area from the fault
        geometry (e.g. due to the projection)
        Therefore, we need to update slip to keep seismic potency (area*slip) 
        unchanged
        """
        cm2m = 0.01
        km2m = 1e3 
        PSarea_m2 = self.PSarea_cm2 * cm2m * cm2m
        self.compute_xy_from_latlon(proj)
        nx, ny = self.nx, self.ny
        dip_rad = self.dip*np.pi/180

        
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
        
        # consider dip angle variation in 3D geometry
        factor_area = dx[:, :] * dy[:, :] / (PSarea_m2*np.cos(dip_rad))
        slip1 = self.slip1 * factor_area
        slip2 = self.slip2 * factor_area
        
        print(
            f"done correcting slip for area. \
The correcting factor ranges between {np.amin(factor_area)} and {np.amax(factor_area)}"
        )   
        
        ratio_potency1 = np.sum(np.abs(slip1)) * self.PSarea_cm2/ (np.sum(np.abs(self.slip1)) * self.PSarea_cm2)
        ratio_potency2 = np.sum(np.abs(slip2)) * self.PSarea_cm2/ (np.sum(np.abs(self.slip2)) * self.PSarea_cm2)
        
        print(f"seismic potency ratio - slip direction 1(area correction): {ratio_potency1}" )
        print(f"seismic potency ratio - slip direction 2(area correction): {ratio_potency2}" )
        return slip1, slip2, factor_area
    
    
    
    def generate_netcdf_fl33(self, fname, mesh_size, write_paraview=False):
        """
        Generate netcdf file to be used with Seissol friction law 33
        - The Asagi requires structured mesh for the input. The following interpoloate
            the fault model to a regular mesh with given size - mesh_size
        :parameter:
            fname: string, file name tag for the output files
            mesh_size: [nx, ny],  structural mesh dimension for ASAGI 
            write_paraview: Boolean, files readable for paraview
        """
        [nx, ny] = mesh_size
        xyzc = np.array([self.x.flatten(), self.y.flatten()]).T
        
        X, Y, Xg, Yg = compute_grid(xyzc, nx, ny)
        dist_to_cloud = compute_distance_to_cloud(Xg, Yg, xyzc)
        
        gDip = interpolate_3D(xyzc, Xg, Yg, dist_to_cloud, self.dip.flatten()[:, None]) * np.pi/180
        gStrike = interpolate_3D(xyzc, Xg, Yg, dist_to_cloud, self.strike.flatten()[:, None]) * np.pi/180

        gSx, gSy, gSz = np.sin(gStrike), np.cos(gStrike), np.zeros(gStrike.shape)      
        gDx, gDy, gDz = np.sin(gStrike)*np.cos(gDip), np.cos(gStrike)*np.cos(gDip), np.sin(gDip)

        
        gSls = interpolate_3D(xyzc, Xg, Yg, dist_to_cloud, self.slip1.flatten()[:, None])
        gSld = interpolate_3D(xyzc, Xg, Yg, dist_to_cloud, self.slip2.flatten()[:, None])
        
        # save output
        self.gSls, self.gSld = gSls, gSld
        self.Xg, self.Yg = Xg, Yg
        
        if not os.path.isdir('FL33'):
            os.mkdir('FL33')
        # Write netcdf
        # variables requires 1D array
        writeNetcdf(
            f"FL33/Sls_{fname}",
            [X, Y],
            ["strike_slip"],
            [gSls],
            paraview_readable=write_paraview,
        )
        writeNetcdf(
            f"FL33/Sld_{fname}",
            [X, Y],
            ["dip_slip"],
            [gSld],
            paraview_readable=write_paraview,
        )
        writeNetcdf(
            f"FL33/Asagi_{fname}",
            [X, Y],
            ["strike_slip","dip_slip"],
            [gSls, gSls],
            paraview_readable=False,
        )
        """
        writeNetcdf(
            f"FL33/Asagi_{fname}",
            [X, Y],
            ["strike_slip","dip_slip", \
             "strike_x", "strike_y", "strike_z", \
             "dip_x", "dip_y", "dip_z"],
            [gSls, gSls,  gSx, gSy, gSz, gDx, gDy, gDz],
            paraview_readable=False,
        )
        """
        return None

    def generate_fault_yaml_fl33(self, prefix, tag):
        """Create fault yaml for FL33/34 
        """
        template_yaml = f"""!Switch
[strike_slip, dip_slip]: !EvalModel
    parameters: [strike_slip, dip_slip]
    model: !Switch
        [strike_slip, dip_slip]: !AffineMap
              matrix:
                ua: [1.0, 0.0, 0.0]
                ub: [0.0, 1.0, 0.0]
              translation:
                ua: 0
                ub: 0
              components: !Any
                - !ASAGI
                    file: Asagi_{tag}.nc
                    parameters: [strike_slip, dip_slip]
                    var: data
                    interpolation: linear
                - !ConstantMap
                  map:
                    strike_slip: 0.0
                    dip_slip:    0.0
    components: !FunctionMap
       map:
          #Note the minus on strike_slip to acknowledge the different convention of SeisSol (T_s>0 means right-lateral)
          strike_slip: return -strike_slip;
          dip_slip: return dip_slip;
          rupture_onset: return 0.0;
          tau_S: return 0.025/1.27; 
          tau_R: return 0.1 - 2.*0.025/1.27;
          rutpure_rise_time: return 0.1;
        """
        if not os.path.isdir('FL33'):
            os.mkdir('FL33')
        fname = f"FL33/{prefix}_{tag}_fault.yaml"
        with open(fname, "w") as fid:
            fid.write(template_yaml)
        print(f"done writing {fname}")
        
        return None
