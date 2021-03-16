import numpy as np
from netCDF4 import Dataset


def writeNetcdf4Paraview(sname, x, y, aName, aData):
    "create a netcdf file readable by paraview (but not by ASAGI)"
    fname = sname + "_paraview.nc"
    print("writing " + fname)
    ####Creating the netcdf file
    nx = x.shape[0]
    ny = y.shape[0]
    rootgrp = Dataset(fname, "w", format="NETCDF4")
    rootgrp.createDimension("u", nx)
    rootgrp.createDimension("v", ny)

    vx = rootgrp.createVariable("u", "f4", ("u",))
    vx[:] = x
    vy = rootgrp.createVariable("v", "f4", ("v",))
    vy[:] = y
    for i in range(len(aName)):
        vTd = rootgrp.createVariable(aName[i], "f4", ("v", "u"))
        vTd[:, :] = aData[i][:, :]
    rootgrp.close()


def writeNetcdf4SeisSol(sname, x, y, aName, aData):
    "create a netcdf file readable by ASAGI (but not by paraview)"
    ########## creating the file for SeisSol
    fname = sname + "_ASAGI.nc"
    print("writing " + fname)
    ####Creating the netcdf file
    nx = x.shape[0]
    ny = y.shape[0]

    rootgrp = Dataset(fname, "w", format="NETCDF4")

    rootgrp.createDimension("u", nx)
    rootgrp.createDimension("v", ny)

    vx = rootgrp.createVariable("u", "f4", ("u",))
    vx[:] = x
    vy = rootgrp.createVariable("v", "f4", ("v",))
    vy[:] = y
    ldata4 = [(name, "f4") for name in aName]
    ldata8 = [(name, "f8") for name in aName]
    mattype4 = np.dtype(ldata4)
    mattype8 = np.dtype(ldata8)
    mat_t = rootgrp.createCompoundType(mattype4, "material")

    # this transform the 4 D array into an array of tuples
    arr = np.stack([aData[i] for i in range(len(aName))], axis=2)
    newarr = arr.view(dtype=mattype8)
    newarr = newarr.reshape(newarr.shape[:-1])
    mat = rootgrp.createVariable("data", mat_t, ("v", "u"))
    mat[:] = newarr
    rootgrp.close()


dx = 200.0
nx = 40
ny = 20
x = np.linspace(0, nx - 1, nx) * dx
y = np.linspace(0, ny - 1, ny) * dx

# Generate a Gaussian shaped slip distribution for illustration purposes
x0, y0 = 3000.0, 1500.0
xg, yg = np.meshgrid(x, y)
d = np.sqrt((xg - x0) ** 2 + (yg - y0) ** 2)
mu = np.sqrt(x0 ** 2 + y0 ** 2)
sigma = 500.0
slip = np.exp(-(d ** 2 / (2.0 * sigma ** 2)))
rake = 40 * np.ones((ny, nx)) * np.pi / 180.0


strike_slip = slip * np.cos(rake)
dip_slip = slip * np.sin(rake)

prefix = "test"
ldataName = ["strike_slip", "dip_slip"]
lgridded_myData = [strike_slip, dip_slip]

writeNetcdf4Paraview(prefix, x, y, ldataName, lgridded_myData)
writeNetcdf4SeisSol(prefix, x, y, ldataName, lgridded_myData)
