import numpy as np
from netCDF4 import Dataset


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
        sdimVarNames = "abcdef"
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
            mattype4 = np.dtype(ldata4)
            mat_t = rootgrp.createCompoundType(mattype4, "material")

            # this transform the n D array into an array of tuples
            arr = np.stack(
                [lData[i].astype(np.float32) for i in range(len(lName))], axis=len(dims)
            )
            arr = np.ascontiguousarray(arr)
            newarr = arr.view(dtype=mattype4)
            newarr = newarr.reshape(newarr.shape[:-1])
            mat = rootgrp.createVariable("data", mat_t, dims)
            mat[:] = newarr
