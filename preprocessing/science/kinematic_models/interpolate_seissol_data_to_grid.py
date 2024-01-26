#!/usr/bin/env python3
import vtk
from vtk.util import numpy_support
import numpy as np
from netCDF4 import Dataset
import seissolxdmf
import time
import argparse


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
            # ldata8 = [(name, "f8") for name in lName]
            mattype4 = np.dtype(ldata4)
            # mattype8 = np.dtype(ldata8)
            mat_t = rootgrp.createCompoundType(mattype4, "material")

            # this transform the nD array into an array of tuples
            arr = np.stack([lData[i] for i in range(len(lName))], axis=len(dims))
            newarr = arr.view(dtype=mattype4)
            newarr = newarr.reshape(newarr.shape[:-1])
            mat = rootgrp.createVariable("data", mat_t, dims)
            mat[:] = newarr


def compute_lIdt(ndt, args_idt):
    if args_idt[0] == -1:
        lIdt = list(range(0, ndt))
    else:
        lIdt = args_idt
        for idt in args.idt:
            if (idt < 0) | (idt >= ndt):
                raise Exception(f"idt {idt} not in range(0,{ndt})")
    return lIdt


class seissolxdmfExtended(seissolxdmf.seissolxdmf):
    def generateVtkObject(self):
        """Filling in vtk arrays with data from hdf5 file."""

        connect = self.ReadConnect()
        nElements, ndim2 = connect.shape

        xyz = self.ReadGeometry()
        points = vtk.vtkPoints()
        if ndim2 == 3:
            print("surface output, assuming the grid is at z=0")
            xyz[:, 2] = 0.0
        points.SetData(numpy_support.numpy_to_vtk(xyz))

        vtkCells = vtk.vtkCellArray()
        connect2 = np.zeros((nElements, ndim2 + 1), dtype=np.int64)
        # number of points in the cell
        connect2[:, 0] = ndim2
        connect2[:, 1:] = connect
        vtkCells.SetCells(nElements, numpy_support.numpy_to_vtkIdTypeArray(connect2))

        if ndim2 == 4:
            unstrGrid3d = vtk.vtkUnstructuredGrid()
            unstrGrid3d.SetPoints(points)
            unstrGrid3d.SetCells(vtk.VTK_TETRA, vtkCells)
            return unstrGrid3d
        elif ndim2 == 3:
            myPolydata = vtk.vtkPolyData()
            myPolydata.SetPoints(points)
            myPolydata.SetPolys(vtkCells)
            return myPolydata
        else:
            raise NotImplementedError


parser = argparse.ArgumentParser(
    description="Read paraview output, and project on structured grid, generating a netcdf file"
)
parser.add_argument("filename", help="SeisSol surface output filename (xdmf)")
parser.add_argument(
    "--add2filename",
    help="string to append to filename of the new file",
    type=str,
    default="",
)
parser.add_argument(
    "--box",
    nargs=1,
    metavar=("dx x0 x1 y0 y1 z0 z1"),
    help="structured grid",
    required=True,
)
parser.add_argument(
    "--idt",
    nargs="+",
    help="list of time step to process (1st = 0); -1 = all",
    type=int,
    default=([-1]),
)
parser.add_argument(
    "--Data",
    nargs="+",
    required=True,
    metavar=("variable"),
    help="name of variable; all for all stored quantities",
)
parser.add_argument(
    "--paraview_readable",
    dest="paraview_readable",
    action="store_true",
    help="generate paraview readable netcdf file",
)

args = parser.parse_args()

sx = seissolxdmfExtended(args.filename)
unstrGrid3d = sx.generateVtkObject()

lidt = compute_lIdt(sx.ndt, args.idt)

if len(args.box[0].split()) == 7:
    dx, x0, x1, y0, y1, z0, z1 = [float(v) for v in args.box[0].split()]
    z = np.arange(z0, z1 + dx, dx)
    assert isinstance(unstrGrid3d, vtk.vtkUnstructuredGrid)
    is2DGrid = False
elif len(args.box[0].split()) == 5:
    dx, x0, x1, y0, y1 = [float(v) for v in args.box[0].split()]
    z = np.array([0])
    z0 = 0.0
    # assert(isinstance(unstrGrid3d, vtk.vtkPolyData))
    is2DGrid = True
else:
    raise f"wrong number of arguments in args.box {args.box}"

x = np.arange(x0, x1 + dx, dx)
y = np.arange(y0, y1 + dx, dx)

if is2DGrid:
    xx, yy = np.meshgrid(x, y)
else:
    yy, zz, xx = np.meshgrid(y, z, x)

# Create grid image volume
image1Size = [x.shape[0], y.shape[0], z.shape[0]]
image1Origin = [x0, y0, z0]
image1Spacing = [dx, dx, dx]

imageData1 = vtk.vtkImageData()
imageData1.SetDimensions(image1Size)
imageData1.SetOrigin(image1Origin)
imageData1.SetSpacing(image1Spacing)


# Perform the interpolation
probeFilter = vtk.vtkProbeFilter()

probeFilter.SetInputData(imageData1)
probeFilter.SpatialMatchOn()

for idt in lidt:
    probedData = []
    for var in args.Data:
        print("update scalars")
        scalars = vtk.vtkFloatArray()

        W = sx.ReadData(var, idt)
        scalars = numpy_support.numpy_to_vtk(
            num_array=W, deep=True, array_type=vtk.VTK_FLOAT
        )
        unstrGrid3d.GetCellData().SetScalars(scalars)
        # Create the CellDataToPointData filter
        cellToPointFilter = vtk.vtkCellDataToPointData()
        cellToPointFilter.SetInputData(unstrGrid3d)
        cellToPointFilter.Update()

        # Get the output grid with point data
        outputGrid = cellToPointFilter.GetOutput()

        # Perform the interpolation
        probeFilter.SetSourceData(outputGrid)
        start = time.time()
        print("start prob filter")
        probeFilter.Update()
        stop = time.time()
        print(f"{var} {idt}: done prob filter in {stop - start} s")

        polyout = probeFilter.GetOutput()
        projData = polyout.GetPointData().GetScalars()
        projDataNp = numpy_support.vtk_to_numpy(projData).reshape(xx.shape)
        probedData.append(projDataNp)

    xyz = [x, y] if is2DGrid else [x, y, z]
    if args.paraview_readable:
        for i, sdata in enumerate(args.Data):
            writeNetcdf(
                f"gridded_{sdata}_{dx:.0f}_{idt}{args.add2filename}",
                xyz,
                [sdata],
                [probedData[i]],
                paraview_readable=True,
            )
    else:
        writeNetcdf(
            f"gridded_asagi_dx{dx:.0f}_{idt}{args.add2filename}",
            xyz,
            args.Data,
            probedData,
            False,
        )
