#!/usr/bin/env python3
import vtk
from vtk.util import numpy_support
import numpy as np
from netCDF4 import Dataset
import seissolxdmf
import time
import argparse
from writeNetcdf import writeNetcdf


def compute_lIdt(ndt: int, args_idt: list) -> list:
    """Compute list of time steps to process"""
    if args_idt[0] == -1:
        return list(range(0, ndt))
    else:
        return args_idt


class SeisSolXdmfExtended(seissolxdmf.seissolxdmf):
    def generate_vtk_ugrid(self) -> vtk.vtkUnstructuredGrid:
        """Generate VTK object from SeisSol XDMF file"""
        connect = self.ReadConnect()
        nElements, ndim2 = connect.shape

        grid_type = {3: vtk.vtkPolyData, 4: vtk.vtkUnstructuredGrid}[ndim2]
        grid = grid_type()

        xyz = self.ReadGeometry()
        points = vtk.vtkPoints()
        if ndim2 == 3:
            print("Surface output, assuming the grid is at z=0")
            xyz[:, 2] = 0.0
        points.SetData(numpy_support.numpy_to_vtk(xyz))
        grid.SetPoints(points)

        vtkCells = vtk.vtkCellArray()
        connect2 = np.zeros((nElements, ndim2 + 1), dtype=np.int64)
        connect2[:, 0] = ndim2
        connect2[:, 1:] = connect
        vtkCells.SetCells(nElements, numpy_support.numpy_to_vtkIdTypeArray(connect2))

        if ndim2 == 4:
            grid.SetCells(vtk.VTK_TETRA, vtkCells)
        elif ndim2 == 3:
            grid.SetPolys(vtkCells)
        else:
            raise NotImplementedError
        return grid


def main():
    parser = argparse.ArgumentParser(
        description="Read paraview output, and project on structured grid, generating a netcdf file"
    )
    parser.add_argument("filename", help="SeisSol output filename (xdmf)")
    parser.add_argument(
        "--add2filename",
        help="string to append to filename of the new file",
        type=str,
        default="",
    )
    parser.add_argument(
        "--box",
        nargs=1,
        metavar=("box_arguments"),
        help=(
            "structured grid, "
            "dx x0 x1 y0 y1 z0 z1 for 3D grid, "
            "dx x0 x1 y0 y1 for 2D grid"
        ),
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
        "--variable",
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

    sx = SeisSolXdmfExtended(args.filename)
    ugrid = sx.generate_vtk_ugrid()

    lidt = compute_lIdt(sx.ndt, args.idt)

    if len(args.box[0].split()) == 7:
        dx, x0, x1, y0, y1, z0, z1 = [float(v) for v in args.box[0].split()]
        z = np.arange(z0, z1 + dx, dx)
        assert isinstance(ugrid, vtk.vtkUnstructuredGrid)
        is2DGrid = False
    elif len(args.box[0].split()) == 5:
        dx, x0, x1, y0, y1 = [float(v) for v in args.box[0].split()]
        z = np.array([0])
        z0 = 0.0
        is2DGrid = True
    else:
        raise ValueError(f"wrong number of arguments in args.box {args.box}")

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
        for var in args.variable:
            print("update scalars")
            scalars = vtk.vtkFloatArray()

            W = sx.ReadData(var, idt)
            scalars = numpy_support.numpy_to_vtk(
                num_array=W, deep=True, array_type=vtk.VTK_FLOAT
            )
            ugrid.GetCellData().SetScalars(scalars)

            # Create the CellDataToPointData filter
            cellToPointFilter = vtk.vtkCellDataToPointData()
            cellToPointFilter.SetInputData(ugrid)
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
            for i, sdata in enumerate(args.variable):
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
                args.variable,
                probedData,
                False,
            )


if __name__ == "__main__":
    main()
