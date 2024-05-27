#!/usr/bin/env python3
from sklearn.decomposition import PCA
from netCDF4 import Dataset
from scipy.interpolate import griddata
import seissolxdmf
import numpy as np
from asagiwriter import writeNetcdf
from scipy.ndimage import gaussian_filter
import os
from typing import List, Tuple, Optional
import argparse

class SeissolxdmfExtended(seissolxdmf.seissolxdmf):
    def ReadFaultTags(self) -> np.ndarray:
        """Read fault tag array"""
        return self.Read1dData("fault-tag", self.nElements, isInt=True).T


class AffineMap:
    def __init__(self, ua: np.ndarray, ub: np.ndarray, ta: np.ndarray, tb: np.ndarray):
        """
        Initialize AffineMap with matrix and translation arguments.

        Parameters:
        ua (np.ndarray): Matrix argument for first axis.
        ub (np.ndarray): Matrix argument for second axis.
        ta (np.ndarray): Translation argument for first axis.
        tb (np.ndarray): Translation argument for second axis.
        """
        self.ua = ua
        self.ub = ub
        self.ta = ta
        self.tb = tb


class Grid2D:
    def __init__(self, u: np.ndarray, v: np.ndarray):
        """
        Initialize Grid2D with u and v matrices.

        Parameters:
        u (np.ndarray): U-axis values.
        v (np.ndarray): V-axis values.
        """
        self.u = u
        self.v = v
        self.ug, self.vg = np.meshgrid(u, v)


def gridto2Dlocal(
    sx: SeissolxdmfExtended,
    coords: np.ndarray,
    lengths: List[float],
    dx: float,
    myAffineMap: AffineMap,
    ldataName: List[str],
    ids: np.ndarray,
    gaussian_kernel: Optional[List[float]],
    taper: Optional[List[float]],
) -> Tuple[Grid2D, List[np.ndarray]]:
    """
    Project fault coordinates to 2D local coordinate system and grid data.

    Parameters:
    sx: seissolxdmf reader,
    coords (np.ndarray): Fault coordinates.
    lengths (List[float]): Lengths for grid calculation.
    dx (float): Grid spacing.
    myAffineMap (AffineMap): AffineMap object containing matrix and translation data.
    ldataName (List[str]): List of data names to be processed.
    ids (np.ndarray): Indices of the data points.
    gaussian_kernel (Optional[List[float]]): Gaussian kernel for smoothing.
    taper (Optional[List[float]]): Taper values for clipping data.

    Returns:
    Tuple[Grid2D, List[np.ndarray]]: Gridded 2D local coordinate system and corresponding data.
    """
    xa = np.dot(coords, myAffineMap.ua) + myAffineMap.ta
    xb = np.dot(coords, myAffineMap.ub) + myAffineMap.tb
    xab = np.vstack((xa, xb)).T

    u = np.arange(min(0.0, np.amin(xa)), max(lengths[0], np.amax(xa)) + dx, dx)
    v = np.arange(min(0.0, np.amin(xb)), max(lengths[1], np.amax(xb)) + dx, dx)
    mygrid = Grid2D(u, v)

    ndt = sx.ReadNdt()

    lgridded_myData = []
    for dataName in ldataName:
        # Read Data
        myData = sx.ReadData(dataName, ndt - 1)[ids]
        # grid data and tapper to 30MPa
        gridded_myData = griddata(xab, myData, (mygrid.ug, mygrid.vg), method="nearest")
        gridded_myData_lin = griddata(
            xab, myData, (mygrid.ug, mygrid.vg), method="linear", fill_value=np.nan
        )
        # using linear interpolation when possible, else nearest neighbor
        ids_in = ~np.isnan(gridded_myData_lin)
        gridded_myData[ids_in] = gridded_myData_lin[ids_in]
        if gaussian_kernel:
            gridded_myData = gaussian_filter(gridded_myData, sigma=gaussian_kernel / dx)

        if taper:
            taper_value = taper * 1e6
            gridded_myData[gridded_myData > taper_value] = taper_value
            gridded_myData[gridded_myData < -taper_value] = -taper_value

        plot_data = False
        if plot_data and dataName == ldataName[0]:
            import matplotlib.pyplot as plt

            plt.pcolormesh(mygrid.ug, mygrid.vg, gridded_myData)
            plt.colorbar()
            plt.axis("equal")
            plt.show()
        lgridded_myData.append(gridded_myData)

    return mygrid, lgridded_myData


def writeAllNetcdf(
    mygrid: Grid2D,
    lgridded_myData: List[np.ndarray],
    sName: str,
    ldataName: List[str],
    paraview_readable: bool = False,
) -> None:
    """
    Write gridded data to NetCDF files.

    Parameters:
    mygrid (Grid2D): Gridded 2D local coordinate system.
    lgridded_myData (List[np.ndarray]): List of gridded data arrays.
    sName (str): prefix for the output files.
    ldataName (List[str]): List of data names.
    paraview_readable (bool): Whether to make the NetCDF files ParaView readable.
    """
    if paraview_readable:
        for i, var in enumerate(ldataName):
            writeNetcdf(
                f"{sName}_{var}",
                [mygrid.u, mygrid.v],
                [var],
                [lgridded_myData[i]],
                paraview_readable=True,
            )
    writeNetcdf(
        f"{sName}_TsTdTn",
        [mygrid.u, mygrid.v],
        ldataName,
        lgridded_myData,
        paraview_readable=False,
    )


def generate_input_files(
    fault_filename: str,
    dx: float,
    gaussian_kernel: Optional[float] = None,
    taper: Optional[float] = None,
    paraview_readable: bool = False,
) -> None:
    """
    Generate input files for the given fault data.

    Parameters:
    fault_filename (str): Filename of the fault data.
    dx (float): Grid spacing.
    gaussian_kernel (Optional[float]): Gaussian kernel for smoothing.
    taper (Optional[float]): Taper values for clipping data.
    paraview_readable (bool): Whether to make the NetCDF files ParaView readable.
    """
    if not os.path.exists("ASAGI_files"):
        os.makedirs("ASAGI_files")

    # Compute fault centroids
    sx = SeissolxdmfExtended(fault_filename)
    tags = sx.ReadFaultTags()
    xyz = sx.ReadGeometry()
    connect = sx.ReadConnect()
    faultCentroids = (1.0 / 3.0) * (
        xyz[connect[:, 0]] + xyz[connect[:, 1]] + xyz[connect[:, 2]]
    )

    unique_tags = np.unique(tags)

    template_yaml = f"""!Any
    components:
    """

    for tag in unique_tags:
        ids = np.where(tags == tag)[0]
        connect_selected = connect[ids, :]
        selected_vertex = list(set(connect_selected.flatten()))
        xyz_selected = xyz[selected_vertex, :]

        # Perform PCA to get principal axes
        pca = PCA(n_components=2)
        points = pca.fit_transform(xyz_selected)
        la, lb = np.amax(points, axis=0) - np.amin(points, axis=0)
        ua, ub = pca.components_
        lower_left = np.argmin(np.sum(points, axis=1))
        xu1 = xyz_selected[lower_left]

        ta = -np.dot(xu1, ua)
        tb = -np.dot(xu1, ub)

        template_yaml += f""" - !GroupFilter
        groups: {tag}
        components: !AffineMap
          matrix:
            ua: [{ua[0]}, {ua[1]}, {ua[2]}]
            ub: [{ub[0]}, {ub[1]}, {ub[2]}]
          translation:
            ua: {ta}
            ub: {tb}
          components: !Any
            - !ASAGI
                file: ASAGI_files/fault{tag}_TsTdTn.nc
                parameters: [Ts0, Td0, Pn0]
                var: data
                interpolation: linear
            - !ConstantMap
              map:
                Ts0: 0.0
                Td0: 0.0\n"""
        myAffineMap = AffineMap(ua, ub, ta, tb)
        ldataName = ["Ts0", "Td0", "Pn0"]
        grid, lgridded_myData = gridto2Dlocal(
            sx,
            faultCentroids[ids],
            [la, lb],
            dx,
            myAffineMap,
            ldataName,
            ids,
            gaussian_kernel,
            taper,
        )

        fn = f"fault{tag}"
        ldataName = ["T_s", "T_d", "T_n"]
        writeAllNetcdf(
            grid, lgridded_myData, f"ASAGI_files/{fn}", ldataName, paraview_readable
        )

    fname = "yaml_files/Ts0Td0.yaml"
    with open(fname, "w") as fid:
        fid.write(template_yaml)
    print(f"done writing {fname}")


def main() -> None:
    """
    Main function to parse arguments and generate input files
    by projecting 3D fault output onto 2D grids for ASAGI.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Project 3D fault output onto 2D grids to be read with ASAGI. "
            "One grid per fault tag."
        )
    )
    parser.add_argument("fault_filename", help="Fault.xdmf filename")
    parser.add_argument(
        "--dx",
        nargs=1,
        help="Grid sampling",
        type=float,
        default=[100.0],
    )
    parser.add_argument(
        "--gaussian_kernel",
        metavar="sigma_m",
        nargs=1,
        help="Apply a Gaussian kernel to smooth out input stresses",
        type=float,
    )
    parser.add_argument(
        "--taper",
        nargs=1,
        help="Taper stress value (MPa)",
        type=float,
    )
    parser.add_argument(
        "--paraview_readable",
        dest="paraview_readable",
        action="store_true",
        help="Write NetCDF files readable by ParaView",
        default=False,
    )

    args = parser.parse_args()
    generate_input_files(
        args.fault_filename,
        args.dx[0],
        args.gaussian_kernel[0] if args.gaussian_kernel else None,
        args.taper[0] if args.taper else None,
        args.paraview_readable,
    )


if __name__ == "__main__":
    main()
