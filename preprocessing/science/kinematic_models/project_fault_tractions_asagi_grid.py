#!/usr/bin/env python3
from sklearn.decomposition import PCA
from netCDF4 import Dataset
from scipy.interpolate import griddata
import seissolxdmf
import argparse
import numpy as np
from writeNetcdf import writeNetcdf
import os


class seissolxdmfExtended(seissolxdmf.seissolxdmf):
    def ReadFaultTags(self):
        """Read partition array"""
        return self.Read1dData("clustering", self.nElements, isInt=True).T


class AffineMap:
    def __init__(self, ua, ub, ta, tb):
        # matrix arguments
        self.ua = ua
        self.ub = ub
        # translation arguments
        self.ta = ta
        self.tb = tb


class Grid2D:
    def __init__(self, u, v):
        # matrix arguments
        self.u = u
        self.v = v
        # translation arguments
        self.ug, self.vg = np.meshgrid(u, v)


def Gridto2Dlocal(coords, lengths, myAffineMap, fault_fname, ldataName, ids):
    # project fault coordinates to 2D local coordinate system
    xa = np.dot(coords, myAffineMap.ua) + myAffineMap.ta
    xb = np.dot(coords, myAffineMap.ub) + myAffineMap.tb
    xab = np.vstack((xa, xb)).T

    u = np.arange(min(0.0, np.amin(xa)), max(lengths[0], np.amax(xa)) + dx, dx)
    v = np.arange(min(0.0, np.amin(xb)), max(lengths[1], np.amax(xb)) + dx, dx)
    mygrid = Grid2D(u, v)

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

        if args.taper:
            taper_value = args.taper[0] * 1e6
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


def WriteAllNetcdf(mygrid, lgridded_myData, sName, ldataName):
    """
    for i, var in enumerate(ldataName):
        writeNetcdf(
            f"{sName}_{var}",
            [mygrid.u, mygrid.v],
            [var],
            [lgridded_myData[i]],
            paraview_readable=True,
        )
    """
    writeNetcdf(
        f"{sName}_TsTdTn",
        [mygrid.u, mygrid.v],
        ldataName,
        lgridded_myData,
        paraview_readable=False,
    )


# parsing python arguments
parser = argparse.ArgumentParser(
    description=(
        "project 3d fault output onto 2d grids to be read with Asagi. One grid per"
        " fault tag"
    )
)
parser.add_argument("fault", help="fault.xdmf filename")
parser.add_argument(
    "--dx",
    nargs=1,
    help="grid smapling",
    type=float,
    default=([50]),
)
parser.add_argument(
    "--taper",
    nargs=1,
    help="tapper stress value (MPa)",
    type=float,
)

args = parser.parse_args()
dx = args.dx[0]

if not os.path.exists("ASAGI_files"):
    os.makedirs("ASAGI_files")

# Compute fault centroids
sx = seissolxdmfExtended(args.fault)
tags = sx.ReadFaultTags()
xyz = sx.ReadGeometry()
connect = sx.ReadConnect()
faultCentroids = (1.0 / 3.0) * (
    xyz[connect[:, 0]] + xyz[connect[:, 1]] + xyz[connect[:, 2]]
)

unique_tags = np.unique(tags)
# Read ndt and nElements
ndt = sx.ReadNdt()

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
                Td0: 0.0
"""
    myAffineMap = AffineMap(ua, ub, ta, tb)
    ldataName = ["Ts0", "Td0", "Pn0"]
    grid, lgridded_myData = Gridto2Dlocal(
        faultCentroids[ids], [la, lb], myAffineMap, args.fault, ldataName, ids
    )

    fn = f"fault{tag}"
    ldataName = ["T_s", "T_d", "T_n"]
    WriteAllNetcdf(grid, lgridded_myData, f"ASAGI_files/{fn}", ldataName)

fname = f"Ts0Td0.yaml"
with open(fname, "w") as fid:
    fid.write(template_yaml)
print(f"done writing {fname}")
