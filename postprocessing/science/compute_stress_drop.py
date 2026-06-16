# compute static stress drop from SeisSol fault output

import numpy as np
import seissolxdmf as sx
import argparse

parser = argparse.ArgumentParser(
    description="compute the static stress drop from the fault xdmf file"
)
parser.add_argument("filename", help="path+filename-fault.xdmf")
parser.add_argument(
    "--slip_threshold",
    nargs=1,
    metavar=("slip threshold"),
    default=([0.01]),
    help="slip threshold to define the rupture area",
    type=float,
)
parser.add_argument(
    "--time_indices",
    nargs=2,
    metavar=("t0", "t1"),
    default=([0, -1]),
    help="provide time steps that include the target rupture episode (e.g., when the file contains multiple events)",
    type=float,
)
args = parser.parse_args()
timeIndices = args.time_indices


def StressDropAllCells(sx, timeIndices=timeIndices):
    shearStress_0 = np.sqrt(
        sx.ReadData("Ts0", timeIndices[0]) ** 2
        + sx.ReadData("Td0", timeIndices[0]) ** 2
    )
    shearStress_end = np.sqrt(
        sx.ReadData("Ts0", timeIndices[1]) ** 2
        + sx.ReadData("Td0", timeIndices[1]) ** 2
    )
    stressDropAllCells = shearStress_0 - shearStress_end
    return stressDropAllCells


def ComputeRuptureAreaCells(sx):
    xyz = sx.ReadGeometry()
    connect = sx.ReadConnect()
    cross = np.cross(
        xyz[connect[:, 1], :] - xyz[connect[:, 0], :],
        xyz[connect[:, 2], :] - xyz[connect[:, 0], :],
    )
    area0 = 0.5 * np.apply_along_axis(np.linalg.norm, 1, cross)
    return area0

print("Computing stress drops...")
sx = sx.seissolxdmf(args.filename)
if timeIndices[1] == -1:
    timeIndices[1] = sx.ReadNdt()

ASl = sx.ReadData("ASl", timeIndices[1]) - sx.ReadData("ASl", timeIndices[0])
idx = np.where(ASl > args.slip_threshold[0])
print(
    f"Rupture area is defined as the area where slip exceeds {args.slip_threshold[0]} m"
)

ruptureAreaCells = ComputeRuptureAreaCells(sx, ASl)
ruptureArea = np.sum(ruptureAreaCells[idx])
print(f"Ruptured area: {ruptureArea*1e-6:.3f} km^2")

stressDropAllCells = StressDropAllCells(sx, timeIndices=timeIndices)
averageStressDrop = np.average(stressDropAllCells[idx], weights=ruptureAreaCells[idx])
print(f"Average stress drop: {averageStressDrop*1e-6:.3f} MPa")

slipWeightedStressDrop = np.average(stressDropAllCells, weights=ASl)
print(f"Slip-weighted stress drop: {slipWeightedStressDrop*1e-6:.3f} MPa")
