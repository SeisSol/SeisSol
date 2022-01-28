import os
import os.path
import seissolxdmf
import shutil
import recreateXdmf
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="resample output file and write as binary files")
parser.add_argument("xdmfFilename", help="xdmf output file")
parser.add_argument("--add2prefix", help="string to append to prefix for new file", type=str, default="_resampled")
parser.add_argument("--Data", nargs="+", metavar=("variable"), default=(""), help="Data to resample (example SRs)")
parser.add_argument("--downsample", help="write one out of n output", type=int)
parser.add_argument("--precision", type=str, choices=["float", "double"], default="float", help="precision of output file")
parser.add_argument("--backend", type=str, choices=["hdf5", "raw"], default="hdf5", help="backend used: raw (.bin file), hdf5 (.h5)")
parser.add_argument("--last", dest="last", default=False, action="store_true", help="output last time step")
parser.add_argument("--idt", nargs="+", help="list of time step to write (ex $(seq 7 3 28))", type=int)
parser.add_argument(
    "--xfilter",
    nargs=2,
    metavar=("xmin", "xmax"),
    help="output only cells with x center coordinate in range xmin xmax",
    type=float,
)
parser.add_argument(
    "--yfilter",
    nargs=2,
    metavar=("ymin", "ymax"),
    help="output only cells with y center coordinate in range ymin ymax",
    type=float,
)
parser.add_argument(
    "--zfilter",
    nargs=2,
    metavar=("zmin", "zmax"),
    help="output only cells with z center coordinate in range zmin zmax",
    type=float,
)
args = parser.parse_args()


class seissolxdmfExtended(seissolxdmf.seissolxdmf):
    def ReadData(self, dataName, idt=-1):
        if dataName == "SR":
            SRs = super().ReadData("SRs", idt)
            SRd = super().ReadData("SRd", idt)
            return np.sqrt(SRs ** 2 + SRd ** 2)
        else:
            return super().ReadData(dataName, idt)


sx = seissolxdmfExtended(args.xdmfFilename)
xyz = sx.ReadGeometry()
connect = sx.ReadConnect()

spatial_filtering = (args.xfilter or args.yfilter) or args.zfilter
nElements = connect.shape[0]

if spatial_filtering:
    print("Warning: spatial filtering significantly slows down this script")
    ids = range(0, sx.nElements)
    xyzc = (xyz[connect[:, 0], :] + xyz[connect[:, 1], :] + xyz[connect[:, 2], :]) / 3.0

    def filter_cells(coords, filter_range):
        m = 0.5 * (filter_range[0] + filter_range[1])
        d = 0.5 * (filter_range[1] - filter_range[0])
        return np.where(np.abs(coords[:] - m) < d)[0]

    if args.xfilter:
        id0 = filter_cells(xyzc[:, 0], args.xfilter)
        ids = np.intersect1d(ids, id0) if len(ids) else id0
    if args.yfilter:
        id0 = filter_cells(xyzc[:, 1], args.yfilter)
        ids = np.intersect1d(ids, id0) if len(ids) else id0
    if args.zfilter:
        id0 = filter_cells(xyzc[:, 2], args.zfilter)
        ids = np.intersect1d(ids, id0) if len(ids) else id0

    if len(ids):
        connect = connect[ids, :]
        nElements = connect.shape[0]
        if nElements != sx.nElements:
            print(f"extracting {nElements} cells out of {sx.nElements}")
        else:
            spatial_filtering = False
    else:
        raise ValueError("all elements are outside filter range")
ndt = sx.ndt

if args.last:
    indices = [ndt - 1]
    if args.idt or args.downsample:
        print("last option cannot be used together with idt and downsample options")
        exit()
else:
    if args.idt and args.downsample:
        print("idt and downsample options cannot be used together")
        exit()
    elif not args.downsample:
        indices = args.idt
    else:
        indices = range(0, ndt, args.downsample)

# Check if input is in hdf5 format or not
dataLocation, data_prec, MemDimension = sx.GetDataLocationPrecisionMemDimension("partition")
splitArgs = dataLocation.split(":")
if len(splitArgs) == 2:
    isHdf5 = True
else:
    isHdf5 = False

if args.precision == "double":
    myDtype = "float64"
    myprec = 8
else:
    myDtype = "float32"
    myprec = 4

if args.backend == "raw":
    write2Binary = True
else:
    write2Binary = False

prefix = os.path.splitext(args.xdmfFilename)[0]
prefix_new = recreateXdmf.generate_new_prefix(prefix, args.add2prefix)

# Create folders if necessary
if write2Binary:
    if not os.path.exists(prefix_new + "_cell"):
        os.mkdir(prefix_new + "_cell")
    if not os.path.exists(prefix_new + "_cell/mesh0/"):
        os.mkdir(prefix_new + "_cell/mesh0/")
    if not os.path.exists(prefix_new + "_vertex"):
        os.mkdir(prefix_new + "_vertex")
    if not os.path.exists(prefix_new + "_vertex/mesh0/"):
        os.mkdir(prefix_new + "_vertex/mesh0/")


# Write geometry and connect
if write2Binary:
    fn2 = prefix_new + "_cell/mesh0/connect.bin"
    if isHdf5:
        # write Connect
        output_file = open(fn2, "wb")
        connect.tofile(output_file)
        output_file.close()
    else:
        shutil.copy2(os.path.splitext(args.xdmfFilename)[0] + "_cell/mesh0/connect.bin", fn2)
    print("done writing " + fn2)
    # write geometry
    fn3 = prefix_new + "_vertex/mesh0/geometry.bin"
    output_file = open(fn3, "wb")
    xyz.tofile(output_file)
    output_file.close()
    print("done writing " + fn3)
else:
    import h5py

    # write geometry to hdf5 format
    h5fv = h5py.File(prefix_new + "_vertex.h5", "w")
    h5fv.create_dataset("/mesh0/geometry", data=xyz)
    h5fv.close()
    print("done writing " + prefix_new + "_vertex.h5")
    # write connect to hdf5 format
    h5fc = h5py.File(prefix_new + "_cell.h5", "w")
    h5fc.create_dataset("/mesh0/connect", data=connect)


# Write data items
for ida, sdata in enumerate(args.Data):
    if write2Binary:
        fname2 = prefix_new + "_cell/mesh0/" + args.Data[ida] + ".bin"
        output_file = open(fname2, "wb")
    else:
        dset = h5fc.create_dataset("/mesh0/" + args.Data[ida], (len(indices), nElements), dtype=myDtype)
    # read only one row
    print(sdata, end=" ", flush=True)
    for kk, i in enumerate(indices):
        if (kk % 10 == 0) and kk > 0:
            print(kk)
        else:
            print(kk, end=" ", flush=True)
        if i >= ndt:
            print("ignoring index %d>=ndt=%d" % (i, ndt))
            continue
        if spatial_filtering:
            myData = sx.ReadData(args.Data[ida], idt=i)[ids]
        else:
            myData = sx.ReadData(args.Data[ida], idt=i)
        if write2Binary:
            myData.astype(myDtype).tofile(output_file)
        else:
            dset[kk, :] = myData[:]
    if write2Binary:
        output_file.close()
        print("done writing " + fname2)

if not write2Binary:
    h5fc.close()
    print("done writing " + prefix_new + "_cell.h5")

# Now recreate the Xdmf

prefix = os.path.splitext(args.xdmfFilename)[0]

# Read all parameters from the xdmf file of SeisSol (if any)

dt = recreateXdmf.ReadDtFromXdmf(args.xdmfFilename)
nvertex = recreateXdmf.ReadNvertexFromXdmf(args.xdmfFilename)
ndt, nmem = recreateXdmf.ReadNdtNmemFromXdmf(args.xdmfFilename)
recreateXdmf.recreateXdmf(prefix, prefix_new, nvertex, nElements, nElements, dt, indices, args.Data, not write2Binary, myprec, args.add2prefix)
