import os
import os.path
import seissolxdmf
import shutil
from recreateXdmf import *

import argparse

parser = argparse.ArgumentParser(description="resample output file and write as binary files")
parser.add_argument("xdmfFilename", help="xdmf output file")
parser.add_argument("--add2prefix", help="string to append to prefix for new file", type=str, default="_resampled")
parser.add_argument("--Data", nargs="+", metavar=("variable"), default=(""), help="Data to resample (example SRs)")
parser.add_argument("--downsample", help="write one out of n output", type=int)
parser.add_argument("--float", dest="float", default=False, action="store_true", help="convert to float")
parser.add_argument("--hdf5", dest="hdf5", default=False, action="store_true", help="output as hdf5")
parser.add_argument("--idt", nargs="+", help="list of time step to write (ex $(seq 7 3 28))", type=int)
args = parser.parse_args()

sx = seissolxdmf.seissolxdmf(args.xdmfFilename)
connect = sx.ReadConnect()
nElements = sx.nElements
ndt = sx.ndt

if (args.idt != None) and (args.downsample != None):
    print("idt and downsample options cannot be used together")
    exit()
elif args.downsample == None:
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

if args.float:
    myDtype = "float32"
    myprec = 4
else:
    myDtype = "float64"
    myprec = 8

if args.hdf5:
    write2Binary = False
else:
    write2Binary = True

prefix = os.path.splitext(args.xdmfFilename)[0]
prefix_new = generate_new_prefix(prefix, args.add2prefix)

############Create folders if necessary#############
if write2Binary:
    if not os.path.exists(prefix_new + "_cell"):
        os.mkdir(prefix_new + "_cell")
    if not os.path.exists(prefix_new + "_cell/mesh0/"):
        os.mkdir(prefix_new + "_cell/mesh0/")
    if not os.path.exists(prefix_new + "_vertex"):
        os.mkdir(prefix_new + "_vertex")
    if not os.path.exists(prefix_new + "_vertex/mesh0/"):
        os.mkdir(prefix_new + "_vertex/mesh0/")


#############write geometry and connect#############
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
    geometry = sx.ReadGeometry()
    fn3 = prefix_new + "_vertex/mesh0/geometry.bin"
    output_file = open(fn3, "wb")
    geometry.tofile(output_file)
    output_file.close()
    print("done writing " + fn3)
else:
    import h5py
    # write geometry to hdf5 format
    h5fv = h5py.File(prefix_new + "_vertex.h5", "w")
    geometry = sx.ReadGeometry()
    h5fv.create_dataset("/mesh0/geometry", data=geometry)
    h5fv.close()
    print("done writing " + prefix_new + "_vertex.h5")
    # write connect to hdf5 format
    h5fc = h5py.File(prefix_new + "_cell.h5", "w")
    h5fc.create_dataset("/mesh0/connect", data=connect)


#############write data items#######################
for ida, sdata in enumerate(args.Data):
    if write2Binary:
        fname2 = prefix_new + "_cell/mesh0/" + args.Data[ida] + ".bin"
        output_file = open(fname2, "wb")
    else:
        dset = h5fc.create_dataset("/mesh0/" + args.Data[ida], (len(indices), nElements), dtype=myDtype)
    # read only one row
    for kk, i in enumerate(indices):
        if i >= ndt:
            print("ignoring index %d>=ndt=%d" % (i, ndt))
            continue
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

#####Now recrete the Xdmf ##################

prefix = os.path.splitext(args.xdmfFilename)[0]

### Read all parameters from the xdmf file of SeisSol (if any)

ncells = ReadNcellsFromXdmf(args.xdmfFilename)
dt = ReadDtFromXdmf(args.xdmfFilename)
nvertex = ReadNvertexFromXdmf(args.xdmfFilename)
ndt, nmem = ReadNdtNmemFromXdmf(args.xdmfFilename)
recreateXdmf(prefix, prefix_new, nvertex, ncells, ncells, dt, indices, args.Data, not write2Binary, myprec, args.add2prefix)
