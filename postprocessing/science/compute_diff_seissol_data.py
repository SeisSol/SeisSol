#!/usr/bin/env python3

import h5py
import numpy as np
import argparse
import os
import seissolxdmf as sx

parser = argparse.ArgumentParser(
    description="make difference between 2 output files: f2-f1"
)
parser.add_argument("xdmf_filename1", help="filename1")
parser.add_argument("xdmf_filename2", help="filename2")
parser.add_argument(
    "--idt",
    nargs="+",
    required=True,
    help="list of time step to differenciate (1st = 0); -1 = all",
    type=int,
)
parser.add_argument(
    "--Data",
    nargs="+",
    required=True,
    metavar=("variable"),
    help="Data to differenciate (example SRs)",
)
parser.add_argument(
    "--ratio",
    dest="ratio",
    default=False,
    action="store_true",
    help="compute relative ratio (f1-f2)/f1 instead of f2-f1",
)

args = parser.parse_args()


def CreateXdmf(fn, nNodes, ntriangles, aDataName, dt, node_per_element):
    topology = "Tetrahedron" if node_per_element == 4 else "Triangle"
    xdmf = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
 <Domain>"""
    for i, idt in enumerate(args.idt):
        xdmf = (
            xdmf
            + f"""
  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">
   <Grid Name="step_{idt}" GridType="Uniform"><!-- mesh id: 0, mesh step: 0 -->
    <Topology TopologyType="{topology}" NumberOfElements="{ntriangles}">
     <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="{ntriangles} {node_per_element}">{fn}.h5:/mesh0/connect</DataItem>
    </Topology>
    <Geometry name="geo" GeometryType="XYZ" NumberOfElements="{nNodes}">
     <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="{nNodes} 3">{fn}.h5:/mesh0/geometry</DataItem>
    </Geometry>
    <Time Value="{idt*dt}"/>"""
        )
        for dataName in aDataName:
            xdmf = (
                xdmf
                + f"""
    <Attribute Name="{dataName}" Center="Cell">
     <DataItem ItemType="HyperSlab" Dimensions="{ntriangles}">
      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">{i} 0 1 1 1 {ntriangles}</DataItem>
      <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="{i+1} {ntriangles}">{fn}.h5:/mesh0/{dataName}</DataItem>
     </DataItem>
    </Attribute>"""
            )
        xdmf = (
            xdmf
            + """
   </Grid>
  </Grid>"""
        )
    xdmf = (
        xdmf
        + """
 </Domain>
</Xdmf>
"""
    )
    fid = open(fn + ".xdmf", "w")
    fid.write(xdmf)
    fid.close()


def write_xdmfh5(fname, aDataName, xyz, connect, aData, dt):
    nNodes = xyz.shape[0]
    ntriangles, node_per_element = connect.shape
    CreateXdmf(fname, nNodes, ntriangles, aDataName, dt, node_per_element)
    # Write h5 file

    h5f = h5py.File(fname + ".h5", "w")
    h5f.create_dataset("mesh0/connect", (ntriangles, node_per_element), dtype="i8")
    h5f["mesh0/connect"][:, :] = connect[:, :]
    h5f.create_dataset("mesh0/geometry", xyz.shape, dtype="d")
    h5f["mesh0/geometry"][:, :] = xyz[:, :]
    for k, dataName in enumerate(aDataName):
        hdname = "mesh0/" + dataName
        h5f.create_dataset(hdname, (len(args.idt), ntriangles), dtype="d")
        for i, idt in enumerate(args.idt):
            h5f[hdname][i, :] = aData[k][idt, :]
    h5f.close()
    print(f"done writing {fname}.xdmf")


def read_reshape2d(sx, dataname):
    """read seissol dataset
    and if there is only one time stamp
    create a second dimension of size 1"""
    myData = sx.ReadData(dataname)
    if len(myData.shape) == 1:
        myData = myData.reshape((1, myData.shape[0]))
    return myData


def fuzzysort(arr, idx, dim=0, tol=1e-6):
    """
    return indexes of sorted points robust to small perturbations of individual components.
    https://stackoverflow.com/questions/19072110/numpy-np-lexsort-with-fuzzy-tolerant-comparisons
    """
    arrd = arr[dim]
    srtdidx = sorted(idx, key=arrd.__getitem__)

    i, ix = 0, srtdidx[0]
    for j, jx in enumerate(srtdidx[1:], start=1):
        if arrd[jx] - arrd[ix] >= tol:
            if j - i > 1:
                srtdidx[i:j] = fuzzysort(arr, srtdidx[i:j], dim + 1, tol)
            i, ix = j, jx

    if i != j:
        srtdidx[i:] = fuzzysort(arr, srtdidx[i:], dim + 1, tol)

    return srtdidx


def lookup_sorted_geom(geom):
    """return the indices to sort the
    geometry array by x, then y, then z
    and the associated inverse look-up table
    """
    ind = fuzzysort(geom.T, list(range(0, geom.shape[0])), tol=1e-4)
    # generate inverse look-up table
    dic = {i: index for i, index in enumerate(ind)}
    ind_inv = np.zeros_like(ind)
    for k, v in dic.items():
        ind_inv[v] = k
    return ind, ind_inv


def read_geom_connect(sx):
    return sx.ReadGeometry(), sx.ReadConnect()


def return_sorted_geom_connect(sx):
    """sort geom array and reindex connect array to match the new geom array"""
    geom, connect = read_geom_connect(sx)
    import pymesh

    nv = geom.shape[0]
    geom, connect, inf = pymesh.remove_duplicated_vertices_raw(geom, connect, tol=1e-4)
    print(f"removed {inf['num_vertex_merged']} duplicates as out {nv}")
    ind, ind_inv = lookup_sorted_geom(geom)
    geom = geom[ind, :]
    connect = np.array([ind_inv[x] for x in connect.flatten()]).reshape(connect.shape)
    # sort along line (then we can use multidim_intersect)
    connect = np.sort(connect, axis=1)
    return geom, connect


def multidim_intersect(arr1, arr2):
    """find indexes of same triangles in 2 connect arrays
    (associated with the same geom array)
    generate 1D arrays of tuples and use numpy function
    https://stackoverflow.com/questions/9269681/intersection-of-2d-numpy-ndarrays
    """
    arr1_view = arr1.view([("", arr1.dtype)] * arr1.shape[1])
    arr2_view = arr2.view([("", arr2.dtype)] * arr2.shape[1])
    intersected, ind1, ind2 = np.intersect1d(arr1_view, arr2_view, return_indices=True)
    ni, n1, n2 = intersected.shape[0], arr1.shape[0], arr2.shape[0]
    print(
        f"{ni} faces in common, n faces connect 1:{n1}, 2:{n2} (diff: {n1-ni}, {n2-ni})"
    )
    return ind1, ind2


def same_geometry(sx1, sx2):
    geom1 = sx1.ReadGeometry()
    geom2 = sx2.ReadGeometry()
    return np.all(np.isclose(geom1, geom2, rtol=1e-3, atol=1e-4))


sx1 = sx.seissolxdmf(args.xdmf_filename1)
sx2 = sx.seissolxdmf(args.xdmf_filename2)

same_geom = same_geometry(sx1, sx2)

if same_geom:
    print("same indexing detected, no need to reindex arrays")
    geom1, connect1 = read_geom_connect(sx1)
    geom2, connect2 = read_geom_connect(sx2)
else:
    geom1, connect1 = return_sorted_geom_connect(sx1)
    geom2, connect2 = return_sorted_geom_connect(sx2)
    if not np.all(np.isclose(geom1, geom2, rtol=1e-3, atol=1e-4)):
        raise ValueError("geometry arrays differ")
    ind1, ind2 = multidim_intersect(connect1, connect2)
    connect1 = connect1[ind1, :]

if args.idt[0] == -1:
    args.idt = list(range(0, sx1.ndt))

aData = []
for dataname in args.Data:
    print(dataname)
    myData1 = read_reshape2d(sx1, dataname)
    myData2 = read_reshape2d(sx2, dataname)
    if same_geom:
        myData = myData1[:, :] - myData2[:, :]
        if args.ratio:
            myData = myData / myData1[:, :]
    else:
        myData = myData1[:, ind1] - myData2[:, ind2]
        if args.ratio:
            myData = myData / myData1[:, ind1]

    for idt in args.idt:
        print(idt, np.amin(myData[idt, :]), np.amax(myData[idt, :]))
    aData.append(myData)
prefix, ext = os.path.splitext(args.xdmf_filename1)
add2prefix = "ratio" if args.ratio else "diff"
fname = f"{add2prefix}_{os.path.basename(prefix)}"

try:
    dt = sx1.ReadTimeStep()
except NameError:
    dt = 0.0

write_xdmfh5(fname, args.Data, geom1, connect1, aData, dt)
