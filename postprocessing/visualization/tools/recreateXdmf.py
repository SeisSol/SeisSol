import argparse
import subprocess
import xml.etree.ElementTree
import os


def next_power_of_2(x):
    return 1 if x == 0 else 2 ** ((x - 1).bit_length())


def generate_new_prefix(prefix, append2prefix):
    if append2prefix == "":
        append2prefix = "_resampled"
    prefix = os.path.basename(prefix)
    lsplit = prefix.split("-")
    if len(lsplit) > 1:
        if lsplit[-1] in ["surface", "low", "fault"]:
            prefix0 = "-".join(lsplit[0:-1])
            prefix_new = prefix0 + append2prefix + "-" + lsplit[-1]
    else:
        prefix_new = prefix + append2prefix
    return prefix_new


def recreateXdmf(prefix, prefix_new, nvertex, ncells, nmem, dt, indices, lsData, tohdf5=False, prec=8, append2prefix="_resampled"):
    full_prefix = prefix
    prefix = os.path.basename(prefix)
    prefix_new = os.path.basename(prefix_new)

    # fault and surface output have 3 colums in connect
    ncolConnect = 4
    scell = "Tetrahedron"
    lsplit = prefix.split("-")
    if len(lsplit) > 1:
        if lsplit[-1] in ["surface", "fault"]:
            ncolConnect = 3
            scell = "Triangle"

    if tohdf5:
        colonOrNothing = ".h5:"
        DataExtension = ""
        DataFormat = "HDF"
    else:
        colonOrNothing = ""
        DataExtension = ".bin"
        DataFormat = "Binary"

    # create and print the new Xdmf file
    header = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
 <Domain>
  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">"""

    score = ""

    for i, ii in enumerate(indices):
        score = (
            score
            + """   <Grid Name="step_%012d" GridType="Uniform">
    <Topology TopologyType="%s" NumberOfElements="%d">
     <DataItem NumberType="Int" Precision="8" Format="%s" Dimensions="%d %d">%s_cell%s/mesh0/connect%s</DataItem>
    </Topology>
    <Geometry name="geo" GeometryType="XYZ" NumberOfElements="%d">
     <DataItem NumberType="Float" Precision="8" Format="%s" Dimensions="%d 3">%s_vertex%s/mesh0/geometry%s</DataItem>
    </Geometry>
    <Time Value="%f"/>"""
            % (
                ii,
                scell,
                ncells,
                DataFormat,
                ncells,
                ncolConnect,
                prefix_new,
                colonOrNothing,
                DataExtension,
                nvertex,
                DataFormat,
                nvertex,
                prefix_new,
                colonOrNothing,
                DataExtension,
                dt * ii,
            )
        )

        for sdata in lsData:
            if sdata == "partition":
                score = (
                    score
                    + """
    <Attribute Name="partition" Center="Cell">
     <DataItem  NumberType="Int" Precision="4" Format="%s" Dimensions="%d">%s_cell%s/mesh0/partition%s</DataItem>
    </Attribute>\n"""
                    % (DataFormat, ncells, prefix_new, colonOrNothing, DataExtension)
                )
            else:
                score = (
                    score
                    + """    <Attribute Name="%s" Center="Cell">
     <DataItem ItemType="HyperSlab" Dimensions="%d">
      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">%d 0 1 1 1 %d</DataItem>
      <DataItem NumberType="Float" Precision="%d" Format="%s" Dimensions="%d %d">%s_cell%s/mesh0/%s%s</DataItem>
     </DataItem>
    </Attribute>\n"""
                    % (sdata, ncells, i, ncells, prec, DataFormat, i + 1, nmem, prefix_new, colonOrNothing, sdata, DataExtension)
                )
        score = score + "   </Grid>\n"

    score = (
        score
        + """  </Grid>
 </Domain>
</Xdmf>"""
    )
    prefix0 = generate_new_prefix(full_prefix, append2prefix)
    fout = open(prefix0 + ".xdmf", "w")
    fout.write(header)
    fout.write("\n")
    fout.write(score)
    fout.write("\n")
    fout.close()
    print("done writing " + prefix0 + ".xdmf")


def ReadNcellsFromXdmf(xdmfFile):
    out = subprocess.check_output(["grep connect " + xdmfFile + " | head -n1"], shell=True)
    e = xml.etree.ElementTree.fromstring(out)
    dimstring = e.attrib["Dimensions"].split()
    return int(dimstring[0])


def ReadDtFromXdmf(xdmfFile):
    out = subprocess.check_output(["grep Time " + xdmfFile + " | head -n3| tail -n1"], shell=True)
    e = xml.etree.ElementTree.fromstring(out)
    dimstring = e.attrib["Value"].split()
    return float(dimstring[0])


def ReadNvertexFromXdmf(xdmfFile):
    out = subprocess.check_output(["grep geometry " + xdmfFile + " | head -n1"], shell=True)
    e = xml.etree.ElementTree.fromstring(out)
    dimstring = e.attrib["Dimensions"].split()
    return int(dimstring[0])


def ReadNdtNmemFromXdmf(xdmfFile):
    out = subprocess.check_output(["grep DataItem " + xdmfFile + " | tail -n2 | head -n1"], shell=True)
    e = xml.etree.ElementTree.fromstring(out)
    dimstring = e.attrib["Dimensions"].split()
    # return (int(dimstring[0])-1, int(dimstring[1]))
    return (int(dimstring[0]), int(dimstring[1]))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="recreate a xdmf file")
    parser.add_argument("prefix", help="prefix including -fault or -surface")
    parser.add_argument("--idt", nargs="+", help="list of time step to differenciate (1st = 0); -1 = all", type=int)
    parser.add_argument("--nvertex", nargs=1, metavar=("nvertex"), help="number of vertex (read if not given)", type=int)
    parser.add_argument("--ncells", nargs=1, metavar=("ncells"), help="number of cells (read if not given)", type=int)
    parser.add_argument("--dt", nargs=1, metavar=("dt"), help="output time step (read if not given)", type=float)
    parser.add_argument("--ndt", nargs=1, metavar=("ndt"), help="number of time steps to output (read if not given)", type=int)
    parser.add_argument("--Data", nargs="+", help="list of data variable to write")
    args = parser.parse_args()

    prefix = args.prefix
    xdmfFile = prefix + ".xdmf"
    ### Read all parameters from the xdmf file of SeisSol (if any)
    if args.ncells != None:
        ncells = args.ncells[0]
        nmem = next_power_of_2(ncells)
        # nmem = next_power_of_2(ncells)*2
    else:
        ncells = ReadNcellsFromXdmf(xdmfFile)

    if args.dt != None:
        dt = args.dt[0]
    else:
        dt = ReadDtFromXdmf(xdmfFile)

    if args.nvertex != None:
        nvertex = args.nvertex[0]
    else:
        nvertex = ReadNvertexFromXdmf(xdmfFile)

    if args.ndt != None:
        ndt = args.ndt[0]
    else:
        ndt, nmem = ReadNdtNmemFromXdmf(xdmfFile)

    recreateXdmf(prefix, prefix, nvertex, ncells, nmem, dt, range(ndt), args.Data)
