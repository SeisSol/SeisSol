#!/usr/bin/env python3
import seissolxdmf
import argparse
import numpy as np
import os


class seissolxdmfExtended(seissolxdmf.seissolxdmf):
    def ReadFaultTags(self):
        """Read partition array"""
        return self.Read1dData("fault-tag", self.nElements, isInt=True).T


# parsing python arguments
parser = argparse.ArgumentParser(
    description="compute average rake segment-wise (based on fault tagging)"
)
parser.add_argument("fault", help="fault.xdmf filename")

args = parser.parse_args()

# Compute fault centroids
sx = seissolxdmfExtended(args.fault)
tags = sx.ReadFaultTags()
unique_tags = np.unique(tags)
xyz = sx.ReadGeometry()
connect = sx.ReadConnect()

cross0 = np.cross(
    xyz[connect[:, 1], :] - xyz[connect[:, 0], :],
    xyz[connect[:, 2], :] - xyz[connect[:, 0], :],
)
face_area = 0.5 * np.apply_along_axis(np.linalg.norm, 1, cross0)

ndt = sx.ReadNdt()

strike_slip = sx.ReadData("Sls", ndt - 1)
dip_slip = sx.ReadData("Sld", ndt - 1)
slip = sx.ReadData("ASl", ndt - 1)
rake = np.arctan2(dip_slip, strike_slip) * 180 / np.pi


template_yaml = """!Any
components:
"""

for tag in unique_tags:
    ids = np.where(tags == tag)[0]
    face_area_sel = face_area[ids] * slip[ids]
    ruptured_area = np.sum(face_area_sel)
    average_rake = np.sum(rake[ids] * face_area_sel) / ruptured_area
    print(average_rake)

    template_yaml += f""" - !GroupFilter
    groups: {tag}
    components: !ConstantMap
              map:
                rake_lowslip: {average_rake}
"""

fname = "yaml_files/rake_lowslip.yaml"
with open(fname, "w") as fid:
    fid.write(template_yaml)
print(f"done writing {fname}")
