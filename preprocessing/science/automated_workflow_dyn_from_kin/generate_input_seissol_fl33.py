#!/usr/bin/env python3
import numpy as np
import glob
import re
import os
import jinja2
import shutil
from scipy.spatial.distance import pdist

# Get the directory of the script
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
input_file_dir = f"{script_directory}/input_files"
templateLoader = jinja2.FileSystemLoader(searchpath=input_file_dir)
templateEnv = jinja2.Environment(loader=templateLoader)

# Regular expression patterns to match vertex lines
vertex_pattern = re.compile(r"VRTX (\d+) ([\d.-]+) ([\d.-]+) ([\d.-]+)")
allv = []
faults = []
ts_files = sorted(glob.glob(f"tmp/*.ts"))

for i, fn in enumerate(ts_files):
    vertices = []
    with open(fn, "r") as file:
        for line in file:
            match = vertex_pattern.match(line)
            if match:
                vertex_id, x, y, z = map(float, match.groups())
                vertices.append([x, y, z])
    allv.extend(vertices)
xyz = np.array(allv)

# Find the maximum distance
max_distance = np.max(pdist(xyz))

# 3200 is the smallest wave speed in the model
end_time = max_distance / 3200.0
template_par = {}
# well in theory we would need to run for end_time, but practically
# a portion of it may be sufficient
template_par["end_time"] = max(25.0, 0.6 * end_time)
template_par["material_fname"] = "yaml_files/usgs_material.yaml"

template = templateEnv.get_template(f"parameters_fl34.tmpl.par")
outputText = template.render(template_par)
fname = "parameters_fl34.par"
with open(fname, "w") as fid:
    fid.write(outputText)
print(f"done creating {fname}")
shutil.copy(
    f"{input_file_dir}/smooth_PREM_material.yaml",
    "yaml_files/smooth_PREM_material.yaml",
)
