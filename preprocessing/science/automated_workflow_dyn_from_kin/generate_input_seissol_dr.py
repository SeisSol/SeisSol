#!/usr/bin/env python3
import numpy as np
import glob
import re
import os
import jinja2
from sklearn.decomposition import PCA
import shutil
from estimate_nucleation_radius import compute_critical_nucleation

if not os.path.exists("yaml_files"):
    os.makedirs("yaml_files")

# Get the directory of the script
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
input_file_dir = f"{script_directory}/input_files"
templateLoader = jinja2.FileSystemLoader(searchpath=input_file_dir)
templateEnv = jinja2.Environment(loader=templateLoader)


aR = [0.7, 0.85, 0.95]
aB = [0.9, 1.0, 1.1]
aC = [0.1, 0.2, 0.3]


def render_file(template_par, template_fname, out_fname, verbose=True):
    template = templateEnv.get_template(template_fname)
    outputText = template.render(template_par)
    fn_tractions = out_fname
    with open(out_fname, "w") as fid:
        fid.write(outputText)
    if verbose:
        print(f"done creating {out_fname}")


hypo_z = np.loadtxt("tmp/hypocenter.txt")[2] * -1e3

list_fault_yaml = []
for B in aB:
    for C in aC:
        for R in aR:
            template_par = {"R": R, "B": B, "C": C}
            fn_fault = f"yaml_files/fault_B{B}_C{C}_R{R}.yaml"
            list_fault_yaml.append(fn_fault)

            render_file(template_par, "fault.tmpl.yaml", fn_fault)

            template_par["end_time"] = 60.0
            template_par["fault_fname"] = fn_fault
            template_par["output_file"] = f"output/dyn_B{B}_C{C}_R{R}"

            fn_param = f"parameters_dyn_B{B}_C{C}_R{R}.par"
            render_file(template_par, "parameters_dyn.tmpl.par", fn_param)

list_nucleation_size = compute_critical_nucleation(
    "output/fl33-fault.xdmf",
    "yaml_files/smooth_PREM_material.yaml",
    "yaml_files/fault_slip.yaml",
    list_fault_yaml,
    -hypo_z,
)
print(list_nucleation_size)

for i, fn in enumerate(list_fault_yaml):
    if list_nucleation_size[i]:
        template_par = {"hypo_z": hypo_z, "r_crit": list_nucleation_size[i]}
        render_file(
            template_par, "forced_rupture.tmpl.yaml", "forced_rupture.yaml", False
        )

        # Append the content of forced_rupture.yaml into fn
        with open("forced_rupture.yaml", "r") as file2:
            content_file2 = file2.read()
        with open(fn, "a") as file1:
            file1.write(content_file2)
    else:
        code = fn.split(".yaml")[0].split("fault_")[1]
        fn_param = f"parameters_dyn_{code}.par"
        print(f"removing {fn} and {fn_param} (nuclation too large)")
        os.remove(fn)
        os.remove(fn_param)


fnames = ["smooth_PREM_material.yaml", "mud.yaml", "rake.yaml", "fault_slip.yaml"]
for fn in fnames:
    shutil.copy(f"{input_file_dir}/{fn}", f"yaml_files/{fn}")
