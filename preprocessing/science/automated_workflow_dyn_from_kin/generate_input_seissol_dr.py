#!/usr/bin/env python3
import numpy as np
import glob
import re
import os
import jinja2
from sklearn.decomposition import PCA
import shutil

if not os.path.exists("yaml_files"):
    os.makedirs("yaml_files")

# Get the directory of the script
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
input_file_dir = f"{script_directory}/input_files"
templateLoader = jinja2.FileSystemLoader(searchpath=input_file_dir)
templateEnv = jinja2.Environment(loader=templateLoader)


aR = [0.7, 0.9, 0.95]
aB = [0.9, 1.0, 1.1]
aC = [0.1, 0.2, 0.4]


def render_file(template_par, template_fname, out_fname):
    template = templateEnv.get_template(template_fname)
    outputText = template.render(template_par)
    fn_tractions = out_fname
    with open(out_fname, "w") as fid:
        fid.write(outputText)
    print(f"done creating {out_fname}")


for B in aB:
    template_par = {"B": B}
    fn_tractions = f"yaml_files/tractions_B{B}.yaml"
    render_file(template_par, "tractions.tmpl.yaml", fn_tractions)

hypo_z = np.loadtxt("tmp/hypocenter.txt")[2] * -1e3

for C in aC:
    template_par = {"C": C, "hypo_z": hypo_z}
    fn_common = f"yaml_files/common2all_C{C}.yaml"
    render_file(template_par, "common2all.tmpl.yaml", fn_common)

for B in aB:
    fn_tractions = f"yaml_files/tractions_B{B}.yaml"
    for C in aC:
        fn_common = f"yaml_files/common2all_C{C}.yaml"
        for R in aR:
            template_par = {
                "R": R,
                "common2all_fname": fn_common,
                "tractions_fname": fn_tractions,
            }
            fn_fault = f"yaml_files/fault_B{B}_C{C}_R{R}.yaml"
            render_file(template_par, "fault.tmpl.yaml", fn_fault)

            template_par["end_time"] = 60.0
            template_par["fault_fname"] = fn_fault
            template_par["output_file"] = f"output/dyn_B{B}_C{C}_R{R}"

            fn_param = f"parameters_dyn_B{B}_C{C}_R{R}.par"
            render_file(template_par, "parameters_dyn.tmpl.par", fn_param)

shutil.copy(
    f"{input_file_dir}/smooth_PREM_material.yaml",
    "yaml_files/smooth_PREM_material.yaml",
)
shutil.copy(f"{input_file_dir}/mud.yaml", "yaml_files/mud.yaml")
shutil.copy(f"{input_file_dir}/rake.yaml", "yaml_files/rake.yaml")
