#!/usr/bin/env python3
import numpy as np
import glob
import re
import os
import jinja2
from sklearn.decomposition import PCA
import shutil
from estimate_nucleation_radius import compute_critical_nucleation
from scipy.stats import qmc
import random


if not os.path.exists("yaml_files"):
    os.makedirs("yaml_files")

# Get the directory of the script
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
input_file_dir = f"{script_directory}/input_files"
templateLoader = jinja2.FileSystemLoader(searchpath=input_file_dir)
templateEnv = jinja2.Environment(loader=templateLoader)

# R, B, C
l_bounds = [0.55, 0.8, 0.1]
u_bounds = [0.95, 1.2, 0.3]

if not os.path.exists("tmp/seed.txt"):
    seed = random.randint(1, 1000000)
    # keep the seed for reproducibility
    with open("tmp/seed.txt", "w+") as fid:
        fid.write(f"{seed}\n")
else:
    with open("tmp/seed.txt", "r") as fid:
        seed = int(fid.readline())
    print("seed read from tmp/seed.txt")


sampler = qmc.LatinHypercube(d=3, seed=seed)
nsample = 50
sample = sampler.random(n=nsample)
pars = qmc.scale(sample, l_bounds, u_bounds)
pars = np.around(pars, decimals=3)


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
for i in range(nsample):
    R, B, C = pars[i, :]
    template_par = {"R": R, "B": B, "C": C, "hypo_z": hypo_z, "r_crit": 3000.0}
    fn_fault = f"yaml_files/fault_B{B}_C{C}_R{R}.yaml"
    list_fault_yaml.append(fn_fault)

    render_file(template_par, "fault.tmpl.yaml", fn_fault)

    template_par["end_time"] = 60.0
    template_par["fault_fname"] = fn_fault
    template_par["output_file"] = f"output/dyn_B{B}_C{C}_R{R}"

    fn_param = f"parameters_dyn_B{B}_C{C}_R{R}.par"
    render_file(template_par, "parameters_dyn.tmpl.par", fn_param)


fnames = ["smooth_PREM_material.yaml", "mud.yaml", "rake.yaml", "fault_slip.yaml"]
for fn in fnames:
    shutil.copy(f"{input_file_dir}/{fn}", f"yaml_files/{fn}")

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
        R, B, C = pars[i, :]
        fn_fault = f"yaml_files/fault_B{B}_C{C}_R{R}.yaml"
        assert fn_fault == fn
        template_par = {
            "R": R,
            "B": B,
            "C": C,
            "hypo_z": hypo_z,
            "r_crit": list_nucleation_size[i],
        }
        render_file(template_par, "fault.tmpl.yaml", fn_fault)
    else:
        code = fn.split(".yaml")[0].split("fault_")[1]
        fn_param = f"parameters_dyn_{code}.par"
        print(f"removing {fn} and {fn_param} (nuclation too large)")
        os.remove(fn)
        os.remove(fn_param)
