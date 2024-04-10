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
import itertools
from compile_scenario_macro_properties import infer_duration
import warnings
import seissolxdmf as sx

if not os.path.exists("yaml_files"):
    os.makedirs("yaml_files")


def compute_max_slip(fn):
    sx0 = sx.seissolxdmf(fn)
    ndt = sx0.ReadNdt()
    ASl = sx0.ReadData("ASl", ndt - 1)
    if np.any(np.isnan(ASl)):
        ASl = sx0.ReadData("ASl", ndt - 2)
    return ASl.max()


usgs_fn = "output/dyn-usgs-fault.xdmf"
max_slip = compute_max_slip(usgs_fn)
assert max_slip > 0

# Get the directory of the script
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
input_file_dir = f"{script_directory}/input_files"
templateLoader = jinja2.FileSystemLoader(searchpath=input_file_dir)
templateEnv = jinja2.Environment(loader=templateLoader)
number_of_segments = len(glob.glob(f"tmp/*.ts"))
print(f"found {number_of_segments} segments")

# mode = "latin_hypercube"
# mode = "grid_search"
mode = "picked_models"
longer_and_more_frequent_output = True

if mode == "latin_hypercube":
    # R, B, C
    l_bounds = [0.55, 0.8, 0.1]
    u_bounds = [0.95, 1.2, 0.3]
    # cohesion is fixed in this appraoch
    list_cohesion = [(0.25, 1)]

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
    column_of_zeros = np.zeros((1, nsample))
    pars = np.insert(pars, 0, column_of_zeros, axis=1)

elif mode == "grid_search":
    # grid parameter space
    paramB = [0.9, 1.0, 1.1, 1.2]
    # paramB = [1.0]
    paramC = [0.1, 0.15, 0.2, 0.25, 0.3]
    # paramC = [0.3]
    paramR = [0.55, 0.6, 0.65, 0.7, 0.8, 0.9]
    # paramR = [0.65]
    list_cohesion = [(0.25, 1)]
    # list_cohesion = [(0.25, 0), (0.25, 1), (0.25, 3)]
    paramCoh = list(range(len(list_cohesion)))
    use_R_segment_wise = True
    if use_R_segment_wise:
        params = [paramCoh, paramB, paramC] + [paramR] * number_of_segments
        # params = [paramCoh, paramB, paramC] + [[0.85]] + [paramR] * (number_of_segments-1)
        assert len(params) == number_of_segments + 3
    else:
        params = [paramCoh, paramB, paramC, paramR]
        assert len(params) == 4
    # Generate all combinations of parameter values
    param_combinations = list(itertools.product(*params))
    # Convert combinations to numpy array and round to desired decimals
    pars = np.around(np.array(param_combinations), decimals=3)
    print(pars)
elif mode == "picked_models":
    list_cohesion = [(0.25, 0), (0.25, 1), (0.25, 2.5)]
    pars = [[0, 0.9, 0.3, 0.65], [1, 1.0, 0.3, 0.65], [2, 1.2, 0.3, 0.65]]
    pars = np.array(pars)
else:
    raise NotImplementedError(f"unkown mode {mode}")

nsample = pars.shape[0]
print(f"parameter space has {nsample} samples")


def render_file(template_par, template_fname, out_fname, verbose=True):
    template = templateEnv.get_template(template_fname)
    outputText = template.render(template_par)
    fn_tractions = out_fname
    with open(out_fname, "w") as fid:
        fid.write(outputText)
    if verbose:
        print(f"done creating {out_fname}")


hypo_z = np.loadtxt("tmp/hypocenter.txt")[2] * -1e3

mr_usgs = np.loadtxt("tmp/moment_rate.mr", skiprows=2)
usgs_duration = infer_duration(mr_usgs[:, 0], mr_usgs[:, 1])


def compute_fault_sampling(usgs_duration):
    if usgs_duration < 15:
        return 0.25
    elif usgs_duration < 30:
        return 0.5
    elif usgs_duration < 60:
        return 1.0
    elif usgs_duration < 200:
        return 2.5
    else:
        return 5.0


fault_sampling = compute_fault_sampling(usgs_duration)


def generate_R_yaml_block(Rvalues):
    if len(Rvalues) == 1:
        return f"""        [R]: !ConstantMap
            map:
              R: {Rvalues[0]}"""

    R_yaml_block = """        [R]: !Any
           components:"""
    for p, Rp in enumerate(Rvalues):
        fault_id = 3 if p == 0 else 64 + p
        R_yaml_block += f"""
            - !GroupFilter
              groups: {fault_id}
              components: !ConstantMap
                 map:
                   R: {Rp}"""
    return R_yaml_block


list_fault_yaml = []
for i in range(nsample):
    row = pars[i, :]
    cohi, B, C = row[0:3]
    cohesion_const, cohesion_lin = list_cohesion[int(cohi)]
    R = row[3:]

    template_par = {
        "R_yaml_block": generate_R_yaml_block(R),
        "cohesion_const": cohesion_const * 1e6,
        "cohesion_lin": cohesion_lin * 1e6,
        "B": B,
        "C": C,
        "min_dc": C * max_slip * 0.15,
        "hypo_z": hypo_z,
        "r_crit": 3000.0,
    }

    sR = "_".join(map(str, R))
    code = f"coh{cohesion_const}_{cohesion_lin}_B{B}_C{C}_R{sR}"
    fn_fault = f"yaml_files/fault_{code}.yaml"
    list_fault_yaml.append(fn_fault)

    render_file(template_par, "fault.tmpl.yaml", fn_fault)

    if longer_and_more_frequent_output:
        template_par["end_time"] = usgs_duration + max(20.0, 0.25 * usgs_duration)
        template_par["terminatorMomentRateThreshold"] = -1
        template_par["surface_output_interval"] = 1.0
    else:
        template_par["end_time"] = usgs_duration + max(5.0, 0.25 * usgs_duration)
        template_par["terminatorMomentRateThreshold"] = 5e17
        template_par["surface_output_interval"] = 5.0
    template_par["fault_fname"] = fn_fault
    template_par["output_file"] = f"output/dyn_{code}"
    template_par["material_fname"] = "yaml_files/usgs_material.yaml"
    template_par["fault_print_time_interval"] = fault_sampling
    fn_param = f"parameters_dyn_{code}.par"
    render_file(template_par, "parameters_dyn.tmpl.par", fn_param)


fnames = ["smooth_PREM_material.yaml", "mud.yaml", "fault_slip.yaml"]
for fn in fnames:
    shutil.copy(f"{input_file_dir}/{fn}", f"yaml_files/{fn}")

if os.path.exists("output/fl33-fault.xdmf"):
    fl33_file = "output/fl33-fault.xdmf"
elif os.path.exists("extracted_output/fl33_extracted-fault.xdmf"):
    fl33_file = "extracted_output/fl33_extracted-fault.xdmf"
else:
    raise FileNotFoundError(
        f"The files output/fl33-fault.xdmf or extracted_output/fl33_extracted-fault.xdmf were not found."
    )

list_nucleation_size = compute_critical_nucleation(
    fl33_file,
    "yaml_files/smooth_PREM_material.yaml",
    "yaml_files/fault_slip.yaml",
    list_fault_yaml,
    -hypo_z,
)
print(list_nucleation_size)
for i, fn in enumerate(list_fault_yaml):
    row = pars[i, :]
    cohi, B, C = row[0:3]
    cohesion_const, cohesion_lin = list_cohesion[int(cohi)]
    R = row[3:]
    sR = "_".join(map(str, R))
    code = f"coh{cohesion_const}_{cohesion_lin}_B{B}_C{C}_R{sR}"
    if list_nucleation_size[i]:
        fn_fault = f"yaml_files/fault_{code}.yaml"
        assert fn_fault == fn
        template_par = {
            "R_yaml_block": generate_R_yaml_block(R),
            "cohesion_const": cohesion_const * 1e6,
            "cohesion_lin": cohesion_lin * 1e6,
            "B": B,
            "C": C,
            "min_dc": C * max_slip * 0.15,
            "hypo_z": hypo_z,
            "r_crit": list_nucleation_size[i],
        }
        render_file(template_par, "fault.tmpl.yaml", fn_fault)
    else:
        fn_param = f"parameters_dyn_{code}.par"
        print(f"removing {fn} and {fn_param} (nucleation too large)")
        os.remove(fn)
        os.remove(fn_param)
