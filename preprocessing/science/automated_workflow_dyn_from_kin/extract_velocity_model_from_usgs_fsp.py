#!/usr/bin/env python3


def read_velocity_model_from_fsp_file(fname):
    import re
    import pandas as pd
    from io import StringIO

    with open(fname, "r") as fid:
        lines = fid.readlines()

    if "FINITE-SOURCE RUPTURE MODEL" not in lines[0]:
        raise ValueError("Not a valid USGS fsp file.")

    def read_param(line, name, dtype=int):
        if name not in line:
            raise ValueError(f"{name} not found in line: {line}")
        else:
            return dtype(line.split(f"{name} =")[1].split()[0])

    def get_to_first_line_starting_with(lines, pattern):
        for i, line in enumerate(lines):
            if line.startswith(pattern):
                return lines[i:]
        raise ValueError(f"{pattern} not found")

    lines = get_to_first_line_starting_with(lines, "% VELOCITY-DENSITY")
    nlayers = read_param(lines[1], "layers")
    text_file = StringIO("\n".join(lines[3 : 5 + nlayers]))

    df = pd.read_csv(text_file, sep="\s+").drop([0])
    df = df.apply(pd.to_numeric, errors="coerce")
    rows_to_remove = df[df["DEPTH"] == 0]
    print("removing row\n", rows_to_remove)
    df = df[df["DEPTH"] > 0].reset_index(drop=True)
    if 'S-VEL' in df:
        # old usgs format
        df.rename(columns={'S-VEL': 'S_VEL'}, inplace=True)
        df.rename(columns={'P-VEL': 'P_VEL'}, inplace=True)

    df["rho"] = 1000.0 * df["DENS"]
    df["mu"] = 1e6 * df["rho"] * df["S_VEL"] ** 2
    df["lambda"] = 1e6 * df["rho"] * (df["P_VEL"] ** 2 - 2.0 * df["S_VEL"] ** 2)
    df.at[0, "DEPTH"] = -10
    print(df)
    return df


def write_yaml_material_file(df):
    to_write = """
!LayeredModel
map: !AffineMap
  matrix:
    z: [0.0, 0.0, 1.0]
  translation:
    z: 0
interpolation: upper
parameters: [rho, mu, lambda, Qp, Qs]
nodes:\n"""
    for index, row in df.iterrows():
        to_write += (
            f"   {-1e3*row['DEPTH']}:"
            f" [{row['rho']},{row['mu']:.10e},{row['lambda']:.10e}, {row['QP']},"
            f" {row ['QS']}]"
        )
        to_write += f" #[{row['P_VEL']}, {row['S_VEL']}]\n"

    fname = "yaml_files/usgs_material.yaml"
    with open(fname, "w") as fid:
        fid.write(to_write)
    print(f"done writing {fname}")


def write_z_rigidity_to_txt(df):
    """to be used in the teleseismic routines"""
    to_write = ""
    first_index = True
    G_prev = False
    eps = 1e-5
    for index, row in df.iterrows():
        if G_prev:
            to_write += f"{-1e3*row['DEPTH']-eps} {G_prev:.10e}\n"
        to_write += f"{-1e3*row['DEPTH']} {row['mu']:.10e}\n"
        G_prev = row["mu"]
    to_write += f"-1e10 {G_prev:.10e}\n"
    fname = "tmp/usgs_rigidity.txt"
    with open(fname, "w") as fid:
        fid.write(to_write)
    print(f"done writing {fname}")


if __name__ == "__main__":
    df = read_velocity_model_from_fsp_file("tmp/complete_inversion.fsp")
    write_yaml_material_file(df)
    write_z_rigidity_to_txt(df)
