#!/usr/bin/env python3
from extract_velocity_model_from_usgs_fsp import (
    write_yaml_material_file,
    write_z_rigidity_to_txt,
)

if __name__ == "__main__":
    import pandas as pd

    df = pd.read_csv("tmp/vel_model_slipnear.txt", sep=" ")
    df["rho"] = 1000.0 * df["DENS"]
    df["mu"] = 1e6 * df["rho"] * df["S_VEL"] ** 2
    df["lambda"] = 1e6 * df["rho"] * (df["P_VEL"] ** 2 - 2.0 * df["S_VEL"] ** 2)
    df["DEPTH_bot"] = df["H"].cumsum()
    df["DEPTH"] = df["DEPTH_bot"].shift(1)
    df.at[0, "DEPTH"] = -10

    print(df)

    write_yaml_material_file(df)
    write_z_rigidity_to_txt(df)
