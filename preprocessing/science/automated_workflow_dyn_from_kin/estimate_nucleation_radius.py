#!/usr/bin/env python3
import numpy as np
import argparse
from scipy import interpolate
import easi
import seissolxdmf
from tqdm import tqdm


class seissolxdmfExtended(seissolxdmf.seissolxdmf):
    def __init__(self, xdmfFilename):
        super().__init__(xdmfFilename)
        self.xyz = self.ReadGeometry()
        self.connect = self.ReadConnect()

    def ReadFaultTags(self):
        """Read partition array"""
        return self.Read1dData("fault-tag", self.nElements, isInt=True).T

    def ComputeCellCenters(self):
        """compute cell center array"""
        return (
            self.xyz[self.connect[:, 0]]
            + self.xyz[self.connect[:, 1]]
            + self.xyz[self.connect[:, 2]]
        ) / 3.0

    def ComputeCellNormals(self):
        """compute cell normal"""
        cross = np.cross(
            self.xyz[self.connect[:, 1], :] - self.xyz[self.connect[:, 0], :],
            self.xyz[self.connect[:, 2], :] - self.xyz[self.connect[:, 0], :],
        )
        norm = np.apply_along_axis(np.linalg.norm, 1, cross)
        return cross / norm.reshape((self.nElements, 1))

    def ComputeCellAreas(self):
        """compute area of each cell"""
        cross = np.cross(
            self.xyz[self.connect[:, 1], :] - self.xyz[self.connect[:, 0], :],
            self.xyz[self.connect[:, 2], :] - self.xyz[self.connect[:, 0], :],
        )
        return 0.5 * np.apply_along_axis(np.linalg.norm, 1, cross)


def k_int(x):
    "K integral in Uenishi 2009, see Galis 2015, eq 20"
    ti = np.linspace(0, 1, 250, endpoint=False)
    f = 1.0 / np.sqrt((1.0 - ti**2) * (1.0 - x**2 * ti**2))  # for t in ti]
    return np.trapz(f, ti)


def e_int(x):
    "E integral in Uenishi 2009, see Galis 2015, eq 21"
    ti = np.linspace(0, 1, 250, endpoint=False)
    f = [np.sqrt((1.0 - x**2 * t**2) / (1.0 - t**2)) for t in ti]
    return np.trapz(f, ti)


def C_Uenishi(x):
    "C(nu) in Uenishi 2009, see Galis 2015, eq 19"
    a = np.sqrt(x * (2.0 - x))
    return (e_int(a) + (1.0 - x) * k_int(a)) / (2.0 - x)


def C_Uenishi_interp():
    "interpolated function based on C_Uenishi (much faster)"
    x = np.arange(0.14, 0.35, 0.001)
    y = np.array([C_Uenishi(xi) for xi in x])
    return interpolate.interp1d(x, y)


def stiffness_average(G1, G2):
    return 2.0 * G1 * G2 / (G1 + G2)


def compute_CG(centers, ids, sx, mat_yaml):
    "compute nu across the fault, and compute C_Uenishi(nu)"
    fault_normal = sx.ComputeCellNormals()[ids]
    nx = fault_normal.shape[0]
    regions = np.ones((nx, 1))
    centers_mp = np.vstack((centers + 0.1 * fault_normal, centers - 0.1 * fault_normal))
    out = easi.evaluate_model(centers_mp, regions, ["lambda", "mu"], mat_yaml)
    lambda_x = stiffness_average(out["lambda"][0:nx], out["lambda"][nx:])
    G = stiffness_average(out["mu"][0:nx], out["mu"][nx:])
    nu = 0.5 * lambda_x / (lambda_x + G)
    f = C_Uenishi_interp()
    return np.array([f(nui) for nui in nu]) * G


def points_in_sphere(points, center, radius):
    distances = np.linalg.norm(points - center, axis=1)
    within_sphere = distances <= radius
    return within_sphere


def compute_slip_area(centers, sx, slip_yaml):
    face_area = sx.ComputeCellAreas()
    tags = sx.ReadFaultTags()
    slip = easi.evaluate_model(centers, tags, ["fault_slip"], slip_yaml)["fault_slip"]
    return np.sum(face_area[slip > 0.05 * np.amax(slip)])


def compute_critical_nucleation(
    fault_xdmf, mat_yaml, slip_yaml, list_fault_yaml, hypo_depth
):
    sx = seissolxdmfExtended(fault_xdmf)
    centers = sx.ComputeCellCenters()
    slip_area = compute_slip_area(centers, sx, slip_yaml)

    center = np.array([0, 0, -hypo_depth])
    ids = points_in_sphere(centers, center, 15e3)
    centers = centers[ids]

    CG = compute_CG(centers, ids, sx, mat_yaml)
    tags = sx.ReadFaultTags()[ids]
    results = []
    for fault_yaml in tqdm(list_fault_yaml):
        print(f"now processing {fault_yaml}")
        out = easi.evaluate_model(
            centers, tags, ["mu_s", "mu_d", "d_c", "T_n"], fault_yaml
        )

        mu_s, mu_d, Dc, normalStress = out["mu_s"], out["mu_d"], out["d_c"], out["T_n"]
        W = -np.abs(mu_s - mu_d) * normalStress / Dc
        # L = 0.624 * CG / np.median(W)
        L = 0.624 * CG / W
        area_crit = np.pi * L**2

        radius = np.arange(0.5e3, 15.0e3, 0.25e3)
        nucRadius = False
        for rad in radius:
            ids = points_in_sphere(centers, center, rad)
            estimatedR = np.median(L[ids])
            if rad > 2.0 * estimatedR:
                ratio_slip_area = 100 * np.pi * rad**2 / slip_area
                print(
                    (
                        "selected_radius, estimated_radius, std, nuc_area/slip_area"
                        " (%) :"
                        f" {rad:.0f} {estimatedR:.0f} {np.std(L[ids]):.0f}"
                        f" {ratio_slip_area:.1f}"
                    ),
                )
                if ratio_slip_area > 15.0:
                    nucRadius = False
                else:
                    nucRadius = rad
                break
        if not nucRadius:
            rad = min(rad, np.sqrt((0.15 / np.pi) * slip_area))
            ratio_slip_area = 100 * np.pi * rad**2 / slip_area
            print(
                (
                    "selected_radius (forced to 15% of slip area),"
                    " estimated_radius, std, nuc_area/slip_area (%) :"
                    f" {rad:.0f} {estimatedR:.0f} {np.std(L[ids]):.0f}"
                    f" {ratio_slip_area:.1f}"
                ),
            )
            nucRadius = rad
        results.append(nucRadius)
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compute critical depletion")
    parser.add_argument(
        "faultfile", help="fault file (vtk or seissol xdmf fault output)"
    )
    parser.add_argument("faultyamlfile", help="yaml file describing fault parameters")
    parser.add_argument("hypocenter_depth", type=float, help="hypocenter depth")
    parser.add_argument(
        "--materialyamlfile",
        help="yaml file desribing Lame parameters",
        default="yaml_files/smooth_PREM_material.yaml",
    )
    parser.add_argument(
        "--slipyamlfile",
        help="yaml file describing fault slip",
        default="yaml_files/fault_slip.yaml",
    )
    args = parser.parse_args()
    fault_xdmf = args.faultfile
    fault_yaml = args.faultyamlfile
    slip_yaml = args.slipyamlfile
    mat_yaml = args.materialyamlfile
    hypo_depth = args.hypocenter_depth

    compute_critical_nucleation(
        fault_xdmf, mat_yaml, slip_yaml, [fault_yaml], hypo_depth
    )
