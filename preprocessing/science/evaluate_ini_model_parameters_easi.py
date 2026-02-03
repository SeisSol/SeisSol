#!/usr/bin/env python3
import argparse

import easi
import numpy as np
import pyvista as pv
import seissolxdmf
import seissolxdmfwriter as sxw


class SeissolxdmfExtended(seissolxdmf.seissolxdmf):
    def __init__(self, xdmfFilename, subdivide_level=0):
        super().__init__(xdmfFilename)

        self.xyz = self.ReadGeometry()
        self.connect = self.ReadConnect().astype(np.int64)
        self.is_volume = self.connect.shape[1] == 4
        self.subdivide_level = subdivide_level

        if subdivide_level > 0:
            n_elements = self.connect.shape[0]
            vertices_per_cell = 4 if self.is_volume else 3

            # Prepare cells array with VTK cell format
            cells = np.hstack(
                [
                    np.full((n_elements, 1), vertices_per_cell, dtype=np.int64),
                    self.connect,
                ]
            ).ravel()

            # Create mesh and subdivide
            if self.is_volume:
                celltypes = np.ascontiguousarray(
                    np.full(n_elements, pv.CellType.TETRA, dtype=np.uint8)
                )
                mesh = pv.UnstructuredGrid(cells, celltypes, self.xyz)
                for _ in range(subdivide_level):
                    mesh = mesh.subdivide_tetra()
            else:
                mesh = pv.PolyData(self.xyz, cells)
                mesh = mesh.subdivide(subdivide_level, "linear")

            # Update geometry and connectivity
            self.xyz = mesh.points
            n = vertices_per_cell + 1
            self.connect = mesh.faces.reshape(-1, n)[:, 1:n].astype(np.int64)

    def ReadFaultTagOrRegion(self):
        """Read fault-tag or region array, optionally subdivided multiple levels."""
        var_name = "group" if self.is_volume else "fault-tag"

        # Read original tags as a 1D integer array
        original_tags = self.Read1dData(var_name, self.nElements, isInt=True)

        # If subdividing, repeat tags according to total number of children
        if self.subdivide_level > 0:
            children_per_parent = 12 if self.is_volume else 4
            total_repeat = children_per_parent**self.subdivide_level
            tags = np.repeat(original_tags, total_repeat)
        else:
            tags = original_tags

        return tags

    def ComputeCellCenters(self):
        """compute cell center array"""
        return self.xyz[self.connect].mean(axis=1)

    def ComputeCellNormals(self, ref_vector):
        v0 = self.xyz[self.connect[:, 0]]
        v1 = self.xyz[self.connect[:, 1]]
        v2 = self.xyz[self.connect[:, 2]]

        un = np.cross(v1 - v0, v2 - v0)
        un /= np.linalg.norm(un, axis=1, keepdims=True)

        # Orient normals
        mysign = np.sign(un @ ref_vector)
        un *= mysign[:, None]
        return un


def compute_tractions(dicStress, un):
    # compute Strike and dip direction
    us = np.stack([-un[:, 1], un[:, 0], np.zeros_like(un[:, 0])], axis=1)
    us /= np.linalg.norm(us, axis=1, keepdims=True)
    ud = np.cross(un, us)

    def compute_traction_vector(dicStress, un):
        Tractions = np.zeros_like(un)
        a = np.stack((dicStress["s_xx"], dicStress["s_xy"], dicStress["s_xz"]), axis=1)
        Tractions[:, 0] = np.sum(a * un, axis=1)
        a = np.stack((dicStress["s_xy"], dicStress["s_yy"], dicStress["s_yz"]), axis=1)
        Tractions[:, 1] = np.sum(a * un, axis=1)
        a = np.stack((dicStress["s_xz"], dicStress["s_yz"], dicStress["s_zz"]), axis=1)
        Tractions[:, 2] = np.sum(a * un, axis=1)
        return Tractions

    Tractions = compute_traction_vector(dicStress, un)

    outDict = {}
    # compute Traction components
    outDict["T_n"] = np.sum(Tractions * un, axis=1)
    outDict["T_s"] = np.sum(Tractions * us, axis=1)
    outDict["T_d"] = np.sum(Tractions * ud, axis=1)
    return outDict


if __name__ == "__main__":
    # parsing python arguments
    parser = argparse.ArgumentParser(
        description=(
            " retrieve/compute model parameter(s) from easi/yaml file and seissol file"
        )
    )
    parser.add_argument("seissol_filename", help="seissol xdmf file (mesh or output)")
    parser.add_argument("yaml_filename", help="easi/yaml filename")
    parser.add_argument(
        "--parameters",
        help=(
            "variable to be read in the yaml file. 'tractions' for initial stress. "
            "Vp and Vs for P and S wave velocities (coma separated string)."
        ),
        required=True,
    )
    parser.add_argument(
        "--ref_vector",
        nargs=1,
        help="reference vector (see seissol parameter file) used to choose fault normal (coma separated string)",
    )
    parser.add_argument(
        "--subdivide_level",
        default=0,
        type=int,
        help="subdivide mesh before evaluating",
    )

    parser.add_argument(
        "--output_prefix",
        help="path and prefix of the output file",
    )
    args = parser.parse_args()
    if not args.output_prefix:
        args.output_prefix = "_".join(args.parameters.split(","))

    sx = SeissolxdmfExtended(args.seissol_filename, args.subdivide_level)
    centers = sx.ComputeCellCenters()
    tags = sx.ReadFaultTagOrRegion()

    parameters = [p.strip() for p in args.parameters.split(",")]
    available = sx.ReadAvailableDataFields()

    param_evaluated = [p for p in parameters if p not in ["tractions", "Vp", "Vs"]]
    out = {}

    # Evaluate directly available fields
    if param_evaluated:
        out.update(
            easi.evaluate_model(centers, tags, param_evaluated, args.yaml_filename)
        )
        parameters = [x for x in parameters if x not in param_evaluated]

    if "Vp" in parameters or "Vs" in parameters:
        base = easi.evaluate_model(
            centers, tags, ["rho", "mu", "lambda"], args.yaml_filename
        )

        if "Vp" in parameters:
            out["Vp"] = np.sqrt((base["lambda"] + 2.0 * base["mu"]) / base["rho"])

        if "Vs" in parameters:
            out["Vs"] = np.sqrt(base["mu"] / base["rho"])

    if "tractions" in parameters:
        if sx.is_volume:
            raise ValueError(
                "the xdmf given file is a volume and 'tractions' given in parameters"
            )

        direct_tractions = ["T_s", "T_d", "T_n"]
        stress_components = ["s_xx", "s_yy", "s_zz", "s_xy", "s_xz", "s_yz"]
        try:
            base = easi.evaluate_model(
                centers, tags, direct_tractions, args.yaml_filename
            )
            out.update(base)
        except ValueError:
            print(f"[T_s, T_d, T_n] not found in {args.yaml_filename}, using s_ij")
            if not args.ref_vector:
                raise ValueError(
                    "ref_vector has to be defined for computing tractions from stress"
                )
            ref_vector = [float(v) for v in args.ref_vector[0].split(",")]
            normals = sx.ComputeCellNormals(ref_vector)
            dicStress = easi.evaluate_model(
                centers,
                tags,
                stress_components,
                args.yaml_filename,
            )
            out_tractions = compute_tractions(dicStress, normals)
            out.update(out_tractions)

    sxw.write(
        args.output_prefix,
        sx.xyz,
        sx.connect,
        out,
        {"0": 0},
        reduce_precision=True,
    )
