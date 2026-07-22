#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import argparse
import json
import sys

import numpy as np
import seissolxdmf as sx

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two meshes.")
    parser.add_argument("mesh", type=str)
    parser.add_argument("mesh_ref", type=str, nargs="?", default=None)
    parser.add_argument(
        "--list-quantities",
        action="store_true",
        help="Print the output's data-field names as a JSON array and exit "
        "(mesh_ref is not needed).",
    )
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument(
        "--category",
        type=str,
        default="mesh",
        help="Label stored in the JSON summary (e.g. volume/fault/surface).",
    )
    parser.add_argument(
        "--geom-epsilon",
        type=float,
        default=1e-10,
        help="Tolerance for the geometry-equality check between output and "
        "reference. Geometry is currently always written in double, so the "
        "strict default is correct; the knob exists for future cases where "
        "output and reference geometry might differ.",
    )
    parser.add_argument(
        "--report-json",
        type=str,
        default=None,
        help="Write a machine-readable summary of the achieved errors to this "
        "path (always written, regardless of pass/fail). Does not affect the "
        "exit code.",
    )

    args = parser.parse_args()

    import meshcompare

    if args.list_quantities:
        print(json.dumps(meshcompare.list_quantities(args.mesh)))
        sys.exit(0)

    if args.mesh_ref is None:
        parser.error("mesh_ref is required unless --list-quantities is given")

    meshcompare.compare(
        args.mesh,
        args.mesh_ref,
        args.epsilon,
        report_json=args.report_json,
        category=args.category,
        geom_epsilon=args.geom_epsilon,
    )
