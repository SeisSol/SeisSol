#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import argparse
import sys

import numpy as np
import seissolxdmf as sx

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two meshes.")
    parser.add_argument("mesh", type=str)
    parser.add_argument("mesh_ref", type=str)
    parser.add_argument("--epsilon", type=float, default=0.01)

    args = parser.parse_args()

    import meshcompare

    meshcompare.compare(args.mesh, args.mesh_ref, args.epsilon)
