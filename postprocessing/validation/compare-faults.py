#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2022 SeisSol Group
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
    parser = argparse.ArgumentParser(description="Compare two faults.")
    parser.add_argument("fault", type=str)
    parser.add_argument("fault_ref", type=str)
    parser.add_argument("--epsilon", type=float, default=0.01)

    args = parser.parse_args()

    import meshcompare

    meshcompare.compare(args.fault, args.fault_ref, args.epsilon)
