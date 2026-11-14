#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Write the point sources of the verification cases of the copies of the cells.

``source-s<N>.fsrm`` holds one source for every simulation of a run of N
simulations, each with a position, an orientation, a size, an onset and a
moment rate of its own, so that no two simulations carry the same wave.
``source-copies-s<N>.fsrm`` holds the same sources twice over: a run of 2N
simulations deals them out to its simulations in turn, so that its simulations
s and N + s carry the wave of the simulation s of the run of N simulations.

The sources push as well as they turn, so that the density of the material
around them scales them. Every position lies inside a cell of the mesh box48,
away from its faces.

    python3 sources.py cases 1 8 32
"""

import argparse
import math
import random
from pathlib import Path

# the step width of the cases, which the moment rates are sampled with
TIMESTEP = 6.103515625e-05
SAMPLES = 16

# box48: the cube [-50, 50]^3, cut into eight cubes and each of them into six
# tetrahedra along its diagonal, whose faces lie where two coordinates are
# equal relative to the cube
EDGE = 50.0
MARGIN = 0.05


def inside(point):
    """Whether a point lies inside a cell of box48, at least MARGIN cube edges from its faces."""
    local = [(value + 2 * EDGE) % EDGE / EDGE for value in point]
    if any(value < MARGIN or value > 1 - MARGIN for value in local):
        return False
    return all(
        abs(local[i] - local[j]) >= MARGIN for i in range(3) for j in range(i + 1, 3)
    )


def source(generator, index):
    """The parameters and the sampled moment rate of a source."""
    while True:
        point = [round(generator.uniform(-48.0, 48.0), 3) for _ in range(3)]
        if inside(point):
            break
    strike = round(generator.uniform(0.0, 360.0), 1)
    dip = round(generator.uniform(10.0, 90.0), 1)
    rake = round(generator.uniform(-180.0, 180.0), 1)
    area = round(generator.uniform(0.5, 2.0), 3)
    onset = round(generator.uniform(0.0, 4.0) * TIMESTEP, 9)
    centre = 4.0 + index % 3
    width = 1.0 + 0.25 * (index % 4)
    rate = [
        1.0e9 * math.exp(-(((sample - centre) / width) ** 2))
        for sample in range(SAMPLES)
    ]
    return point + [strike, dip, rake, area, onset], rate


def write(path, sources):
    lines = [
        "moment tensor",
        "0.0 1.0 0.0",
        "1.0 0.0 0.0",
        "0.0 0.0 0.0",
        "solid velocity component",
        "0.0 0.0 1.0",
        "number of point sources",
        str(len(sources)),
        "point source parameters",
    ]
    lines += [
        " ".join(repr(value) for value in parameters) for parameters, _ in sources
    ]
    lines += ["source function: dt, number of samples", f"{TIMESTEP!r} {SAMPLES}"]
    lines += ["samples"]
    lines += [" ".join(f"{value:.9e}" for value in rate) for _, rate in sources]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("directory", type=Path)
    parser.add_argument("simulations", type=int, nargs="+")
    args = parser.parse_args()
    args.directory.mkdir(parents=True, exist_ok=True)
    for count in args.simulations:
        # the sources of a count do not depend on the other counts written
        generator = random.Random(count)
        sources = [source(generator, index) for index in range(count)]
        write(args.directory / f"source-s{count}.fsrm", sources)
        write(args.directory / f"source-copies-s{count}.fsrm", sources + sources)


if __name__ == "__main__":
    main()
