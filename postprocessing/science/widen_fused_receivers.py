#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Convert the text receivers of a fused run to the layout SeisSol wrote before.

SeisSol writes a receiver of a fused run as a row per sample and simulation, with
the simulation, counted from zero, in the column SimulationIndex. Before, it wrote
a row per sample and a column per quantity and simulation, with the simulation in
the name: counted from zero and appended in the off-fault receivers (v10, v11,
...), counted from one after a dash in the on-fault ones (SRs-1, SRs-2, ...),
whose header had a line P_0<n>, T_s<n> and T_d<n> per simulation. This script
writes that layout, for tools that read it.

The values are carried over as they are written, without reading them as numbers
and printing them anew. A file of a run of a single simulation, which has no
SimulationIndex, is copied as it is.

usage: widen_fused_receivers.py OUTPUT_DIRECTORY RECEIVER_FILE [RECEIVER_FILE ...]
"""

import argparse
import re
import shutil
from pathlib import Path

STRESS_LINE = re.compile(r"#\s*(P_0|T_s|T_d)\s+(.*)$")


def widen(lines, fault):
    """The lines of a receiver file in the layout of a column per simulation."""
    names = re.findall(r'"([^"]*)"', lines[1])
    column = names.index("SimulationIndex")
    quantities = [name for i, name in enumerate(names) if i not in (0, column)]

    comments = []
    rows = []
    for line in lines[2:]:
        if line.startswith("#"):
            comments.append(line)
        elif line.strip():
            rows.append(line.split())
    simulations = sorted({int(row[column]) for row in rows})

    def suffix(simulation):
        return f"-{simulation + 1}" if fault else str(simulation)

    wide = ["Time"] + [
        quantity + suffix(simulation)
        for simulation in simulations
        for quantity in quantities
    ]
    separator = " ," if fault else ","
    output = [lines[0], "VARIABLES = " + separator.join(f'"{name}"' for name in wide)]

    for comment in comments:
        match = STRESS_LINE.match(comment)
        if fault and match:
            # a value per simulation, which took a line of its own
            for simulation, value in zip(simulations, match.group(2).split()):
                output.append(f"# {match.group(1)}{simulation + 1}\t{value}")
        else:
            output.append(comment)

    for start in range(0, len(rows), len(simulations)):
        sample = rows[start : start + len(simulations)]
        if [int(row[column]) for row in sample] != simulations or any(
            row[0] != sample[0][0] for row in sample
        ):
            raise ValueError(
                f"the simulations of a sample are incomplete at row {start}"
            )
        values = [sample[0][0]] + [
            value
            for row in sample
            for i, value in enumerate(row)
            if i not in (0, column)
        ]
        if fault:
            output.append("".join(value + "\t" for value in values))
        else:
            output.append("".join("  " + value for value in values))
    return output


def main():
    parser = argparse.ArgumentParser(
        description="Convert the text receivers of a fused run to a column per "
        "simulation, as SeisSol wrote them before."
    )
    parser.add_argument("output", type=Path, help="directory to write the files to")
    parser.add_argument("files", type=Path, nargs="+", help="receiver files")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    for path in args.files:
        target = args.output / path.name
        if target.resolve() == path.resolve():
            parser.error(
                f"{path} would be overwritten; choose another output directory"
            )
        lines = path.read_text().splitlines()
        if len(lines) < 2 or "SimulationIndex" not in lines[1]:
            shutil.copyfile(path, target)
            continue
        fault = "faultreceiver" in path.name
        target.write_text("\n".join(widen(lines, fault)) + "\n")


if __name__ == "__main__":
    main()
