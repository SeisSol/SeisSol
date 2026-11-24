#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Include ratchet for the layout of src/ into Core/, WavePropagation/, Features/ and Simulation/.

Every file below src/ belongs to a unit and has a role, both given by its path; a file that has
not moved yet counts with the path that .ci/restructure.py gives it.

  Core/<X>/...                       core of unit Core/<X>
  WavePropagation/...                the wave propagation: always there, as one of its equations
                                     and schemes; the base of the features
  Features/<F>/...                   feature F: can be left out, and the simulation still runs
  .../Physics/...                    physics of the wave propagation or of a feature
  .../{Solver,Output}/...            solver of it
  .../Setup/...                      setup of it
  .../<file>                         root of it (configuration, basic types)
  Simulation/..., Proxy/..., src/*   composition

An include breaks a rule if
  - the core includes the wave propagation, a feature or the composition,
  - the physics includes more of the core than Core/Common and Core/Numerical, or the solver or
    setup, or an MPI, YAML, easi, device, or file format library,
  - the solver includes the setup, or the root the solver or setup,
  - the wave propagation or a feature includes the composition,
  - the wave propagation includes a feature, or a feature another feature,
  - or a unit of the core includes one that comes later in CORE_ORDER (the core has no cycles).
A feature may include the wave propagation (its physics only the physics and the root). The
composition may include anything. Together, the rules leave no cycles between the units.

Every rule break must be listed in the baseline (.ci/layering-baseline.txt), and the baseline must
not list rule breaks that are gone: after removing some, run with --update. New rule breaks only
go into the baseline with --update --accept-new.

usage: layering.py [--update [--accept-new]] [--report]
"""

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from includegraph import Sources, read_includes, repo_root  # noqa: E402

try:
    from restructure import target  # noqa: E402
except ImportError:  # the layout is complete

    def target(path, batches=None):
        return path


BASELINE = ".ci/layering-baseline.txt"
# (split, so that REUSE does not take the strings for the license of this file)
LICENSE_HEADER = (
    "# SPDX"
    "-FileCopyrightText: 2026 SeisSol Group\n"
    "#\n"
    "# SPDX"
    "-License-Identifier: BSD-3-Clause\n"
    "# SPDX"
    "-LicenseComments: Full text under /LICENSE and /LICENSES/\n"
    "#\n"
    "# SPDX"
    "-FileContributor: Author lists in /AUTHORS and /CITATION.cff\n"
)
# the units of the core, from the bottom; a unit may include itself and the units before it
CORE_ORDER = (
    ("src", "Core/Common"),
    ("Core/Data",),
    ("Core/Parallel",),
    ("Core/Geometry",),
    ("Core/Numerical",),
    ("Core/Memory",),
    ("Core/Kernels",),
    ("Core/Modules",),
    ("Core/IO",),
    ("Core/Monitoring",),
    ("Core/Run",),
    ("Core/TimeStepping",),
)
CORE_LEVEL = {unit: level for level, units in enumerate(CORE_ORDER) for unit in units}
WAVEPROPAGATION = "WavePropagation"
ROLES = {"Physics": "physics", "Solver": "solver", "Output": "solver", "Setup": "setup"}
PHYSICS_CORE = ("Core/Common", "Core/Numerical")
RUNTIME_CORE = (
    "Core/Parallel",
    "Core/Memory",
    "Core/Monitoring",
    "Core/Modules",
    "Core/IO",
    "Core/Kernels",
    "Core/Run",
    "Core/TimeStepping",
)
LIBRARIES = ("mpi.h", "yaml-cpp/", "easi/", "Device/", "asagi.h", "hdf5.h", "netcdf.h")


def classify(path):
    """(unit, role) of a path below src/ in the target layout."""
    parts = path.split("/")
    if parts[0] == "Core" and len(parts) > 2:
        return "/".join(parts[:2]), "core"
    if parts[0] == WAVEPROPAGATION and len(parts) > 1:
        role = ROLES.get(parts[1], "root") if len(parts) > 2 else "root"
        return WAVEPROPAGATION, role
    if parts[0] == "Features" and len(parts) > 2:
        role = ROLES.get(parts[2], "root") if len(parts) > 3 else "root"
        return "/".join(parts[:2]), role
    if len(parts) == 1:
        return "src", "core"
    if parts[0] in ("Simulation", "Proxy"):
        return parts[0], "composition"
    return parts[0], "unplaced"


def kind(unit, role):
    """core, wavepropagation, feature, composition or unplaced."""
    if role in ("core", "composition", "unplaced"):
        return role
    return "wavepropagation" if unit == WAVEPROPAGATION else "feature"


def rule_break(source, destination):
    """The rule that an include from `source` to `destination` ((unit, role) each) breaks."""
    source_unit, source_role = source
    destination_unit, destination_role = destination
    if source_role in ("composition", "unplaced") or destination_role == "unplaced":
        return None
    source_kind, destination_kind = kind(*source), kind(*destination)
    if source_kind == "core":
        if destination_kind == "core":
            if CORE_LEVEL.get(destination_unit, 0) > CORE_LEVEL.get(source_unit, 0):
                return "core-includes-later-core"
            return None
        return f"core-includes-{destination_kind}"
    if destination_kind == "composition":
        return f"{source_kind}-includes-composition"
    if destination_kind == "core":
        if source_role == "physics" and destination_unit not in PHYSICS_CORE:
            return "physics-includes-runtime"
        return None
    if destination_unit != source_unit:
        if destination_kind == "feature":
            return f"{source_kind}-includes-feature"
        # a feature includes the wave propagation
        if source_role == "physics" and destination_role not in ("physics", "root"):
            return "physics-includes-solver"
        return None
    if source_role in ("physics", "root") and destination_role in ("solver", "setup"):
        return f"{source_role}-includes-{destination_role}"
    if source_role == "solver" and destination_role == "setup":
        return "solver-includes-setup"
    return None


class Graph:
    def __init__(self, root):
        self.root = root
        sources = Sources(root)
        self.files = sorted(sources.src)
        self.target = {path: target(path) for path in self.files}
        self.edges = defaultdict(
            set
        )  # file -> files (both relative to src/, current paths)
        self.libraries = defaultdict(set)  # file -> external includes of LIBRARIES
        for path in self.files:
            for _, bracket, spelled in read_includes(root, "src/" + path):
                resolved = None
                if bracket == '"':
                    resolved, _ = sources.resolve("src/" + path, spelled)
                if resolved is not None:
                    self.edges[path].add(resolved)
                elif spelled.startswith(LIBRARIES):
                    self.libraries[path].add(spelled)

    def where(self, path):
        return classify(self.target[path])

    def rule_breaks(self):
        breaks = set()
        for path in self.files:
            source = self.where(path)
            for included in self.edges[path]:
                rule = rule_break(source, self.where(included))
                if rule is not None:
                    breaks.add((rule, self.target[path], self.target[included]))
            if source[1] == "physics":
                for library in self.libraries[path]:
                    breaks.add(("physics-includes-library", self.target[path], library))
        return breaks

    def closure(self, path):
        seen, stack = set(), [path]
        while stack:
            for included in self.edges[stack.pop()]:
                if included not in seen:
                    seen.add(included)
                    stack.append(included)
        return seen


def strongly_connected(edges):
    """The strongly connected components of a graph {node: {nodes}} (iterative Tarjan)."""
    nodes = set(edges) | {node for targets in edges.values() for node in targets}
    index, low, on_stack, stack, components = {}, {}, set(), [], []
    counter = 0
    for start in sorted(nodes):
        if start in index:
            continue
        work = [(start, iter(sorted(edges.get(start, ()))))]
        index[start] = low[start] = counter
        counter += 1
        stack.append(start)
        on_stack.add(start)
        while work:
            node, children = work[-1]
            advanced = False
            for child in children:
                if child not in index:
                    index[child] = low[child] = counter
                    counter += 1
                    stack.append(child)
                    on_stack.add(child)
                    work.append((child, iter(sorted(edges.get(child, ())))))
                    advanced = True
                    break
                if child in on_stack:
                    low[node] = min(low[node], index[child])
            if advanced:
                continue
            work.pop()
            if work:
                parent = work[-1][0]
                low[parent] = min(low[parent], low[node])
            if low[node] == index[node]:
                component = set()
                while True:
                    member = stack.pop()
                    on_stack.discard(member)
                    component.add(member)
                    if member == node:
                        break
                components.append(frozenset(component))
    return components


def read_baseline(root):
    path = os.path.join(root, BASELINE)
    if not os.path.exists(path):
        return set()
    entries = set()
    with open(path, encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith("#"):
                entries.add(tuple(line.split()))
    return entries


def write_baseline(root, breaks):
    with open(os.path.join(root, BASELINE), "w", encoding="utf-8") as file:
        file.write(
            LICENSE_HEADER + "#\n"
            "# The include rule breaks that .ci/layering.py tolerates (rule, includer, included;\n"
            "# paths in the target layout). Only ever shrinks: run .ci/layering.py --update after\n"
            "# removing rule breaks.\n"
        )
        for entry in sorted(breaks):
            file.write(" ".join(entry) + "\n")


def report(graph, breaks):
    lines = []
    counts = defaultdict(int)
    for rule, _, _ in breaks:
        counts[rule] += 1
    lines.append(f"rule breaks: {len(breaks)}")
    for rule in sorted(counts):
        lines.append(f"  {rule}: {counts[rule]}")

    unit_edges = defaultdict(set)
    for path in graph.files:
        for included in graph.edges[path]:
            source, destination = graph.where(path)[0], graph.where(included)[0]
            if source != destination:
                unit_edges[source].add(destination)
    units = {graph.where(path)[0] for path in graph.files}
    components = sorted(
        (c for c in strongly_connected(unit_edges) if len(c) > 1), key=len, reverse=True
    )
    largest = len(components[0]) if components else 1
    lines.append(f"units: {len(units)}, largest cycle: {largest} units")
    for component in components:
        lines.append(f"  cycle: {' '.join(sorted(component))}")

    physics = [path for path in graph.files if graph.where(path)[1] == "physics"]
    impure = 0
    for path in physics:
        reached = graph.closure(path) | {path}
        if any(
            graph.where(other)[0] in RUNTIME_CORE
            or graph.where(other)[1] in ("composition", "unplaced")
            or graph.libraries[other]
            for other in reached
        ):
            impure += 1
    lines.append(
        f"physics files whose includes reach the runtime, the composition or a library: "
        f"{impure} of {len(physics)}"
    )

    yaml = {
        path
        for path in graph.files
        if any("yaml" in lib for lib in graph.libraries[path])
    }
    cpps = [path for path in graph.files if path.endswith((".cpp", ".cu"))]
    reaching = sum(1 for path in cpps if (graph.closure(path) | {path}) & yaml)
    lines.append(f"translation units that include yaml-cpp: {reaching} of {len(cpps)}")

    unplaced = sorted(
        {
            graph.where(path)[0]
            for path in graph.files
            if graph.where(path)[1] == "unplaced"
        }
    )
    if unplaced:
        lines.append(f"directories outside the layout: {' '.join(unplaced)}")
    print("\n".join(lines))
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a", encoding="utf-8") as file:
            file.write("## Layering\n\n```\n" + "\n".join(lines) + "\n```\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--update", action="store_true", help="rewrite the baseline")
    parser.add_argument(
        "--accept-new", action="store_true", help="with --update: add new rule breaks"
    )
    parser.add_argument("--report", action="store_true", help="print the metrics")
    args = parser.parse_args()
    root = repo_root()
    graph = Graph(root)
    breaks = graph.rule_breaks()
    baseline = read_baseline(root)
    new = sorted(breaks - baseline)
    gone = sorted(baseline - breaks)
    if args.report:
        report(graph, breaks)
    if args.update:
        if new and not args.accept_new:
            print(f"{len(new)} new rule breaks; add them with --accept-new:")
            for entry in new:
                print("  " + " ".join(entry))
            return 1
        write_baseline(root, breaks)
        print(f"baseline: {len(breaks)} rule breaks ({len(gone)} gone, {len(new)} new)")
        return 0
    status = 0
    if new:
        print(
            f"{len(new)} new include rule breaks (see the docstring of .ci/layering.py):"
        )
        for rule, includer, included in new:
            print(f"  {rule}: {includer} includes {included}")
        status = 1
    if gone:
        print(
            f"{len(gone)} rule breaks of the baseline are gone; "
            "run .ci/layering.py --update and commit the baseline"
        )
        status = 1
    if status == 0:
        print(f"layering: {len(breaks)} known rule breaks, none new")
    return status


if __name__ == "__main__":
    sys.exit(main())
