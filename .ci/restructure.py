#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Moves the sources of SeisSol into the layout of Core/, WavePropagation/, Features/ and
Simulation/.

A batch is a list of rules (old prefix, new prefix) on the paths below src/ as they are when the
batch runs; a path takes the rule with the longest matching prefix. The batches run in order, and
`target()` gives the path a file has after all of them (layering.py checks the includes by it).

The script can be replayed on any branch: it moves the files that exist there (files that only
the branch has go by the rules of their directories), rewrites the includes in src/, app/ and
tests/, updates the CMake source lists, and rewrites mentions of the moved paths elsewhere
(codegen/, the CMake and clang-tidy configuration, docs/, comments). Afterwards, run
`.ci/header.py --fix src app tests` (the include guards follow the paths) and clang-format, e.g.
with `pre-commit run --all-files`; both rewrite files, so never build at the same time.

usage:
  restructure.py plan [--batch NAME ...] [--tree]
  restructure.py apply --batch NAME [--batch NAME ...]
  restructure.py target PATH ...
"""

import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from includegraph import Sources, read_includes, repo_root  # noqa: E402

WP = "WavePropagation"
DR = "Features/DynamicRupture"
PL = "Features/Plasticity"
ST = "Features/SourceTerm"
IC = "Features/InitialConditions"
TM = "Features/TimeMirror"
SIM = "Simulation"

EQUATIONS = (
    "acoustic",
    "anisotropic",
    "elastic",
    "poroelastic",
    "viscoacoustic",
    "viscoelastic",
)

CORE = [
    (f"{name}/", f"Core/{name}/")
    for name in (
        "Common",
        "Numerical",
        "Geometry",
        "Parallel",
        "Memory",
        "Monitoring",
        "Modules",
        "IO",
    )
]

LAYOUT = [
    # core: the low-level headers of Initializer/, the kernel infrastructure, the time stepping
    ("Initializer/BasicTypedefs.h", "Core/Common/BasicTypedefs.h"),
    ("Initializer/PreProcessorMacros.h", "Core/Common/PreProcessorMacros.h"),
    ("Initializer/InputAux.h", "Core/Common/InputAux.h"),
    ("Solver/MultipleSimulations.h", "Core/Common/MultipleSimulations.h"),
    ("Initializer/CellLocalInformation.h", "Core/Memory/CellLocalInformation.h"),
    ("Initializer/LtsSetup.h", "Core/Memory/LtsSetup.h"),
    ("Initializer/Typedefs.h", "Core/Memory/Typedefs.h"),
    ("Initializer/DeviceGraph.h", "Core/Parallel/DeviceGraph.h"),
    ("Initializer/FaceMap.", "Core/Geometry/FaceMap."),
    ("Initializer/PointMapper.", "Core/Geometry/PointMapper."),
    (
        "Initializer/BatchRecorders/DataTypes/",
        "Core/Kernels/BatchRecorders/DataTypes/",
    ),
    (
        "Initializer/BatchRecorders/Recorders.h",
        "Core/Kernels/BatchRecorders/Recorders.h",
    ),
    ("Kernels/Common.h", "Core/Kernels/Common.h"),
    ("Kernels/Kernel.h", "Core/Kernels/Kernel.h"),
    ("Kernels/MemoryOps.h", "Core/Kernels/MemoryOps.h"),
    ("Kernels/Touch.", "Proxy/Touch."),
    *[
        (f"Kernels/DeviceAux/{name}.h", f"Core/Parallel/DeviceAux/{name}.h")
        for name in ("ClockAux", "GatherAux", "HostCopyAux")
    ],
    *[
        (
            f"Kernels/DeviceAux/{backend}/{name}.cpp",
            f"Core/Parallel/DeviceAux/{backend}/{name}.cpp",
        )
        for backend in ("cudahip", "sycl")
        for name in ("ClockAux", "GatherAux", "HostCopyAux")
    ],
    ("Solver/TimeStepping/", "Core/TimeStepping/"),
    ("Solver/Settings.h", "Core/Memory/SimulationSettings.h"),
    (
        "Core/Monitoring/ActorStateStatistics.",
        "Core/TimeStepping/Actor/ActorStateStatistics.",
    ),
    ("Core/Monitoring/Unit.", "Core/Common/Unit."),
    ("Core/Monitoring/Instrumentation.h", "Core/Common/Instrumentation.h"),
    ("Core/Monitoring/Trace.", "Core/Parallel/Trace."),
    ("Initializer/Clustering/", "Core/TimeStepping/Clustering/"),
    ("Initializer/TimeStepping/", "Core/TimeStepping/Layout/"),
    (
        "ResultWriter/ClusteringWriter.",
        "Core/TimeStepping/Clustering/ClusteringWriter.",
    ),
    ("ResultWriter/MiniSeisSolWriter.", "Core/Monitoring/MiniSeisSolWriter."),
    ("ResultWriter/ThreadsPinningWriter.", "Core/Monitoring/ThreadsPinningWriter."),
    ("ResultWriter/AsyncIO.", "Core/IO/AsyncIO."),
    ("ResultWriter/PickpointWriter.", "Core/IO/PickpointWriter."),
    # the composition
    ("Initializer/InitProcedure/", f"{SIM}/InitProcedure/"),
    ("Initializer/Parameters/", f"{SIM}/Parameters/"),
    # the files of the simulations of a run, from which the parts read their models
    ("Initializer/Parameters/SimulationFiles.h", "Core/Data/SimulationFiles.h"),
    ("Initializer/MemoryManager.", f"{SIM}/MemoryManager."),
    ("Initializer/ParameterDB.", f"{SIM}/ParameterDB."),
    ("Initializer/Preflight", f"{SIM}/Preflight"),
    ("Reader/", "Core/IO/Easi/"),
    # the data table and the easi reader of davschneller/easi-analytic, for a replay there
    ("Reader/Scripting/DataTable.", "Core/Data/DataTable."),
    ("Reader/Scripting/DataReader.", "Core/Data/DataReader."),
    ("Reader/Scripting/EasiReader.", "Core/IO/Easi/EasiReader."),
    ("Reader/Scripting/ReaderBuilder.", "Core/IO/Easi/ReaderBuilder."),
    ("SeisSol.", f"{SIM}/SeisSol."),
    ("Solver/Simulator.", f"{SIM}/Simulator."),
    ("Solver/Estimator.", f"{SIM}/Estimator."),
    ("Kernels/BuiltConfigs.cpp", f"{SIM}/BuiltConfigs.cpp"),
    ("ResultWriter/EnergyOutput.", f"{SIM}/EnergyOutput."),
    # the wave propagation: the base of the features, with the equations and the schemes
    ("Equations/", f"{WP}/Physics/"),
    *[(f"Equations/{name}/Model/", f"{WP}/Physics/{name}/") for name in EQUATIONS],
    ("Model/", f"{WP}/Physics/"),
    # the bridges between the equations and the schemes
    ("Equations/Setup.h", f"{WP}/Solver/Setup.h"),
    ("Model/Concepts.h", f"{WP}/Solver/Concepts.h"),
    ("Kernels/LinearCK/", f"{WP}/Solver/LinearCK/"),
    ("Kernels/LinearCKAnelastic/", f"{WP}/Solver/LinearCKAnelastic/"),
    ("Kernels/STP/", f"{WP}/Solver/STP/"),
    *[
        (f"Kernels/{name}", f"{WP}/Solver/{name}")
        for name in (
            "Interface.h",
            "Data.h",
            "Solver.h",
            "SolverSelector.h",
            "Local.h",
            "Neighbor.h",
            "Time.h",
            "Spacetime.h",
            "TimeCommon.",
            "ConfigBoundary.",
            "AnalyticalBoundary.h",
        )
    ],
    (
        "Initializer/BatchRecorders/LocalIntegrationRecorder.cpp",
        f"{WP}/Solver/LocalIntegrationRecorder.cpp",
    ),
    (
        "Initializer/BatchRecorders/NeighIntegrationRecorder.cpp",
        f"{WP}/Solver/NeighIntegrationRecorder.cpp",
    ),
    ("Initializer/Model/CellLocalMatrices.", f"{WP}/Setup/CellLocalMatrices."),
    ("Initializer/Model/BoundaryMappings.", f"{WP}/Setup/BoundaryMappings."),
    ("Initializer/BoundaryHelper.h", f"{WP}/Solver/BoundaryHelper.h"),
    ("Initializer/BoundarySetup.h", f"{WP}/Setup/BoundarySetup.h"),
    ("Kernels/Receiver.", f"{WP}/Output/Receiver."),
    ("ResultWriter/ReceiverWriter.", f"{WP}/Output/ReceiverWriter."),
    ("ResultWriter/DofSync.", f"{WP}/Output/DofSync."),
    ("Solver/FreeSurfaceIntegrator.", f"{WP}/Output/FreeSurfaceIntegrator."),
    # dynamic rupture
    ("DynamicRupture/Misc.", f"{DR}/Misc."),
    ("DynamicRupture/Typedefs.h", f"{DR}/Typedefs.h"),
    ("DynamicRupture/Factory.", f"{DR}/Setup/Factory."),
    ("DynamicRupture/FrictionLaws/", f"{DR}/Solver/FrictionLaws/"),
    (
        "DynamicRupture/FrictionLaws/RateAndStateCommon.h",
        f"{DR}/Physics/RateAndStateCommon.h",
    ),
    ("DynamicRupture/FrictionLaws/TPCommon.h", f"{DR}/Physics/TPCommon.h"),
    ("DynamicRupture/FrictionLaws/Dual.h", f"{DR}/Physics/Dual.h"),
    ("DynamicRupture/Initializer/", f"{DR}/Setup/Initializer/"),
    ("DynamicRupture/Output/", f"{DR}/Output/"),
    ("Kernels/DynamicRupture.", f"{DR}/Solver/DynamicRupture."),
    (
        "Solver/TimeStepping/Compute/DynamicRuptureCluster.",
        f"{DR}/Solver/DynamicRuptureCluster.",
    ),
    ("Initializer/Model/DynamicRupture", f"{DR}/Setup/DynamicRupture"),
    (
        "Initializer/BatchRecorders/DynamicRuptureRecorder.cpp",
        f"{DR}/Solver/DynamicRuptureRecorder.cpp",
    ),
    # plasticity
    ("Model/Plasticity", f"{PL}/Physics/Plasticity"),
    ("Kernels/Plasticity.", f"{PL}/Solver/Plasticity."),
    ("Kernels/DeviceAux/PlasticityAux.h", f"{PL}/Solver/DeviceAux/PlasticityAux.h"),
    *[
        (
            f"Kernels/DeviceAux/{backend}/PlasticityAux.cpp",
            f"{PL}/Solver/DeviceAux/{backend}/PlasticityAux.cpp",
        )
        for backend in ("cudahip", "sycl")
    ],
    (
        "Initializer/BatchRecorders/PlasticityRecorder.cpp",
        f"{PL}/Solver/PlasticityRecorder.cpp",
    ),
    ("Core/Memory/Descriptor/Plasticity.h", f"{PL}/Solver/Storage.h"),
    ("Core/Memory/Descriptor/DynamicRupture.h", f"{DR}/Solver/Storage.h"),
    # point sources
    ("SourceTerm/", f"{ST}/Setup/"),
    ("SourceTerm/Typedefs.h", f"{ST}/Typedefs.h"),
    ("SourceTerm/PointSource.", f"{ST}/Physics/PointSource."),
    ("Kernels/PointSourceCluster", f"{ST}/Solver/PointSourceCluster"),
    # initial conditions
    ("Physics/InitialField.h", f"{WP}/Physics/InitialField.h"),
    ("Physics/Scenario/", f"{IC}/Physics/Scenario/"),
    ("Initializer/InitialFieldProjection.", f"{IC}/Setup/InitialFieldProjection."),
    ("ResultWriter/AnalysisWriter.", f"{IC}/Output/AnalysisWriter."),
    # the instantaneous time mirror
    (
        "Physics/InstantaneousTimeMirrorManager.",
        f"{TM}/Solver/InstantaneousTimeMirrorManager.",
    ),
]

BATCHES = [("core", CORE), ("layout", LAYOUT)]

# includes that the build only resolved through the include directory of the equation of the
# first configuration (removed by the layout batch): what they name instead
ALIASES = {"layout": {"Model/Datastructures.h": "Equations/Datastructures.h"}}

# CMake blocks the batches remove from src/sources.cmake (from the first line to the next endif())
REMOVED_CMAKE_BLOCKS = {
    "layout": ["# the material headers of the equation of the first configuration"],
}

# directories whose CMakeLists.txt holds an OBJECT library; Core/ and Features/ group them
UNIT_GROUPS = ("Core", "Features")

# where the script rewrites mentions of moved paths (besides the includes)
TEXT_DIRS = ("src", "app", "tests", "codegen", "cmake", "docs", ".github", ".ci")
TEXT_FILES = ("CMakeLists.txt", ".clang-tidy")
TEXT_SKIP = (".ci/restructure.py", ".ci/layering-baseline.txt")
TEXT_SUFFIXES = (
    ".h",
    ".hpp",
    ".cuh",
    ".cpp",
    ".cu",
    ".in",
    ".py",
    ".cmake",
    ".txt",
    ".rst",
    ".md",
    ".par",
    ".yml",
    ".yaml",
    ".sh",
    ".tidy",
)

# (split, so that REUSE does not take the strings for the license of this file)
CMAKE_HEADER = (
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


def rule_target(rules, path):
    best = None
    for old, new in rules:
        if path.startswith(old) and (best is None or len(old) > len(best[0])):
            best = (old, new)
    if best is None:
        return path
    return best[1] + path[len(best[0]) :]


def target(path, batches=None):
    """The path below src/ that `path` has after the given batches (default: all)."""
    for name, rules in BATCHES:
        if batches is None or name in batches:
            path = rule_target(rules, path)
    return path


def git(root, *args):
    subprocess.run(["git", "-C", root, *args], check=True)


def unit_of(path):
    """The directory of the OBJECT library that compiles `path` (relative to src/), or None for
    the files of src/ itself (src/sources.cmake lists them)."""
    parts = path.split("/")
    if parts[0] in UNIT_GROUPS and len(parts) > 2:
        return "/".join(parts[:2])
    if len(parts) > 1:
        return parts[0]
    return None


def library_name(unit):
    return "seissol-" + unit.split("/")[-1].lower()


class Unit:
    """The CMakeLists.txt of a directory with an OBJECT library."""

    STANDARD = re.compile(r"^\s*target_link_libraries\(\S+ PUBLIC \S+\)\s*$")

    def __init__(self, unit, text):
        self.unit = unit
        self.library = None
        self.sources = []
        self.conditional = []  # [(condition, [sources])]
        self.other = []
        lines = text.splitlines()
        index = 0
        # the license header
        while index < len(lines) and (
            lines[index].startswith("#") or not lines[index].strip()
        ):
            index += 1
        while index < len(lines):
            line = lines[index]
            match = re.match(r"^\s*add_library\((\S+)\s+OBJECT\b(.*)$", line)
            if match:
                self.library = match.group(1)
                rest = match.group(2)
                while True:
                    done = ")" in rest
                    self.sources += rest.replace(")", " ").split()
                    if done:
                        break
                    index += 1
                    rest = lines[index]
                index += 1
                continue
            match = re.match(r"^\s*if\s*\((.*)\)\s*$", line)
            if match:
                block = []
                index += 1
                while not re.match(r"^\s*endif\s*\(\s*\)\s*$", lines[index]):
                    block.append(lines[index])
                    index += 1
                files = []
                for entry in block:
                    inner = re.match(
                        r"^\s*target_sources\(\S+\s+PRIVATE\s+(.*)\)\s*$", entry
                    )
                    if inner is None:
                        raise RuntimeError(
                            f"{unit}/CMakeLists.txt: cannot carry over: {entry}"
                        )
                    files += inner.group(1).split()
                self.conditional.append((match.group(1), files))
                index += 1
                continue
            if not self.STANDARD.match(line) and (line.strip() or self.other):
                self.other.append(line)
            index += 1
        while self.other and not self.other[-1].strip():
            self.other.pop()

    def render(self):
        text = CMAKE_HEADER + "\n"
        text += f"add_library({self.library} OBJECT\n"
        for source in sorted(self.sources):
            text += f"    {source}\n"
        text += ")\n\n"
        for condition, files in self.conditional:
            text += f"if ({condition})\n"
            for source in files:
                text += f"    target_sources({self.library} PRIVATE {source})\n"
            text += "endif()\n\n"
        text += (
            f"target_link_libraries({self.library} PUBLIC seissol-common-properties)\n"
        )
        text += f"target_link_libraries(seissol-lib PUBLIC {self.library})\n"
        if self.other:
            text += "\n" + "\n".join(self.other) + "\n"
        return text


def find_units(root):
    """The directories (relative to src/) with a CMakeLists.txt that holds an OBJECT library."""
    units = {}
    srcdir = os.path.join(root, "src")
    for dirpath, _, names in os.walk(srcdir):
        if "CMakeLists.txt" in names and dirpath != srcdir:
            relative = os.path.relpath(dirpath, srcdir)
            with open(
                os.path.join(dirpath, "CMakeLists.txt"), encoding="utf-8"
            ) as file:
                text = file.read()
            if re.search(r"add_library\(\S+\s+OBJECT", text):
                units[relative] = Unit(relative, text)
    return units


def compute_moves(sources, batches):
    moves = {}
    for path in sorted(sources.src):
        new = target(path, batches)
        if new != path:
            moves[path] = new
    targets = defaultdict(list)
    for old, new in moves.items():
        targets[new].append(old)
    clashes = {new: olds for new, olds in targets.items() if len(olds) > 1}
    clashes.update(
        {new: [new] for new in targets if new in sources.src and new not in moves}
    )
    if clashes:
        for new, olds in sorted(clashes.items()):
            print(f"clash: {new} <- {', '.join(olds)}", file=sys.stderr)
        raise SystemExit(1)
    return moves


def plan(args, root):
    sources = Sources(root)
    batches = args.batch or None
    moves = compute_moves(sources, batches)
    if args.tree:
        counts = defaultdict(int)
        for path in sources.src:
            directory = os.path.dirname(target(path, batches))
            counts[directory] += 1
            while directory:
                directory = os.path.dirname(directory)
                counts.setdefault(directory, 0)
        for directory in sorted(counts, key=lambda d: d.split("/")):
            depth = directory.count("/") + 1 if directory else 0
            name = os.path.basename(directory) or "src"
            print(
                f"{'  ' * depth}{name}/"
                + (f"  {counts[directory]}" if counts[directory] else "")
            )
    else:
        for old, new in sorted(moves.items()):
            print(f"{old} -> {new}")
    print(f"{len(moves)} of {len(sources.src)} files move", file=sys.stderr)
    if batches is None:
        unmapped = sorted(
            path
            for path in sources.src
            if "/" in target(path)
            and target(path).split("/")[0] not in (*UNIT_GROUPS, SIM, "Proxy")
        )
        if unmapped:
            print(f"{len(unmapped)} files stay outside the layout:", file=sys.stderr)
            for path in unmapped:
                print(f"  {path}", file=sys.stderr)


def rewrite_includes(root, sources, moves, aliases):
    """Rewrites the includes of all C++ files (at their new places) to the new paths."""
    count = 0
    for includer in sources.files:
        new_includer = includer
        if includer.startswith("src/") and includer[4:] in moves:
            new_includer = "src/" + moves[includer[4:]]
        # the files have moved already
        includes = read_includes(root, new_includer)
        changes = {}
        for index, bracket, spelled in includes:
            if bracket != '"':
                continue
            old_target, relative = sources.resolve(includer, spelled)
            if (
                old_target is None
                and spelled in aliases
                and aliases[spelled] in sources.src
            ):
                old_target = aliases[spelled]
            if old_target is None:
                continue
            new_target = moves.get(old_target, old_target)
            if new_target == old_target and new_includer == includer:
                continue
            if relative:
                here = os.path.dirname(new_includer[4:])
                if os.path.normpath(os.path.join(here, spelled)) == new_target:
                    continue
            elif new_target == os.path.normpath(spelled):
                continue
            changes[index] = (spelled, new_target)
        if not changes:
            continue
        path = os.path.join(root, new_includer)
        with open(path, encoding="utf-8") as file:
            lines = file.readlines()
        seen = set()
        result = []
        for index, line in enumerate(lines):
            if index in changes:
                spelled, new_target = changes[index]
                line = line.replace(f'"{spelled}"', f'"{new_target}"', 1)
                count += 1
            key = line.split("//")[0].strip()
            if key.startswith("#include") and key in seen:
                continue
            if key.startswith("#include"):
                seen.add(key)
            result.append(line)
        with open(path, "w", encoding="utf-8") as file:
            file.writelines(result)
    return count


def update_cmake(root, old_units, moves, batches):
    """Moves the sources between the OBJECT libraries and regroups the directories."""
    srcdir = os.path.join(root, "src")
    new_units = {}

    def moved_unit_dir(unit):
        # an OBJECT library moves along with its directory, if the directory moves as a whole
        probe = target(f"{unit}/CMakeLists.txt", batches)
        directory = os.path.dirname(probe)
        return directory if unit_of(probe) == directory else None

    for unit, data in old_units.items():
        destination = moved_unit_dir(unit)
        if destination is not None and destination not in new_units:
            carried = Unit.__new__(Unit)
            carried.unit = destination
            carried.library = (
                data.library if destination == unit else library_name(destination)
            )
            carried.sources, carried.conditional, carried.other = [], [], data.other
            new_units[destination] = carried

    def place(old_source, condition=None):
        new_source = moves.get(old_source, old_source)
        unit = unit_of(new_source)
        if unit not in new_units:
            created = Unit.__new__(Unit)
            created.unit = unit
            created.library = library_name(unit)
            created.sources, created.conditional, created.other = [], [], []
            new_units[unit] = created
        relative = os.path.relpath(new_source, unit)
        if condition is None:
            new_units[unit].sources.append(relative)
        else:
            for existing_condition, files in new_units[unit].conditional:
                if existing_condition == condition:
                    files.append(relative)
                    break
            else:
                new_units[unit].conditional.append((condition, [relative]))

    for unit, data in old_units.items():
        for source in data.sources:
            place(os.path.normpath(os.path.join(unit, source)))
        for condition, files in data.conditional:
            for source in files:
                place(os.path.normpath(os.path.join(unit, source)), condition)
        if data.other and moved_unit_dir(unit) is None:
            raise RuntimeError(
                f"{unit}/CMakeLists.txt dissolves, but has: {data.other}"
            )

    # the CMake files of the directories that move as a whole go along
    for unit in old_units:
        destination = moved_unit_dir(unit)
        if destination is not None and destination != unit:
            os.makedirs(os.path.join(srcdir, destination), exist_ok=True)
            git(
                root,
                "mv",
                f"src/{unit}/CMakeLists.txt",
                f"src/{destination}/CMakeLists.txt",
            )
    for unit in old_units:
        if moved_unit_dir(unit) is None and unit not in new_units:
            git(root, "rm", "-q", f"src/{unit}/CMakeLists.txt")

    for unit, data in new_units.items():
        if not data.sources and not data.conditional:
            path = f"src/{unit}/CMakeLists.txt"
            if os.path.exists(os.path.join(root, path)):
                git(root, "rm", "-q", path)
            continue
        path = os.path.join(srcdir, unit, "CMakeLists.txt")
        text = data.render()
        if os.path.exists(path):
            with open(path, encoding="utf-8") as file:
                old = Unit(unit, file.read())
            if (
                sorted(old.sources) == sorted(data.sources)
                and old.conditional == data.conditional
                and old.library == data.library
            ):
                continue
        with open(path, "w", encoding="utf-8") as file:
            file.write(text)
        git(root, "add", f"src/{unit}/CMakeLists.txt")

    # the groups and src/CMakeLists.txt add the directories with a CMakeLists.txt
    for group in UNIT_GROUPS:
        groupdir = os.path.join(srcdir, group)
        if not os.path.isdir(groupdir):
            continue
        children = sorted(
            name
            for name in os.listdir(groupdir)
            if os.path.exists(os.path.join(groupdir, name, "CMakeLists.txt"))
        )
        with open(
            os.path.join(groupdir, "CMakeLists.txt"), "w", encoding="utf-8"
        ) as file:
            file.write(CMAKE_HEADER + "\n")
            file.write(
                "\n\n".join(f"add_subdirectory({child})" for child in children) + "\n"
            )
        git(root, "add", f"src/{group}/CMakeLists.txt")
    top = sorted(
        name
        for name in os.listdir(srcdir)
        if os.path.exists(os.path.join(srcdir, name, "CMakeLists.txt"))
    )
    path = os.path.join(srcdir, "CMakeLists.txt")
    with open(path, encoding="utf-8") as file:
        lines = file.read().splitlines()
    positions = [
        i for i, line in enumerate(lines) if line.startswith("add_subdirectory(")
    ]
    # one empty line between the entries, as src/CMakeLists.txt asks for
    entries = "\n\n".join(f"add_subdirectory({name})" for name in top).split("\n")
    lines[positions[0] : positions[-1] + 1] = entries
    with open(path, "w", encoding="utf-8") as file:
        file.write("\n".join(lines) + "\n")


def remove_cmake_blocks(root, batches):
    path = os.path.join(root, "src", "sources.cmake")
    with open(path, encoding="utf-8") as file:
        lines = file.read().splitlines()
    for batch in batches:
        for first in REMOVED_CMAKE_BLOCKS.get(batch, []):
            if first not in lines:
                continue
            start = lines.index(first)
            end = start
            while not re.match(r"^\s*endif\s*\(\s*\)\s*$", lines[end]):
                end += 1
            del lines[start : end + 1]
            while (
                start < len(lines)
                and start > 0
                and not lines[start].strip()
                and not lines[start - 1].strip()
            ):
                del lines[start]
    with open(path, "w", encoding="utf-8") as file:
        file.write("\n".join(lines) + "\n")


def rewrite_mentions(root, moves):
    """Rewrites mentions of moved paths (relative to src/, optionally with src/ in front) in the
    sources, the configuration and the docs; links (preceded by a slash) stay."""
    if not moves:
        return 0
    names = sorted(moves, key=len, reverse=True)
    pattern = re.compile(
        r"(?<![\w./-])(src/)?("
        + "|".join(re.escape(name) for name in names)
        + r")(?![\w/-])"
    )
    count = 0
    files = [
        path
        for path in subprocess.run(
            ["git", "-C", root, "ls-files", "-z", "--", *TEXT_DIRS, *TEXT_FILES],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.split("\0")
        if path
        and (path.endswith(TEXT_SUFFIXES) or os.path.basename(path) in TEXT_FILES)
        and path not in TEXT_SKIP
        and not (path.startswith("src/") and os.path.basename(path) == "CMakeLists.txt")
    ]
    for path in files:
        full = os.path.join(root, path)
        if not os.path.isfile(full):
            continue
        with open(full, encoding="utf-8", errors="surrogateescape") as file:
            text = file.read()
        new_text, replaced = pattern.subn(
            lambda match: (match.group(1) or "") + moves[match.group(2)], text
        )
        if replaced:
            count += replaced
            with open(full, "w", encoding="utf-8", errors="surrogateescape") as file:
                file.write(new_text)
    return count


def apply(args, root):
    batches = args.batch
    sources = Sources(root)
    moves = compute_moves(sources, batches)
    aliases = {}
    for batch in batches:
        aliases.update(ALIASES.get(batch, {}))
    old_units = find_units(root)
    print(f"moving {len(moves)} files")
    for old, new in sorted(moves.items()):
        os.makedirs(os.path.join(root, "src", os.path.dirname(new)), exist_ok=True)
        git(root, "mv", f"src/{old}", f"src/{new}")
    count = rewrite_includes(root, sources, moves, aliases)
    print(f"rewrote {count} includes")
    update_cmake(root, old_units, moves, batches)
    remove_cmake_blocks(root, batches)
    count = rewrite_mentions(root, moves)
    print(f"rewrote {count} other mentions of moved paths")
    # git does not track directories; the ones the moves left empty go
    for dirpath, _, _ in sorted(os.walk(os.path.join(root, "src")), reverse=True):
        if not os.listdir(dirpath):
            os.rmdir(dirpath)
    print("next: .ci/header.py --fix src app tests, then clang-format (pre-commit)")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    plan_parser = sub.add_parser("plan")
    plan_parser.add_argument(
        "--batch", action="append", choices=[n for n, _ in BATCHES]
    )
    plan_parser.add_argument("--tree", action="store_true")
    apply_parser = sub.add_parser("apply")
    apply_parser.add_argument(
        "--batch", action="append", choices=[n for n, _ in BATCHES], required=True
    )
    target_parser = sub.add_parser("target")
    target_parser.add_argument("path", nargs="+")
    args = parser.parse_args()
    root = repo_root()
    if args.command == "plan":
        plan(args, root)
    elif args.command == "apply":
        apply(args, root)
    else:
        for path in args.path:
            print(target(path))


if __name__ == "__main__":
    main()
