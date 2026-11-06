#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Find the commits from which on the result of a case changes.

The script builds SeisSol at commits of a range, runs one case with each build, and
compares the outputs of the runs with compare-history.py. It runs every stride-th commit
of the range first; wherever two neighbouring runs differ, it runs the commit in the
middle, until each change sits between two adjacent commits. A change is thus assumed to
stay: a result that changes and changes back between two runs goes unnoticed.

    bisect-history.py --source ~/seissol-build-tree --case precomputed-seissol/tpv105 \\
        --range v1.3.1..origin/master --work ~/history-tpv105 \\
        --cmake-arg=-DHOST_ARCH=skx --cmake-arg=-DGEMM_TOOLS_LIST=LIBXSMM

--source is a clone of SeisSol that the script checks out at each commit; it must not be
the clone this script runs from. Of the first-parent commits of the range, only those
that change what is built go into the search (see relevance); --commits takes a list of
its own instead. The build is the one configuration of the case (by default elastic,
order 6, double precision, the Stroud rule on the fault), with the same further CMake
arguments at every commit: keep the host architecture and the GEMM generators fixed, or
they show up as changes. The case directory is copied per run, next to a link to the
meshes directory of precomputed-seissol, so that "../meshes/..." resolves.

An old commit may not build with the tools of today, e.g. a code generator that a newer
Python rejects. --fixup names a shell command that runs in the source after each
checkout (with the commit in BISECT_COMMIT) and patches what needs it; the build log
records what it changed, and the report lists the commits that it changed.

Everything lands in --work: the candidates (commits.tsv), the state of every commit
(state.json, so that a second call resumes), the logs, the outputs, the changes found
(changes.md, changes.json) and a manifest of all runs for compare-history.py, which runs
on it at the end (history/).

Two runs count as equal if no metric of compare-history.py exceeds its tolerance: the
nrms of any station or map (--tolerance, in percent), the largest difference of the
rupture time (--rupture-tolerance, in ms) and the relative difference of the seismic
moment (--moment-tolerance; some versions write the energies with six digits only).
--repeat runs the first commit twice and reports what two runs of one build differ by.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)


def _load_sibling(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


history = _load_sibling("compare_history", "compare-history.py")

# --- which commits can change what is built -------------------------------------------

# paths that do not take part in a build or a run of SeisSol
UNBUILT = re.compile(
    r"^(docs/|Documentation/|\.github/|\.ci/|\.gitlab|tests/|scripts/|postprocessing/"
    r"|preprocessing/|README|CHANGELOG|CITATION|AUTHORS|LICENSE|LICENSES/|\.clang"
    r"|\.pre-commit|\.flake8|\.reuse|REUSE|\.gitignore|\.git-blame|CONTRIBUTING"
    r"|CODE_OF_CONDUCT|Dockerfile|\.dockerignore|\.cmake-format)"
)
TESTS = re.compile(r"(\.t\.h$|/tests?/|_test\.|Tests?/)")
# code that a CPU build of an elastic case with a rate-and-state fault does not run, or
# that only writes what the comparison does not read (the free surface, the flop counts)
OTHER_CODE = re.compile(
    r"(GpuImpl/|\.cuh?$|Sycl|SYCL|sycl|Hip[A-Z.]|/hip|Cuda|CUDA|cuda|DeviceAux"
    r"|OpenMPOffload|[Aa]nisotropic|[Pp]oroelastic|[Vv]iscoelastic|[Aa]coustic|Damage"
    r"|Plasticity|ImposedSlip|LinearSlipWeakening|Checkpoint|pybind|Monitoring/"
    r"|Proxy/|FreeSurface)"
)
RUPTURE = re.compile(
    r"(DynamicRupture|dynamic_rupture|Friction|ThermalPressuri|RateAndState)",
    re.IGNORECASE,
)


def relevance(paths: list[str]) -> str:
    """How a commit that changes these paths can bear on the case.

    "none": it changes nothing that is built; "unlikely": only code that the case does
    not run; "rupture": the dynamic rupture; "core": anything else that is built.
    """
    built = [
        path for path in paths if not UNBUILT.match(path) and not TESTS.search(path)
    ]
    if not built:
        return "none"
    if all(OTHER_CODE.search(path) for path in built):
        return "unlikely"
    if any(RUPTURE.search(path) for path in built):
        return "rupture"
    return "core"


def git(source: str, *args: str, check: bool = True) -> str:
    result = subprocess.run(
        ["git", "-C", source, *args], capture_output=True, text=True, check=False
    )
    if check and result.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)}: {result.stderr.strip()}")
    return result.stdout


@dataclass
class Commit:
    sha: str
    subject: str
    date: str
    relevance: str

    @property
    def short(self) -> str:
        return self.sha[:9]


def candidates(source: str, revisions: str, classes: set[str]) -> list[Commit]:
    """The first-parent commits of a range, from its start on, that are relevant.

    The start of the range is the first candidate, so that the search has a run to
    begin from.
    """
    start, _, _ = revisions.partition("..")
    shas = git(source, "rev-list", "--first-parent", "--reverse", revisions).split()
    if start:
        shas = [git(source, "rev-parse", f"{start}^{{commit}}").strip()] + shas
    commits = []
    for index, sha in enumerate(shas):
        subject = git(source, "log", "-1", "--format=%ad%x09%s", "--date=short", sha)
        date, _, subject = subject.strip().partition("\t")
        if index == 0 and start:
            commits.append(Commit(sha, subject, date, "start"))
            continue
        paths = git(source, "diff", "--name-only", f"{sha}^1", sha).split()
        kind = relevance(paths)
        if kind in classes:
            commits.append(Commit(sha, subject, date, kind))
    return commits


def read_commits(source: str, path: str) -> list[Commit]:
    """Commits listed one per line (a revision, optionally followed by anything)."""
    commits = []
    with open(path) as handle:
        for line in handle:
            revision = line.split("#")[0].split()
            if not revision:
                continue
            sha = git(source, "rev-parse", f"{revision[0]}^{{commit}}").strip()
            date, _, subject = (
                git(source, "log", "-1", "--format=%ad%x09%s", "--date=short", sha)
                .strip()
                .partition("\t")
            )
            commits.append(Commit(sha, subject, date, "listed"))
    return commits


# --- the search -----------------------------------------------------------------------


@dataclass
class SearchResult:
    evaluated: dict  # index -> whether the run succeeded
    changes: list  # (index before, index of the change)
    unresolved: list  # (index before, index after): a change among failed commits


def search(count: int, evaluate, differs, stride: int) -> SearchResult:
    """Locate every change in a sequence of count states that can be evaluated.

    evaluate(i) runs state i and says whether that worked; differs(i, j) compares two
    evaluated states. Every stride-th state, and the last, is evaluated first; then the
    state nearest the middle of each pair of neighbouring evaluated states that differ,
    until the pairs that differ are adjacent, or all states between them failed.
    """
    evaluated: dict[int, bool] = {}

    def run(index):
        if index not in evaluated:
            evaluated[index] = bool(evaluate(index))
        return evaluated[index]

    if count == 0:
        return SearchResult(evaluated, [], [])
    for index in sorted(set(range(0, count, max(1, stride))) | {count - 1}):
        run(index)

    while True:
        good = [i for i in sorted(evaluated) if evaluated[i]]
        open_pairs = [
            (a, b)
            for a, b in zip(good, good[1:])
            if b - a > 1
            and any(k not in evaluated for k in range(a + 1, b))
            and differs(a, b)
        ]
        if not open_pairs:
            break
        for a, b in open_pairs:
            middle = (a + b) / 2.0
            pending = [k for k in range(a + 1, b) if k not in evaluated]
            run(min(pending, key=lambda k: (abs(k - middle), k)))

    good = [i for i in sorted(evaluated) if evaluated[i]]
    changes, unresolved = [], []
    for a, b in zip(good, good[1:]):
        if differs(a, b):
            (changes if b - a == 1 else unresolved).append((a, b))
    return SearchResult(evaluated, changes, unresolved)


# --- comparing two runs ---------------------------------------------------------------


@dataclass
class Tolerances:
    nrms: float = 1e-9  # percent
    rupture_ms: float = 1e-6
    moment: float = 1e-5


def distance(a, b, region, spacing: float = 100.0, cache: dict | None = None) -> dict:
    """The largest differences between two runs (Solutions of compare-history.py).

    cache keeps the fault data of a run on the grid, by the id of the run.
    """
    grid = history.Grid.over(region, spacing)
    cache = {} if cache is None else cache
    for solution in (a, b):
        if solution.id not in cache:
            cache[solution.id] = history.grid_fault(solution, grid)
    gridded = {a.id: cache[a.id], b.id: cache[b.id]}
    settings = history.Settings(math.inf, 0.0, {})
    rows = history.compare_solutions(a, b, gridded, settings)

    def largest(select):
        values = [row["value"] for row in rows if select(row)]
        return max(values) if values else None

    return {
        "station_nrms": largest(
            lambda row: row["category"] in ("fault", "body") and row["metric"] == "nrms"
        ),
        "map_nrms": largest(
            lambda row: row["category"] == "fault-map" and row["metric"] == "nrms"
        ),
        "rt_max_ms": largest(
            lambda row: row["category"] == "rupture" and row["metric"] == "max_ms"
        ),
        "rt_rms_ms": largest(
            lambda row: row["category"] == "rupture" and row["metric"] == "rms_ms"
        ),
        "rt_bias_ms": largest(
            lambda row: row["category"] == "rupture" and row["metric"] == "bias_ms"
        ),
        "rt_area_mismatch": largest(
            lambda row: row["category"] == "rupture"
            and row["metric"] == "area_mismatch"
        ),
        "m0_rel": largest(
            lambda row: row["category"] == "moment" and row["metric"] == "rel"
        ),
        "compared": len(rows),
    }


def differ(measured: dict, tolerances: Tolerances) -> bool:
    if not measured["compared"]:
        return True

    def above(key, limit):
        value = measured.get(key)
        return value is not None and (not math.isfinite(value) or abs(value) > limit)

    return (
        above("station_nrms", tolerances.nrms)
        or above("map_nrms", tolerances.nrms)
        or above("rt_max_ms", tolerances.rupture_ms)
        or above("rt_area_mismatch", 0.0)
        or above("m0_rel", tolerances.moment)
    )


# --- building and running -------------------------------------------------------------

DEFAULT_CMAKE_ARGS = (
    "-DCMAKE_BUILD_TYPE=Release",
    "-DEQUATIONS=elastic",
    "-DORDER=6",
    "-DPRECISION=double",
    "-DNUMBER_OF_MECHANISMS=0",
    "-DDR_QUAD_RULE=stroud",
    "-DTESTING=OFF",
)
# the target of the executable, as it has been named
TARGETS = ("seissol-bin", "SeisSol-bin")
_NOT_THE_EXECUTABLE = re.compile(
    r"(lib|proxy|test|gpucode|offloading|device|\.so|\.a$)"
)


def find_executable(build: str) -> str | None:
    """The SeisSol executable a build produced, by its name."""
    found = []
    for root, _, files in os.walk(build):
        for name in files:
            path = os.path.join(root, name)
            if (
                re.match(r"^(SeisSol_|seissol)", name)
                and not _NOT_THE_EXECUTABLE.search(name.lower())
                and os.access(path, os.X_OK)
            ):
                found.append(path)
    return max(found, key=os.path.getmtime) if found else None


@dataclass
class Campaign:
    source: str
    case: str
    work: str
    cmake_args: list
    jobs: int
    launcher: str
    ranks: int
    threads: int
    build_timeout: float
    run_timeout: float
    output_prefix: str = "output/tpv"
    commits: list = field(default_factory=list)
    state: dict = field(default_factory=dict)
    fixups: list = field(default_factory=list)
    # the paths that the fixups changed at a commit, by its sha
    patched: dict = field(default_factory=dict)

    @property
    def build_dir(self) -> str:
        return os.path.join(self.work, "build")

    def log(self, commit: Commit, kind: str) -> str:
        os.makedirs(os.path.join(self.work, "logs"), exist_ok=True)
        return os.path.join(self.work, "logs", f"{commit.short}.{kind}.log")

    def run_dir(self, commit: Commit, suffix: str = "") -> str:
        return os.path.join(self.work, "runs", commit.short + suffix)

    def output(self, commit: Commit, suffix: str = "") -> str:
        case = os.path.basename(self.case.rstrip("/"))
        return os.path.join(self.run_dir(commit, suffix), case, self.output_prefix)

    def save(self) -> None:
        path = os.path.join(self.work, "state.json")
        with open(path + ".new", "w") as handle:
            json.dump(self.state, handle, indent=1)
        os.replace(path + ".new", path)

    def _call(self, command, log, timeout, cwd=None, env=None) -> bool:
        with open(log, "a") as handle:
            handle.write(f"$ {shlex.join(command)}\n")
            handle.flush()
            try:
                result = subprocess.run(
                    command,
                    cwd=cwd,
                    env=env,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    timeout=timeout,
                    check=False,
                )
            except subprocess.TimeoutExpired:
                handle.write(f"timed out after {timeout} s\n")
                return False
        return result.returncode == 0

    def checkout(self, commit: Commit, log: str) -> bool:
        # the clone is the script's own: what a commit does not track goes, also the
        # submodules that another commit has and this one has not
        steps = (
            ["git", "-C", self.source, "checkout", "--force", "--detach", commit.sha],
            ["git", "-C", self.source, "clean", "-ffdxq"],
            ["git", "-C", self.source, "submodule", "sync", "--recursive"],
            [
                "git",
                "-C",
                self.source,
                "submodule",
                "update",
                "--init",
                "--recursive",
                "--force",
            ],
            [
                "git",
                "-C",
                self.source,
                "submodule",
                "foreach",
                "--recursive",
                "git clean -ffdxq",
            ],
        )
        if not all(self._call(step, log, 3600) for step in steps):
            return False
        return self.fix(commit, log)

    def fix(self, commit: Commit, log: str) -> bool:
        """Run the fixups in the checked-out source and note what they changed."""
        self.patched.pop(commit.sha, None)
        if not self.fixups:
            return True
        env = dict(os.environ, BISECT_COMMIT=commit.sha)
        for command in self.fixups:
            if not self._call(["/bin/sh", "-c", command], log, 3600, self.source, env):
                return False
        status = git(self.source, "status", "--porcelain", "--ignore-submodules=none")
        paths = sorted(
            {line[3:].strip() for line in status.splitlines() if line.strip()}
        )
        with open(log, "a") as handle:
            handle.write(f"fixups changed: {', '.join(paths) or 'nothing'}\n")
        if paths:
            self.patched[commit.sha] = paths
        return True

    def build(self, commit: Commit) -> str | None:
        """Build the executable at a commit; a failed build is tried once more anew."""
        log = self.log(commit, "build")
        if not self.checkout(commit, log):
            return None
        for attempt in ("incremental", "clean"):
            if attempt == "clean":
                shutil.rmtree(self.build_dir, ignore_errors=True)
            # an executable left from another commit must not pass for this one
            while (stale := find_executable(self.build_dir)) is not None:
                os.remove(stale)
            generator = ["-G", "Ninja"] if shutil.which("ninja") else []
            launcher = (
                [
                    "-DCMAKE_CXX_COMPILER_LAUNCHER=ccache",
                    "-DCMAKE_C_COMPILER_LAUNCHER=ccache",
                ]
                if shutil.which("ccache")
                else []
            )
            configured = self._call(
                ["cmake", "-S", self.source, "-B", self.build_dir]
                + generator
                + launcher
                + list(DEFAULT_CMAKE_ARGS)
                + list(self.cmake_args),
                log,
                self.build_timeout,
            )
            if not configured:
                continue
            for target in TARGETS:
                built = self._call(
                    ["cmake", "--build", self.build_dir, "--target", target]
                    + ["-j", str(self.jobs)],
                    log,
                    self.build_timeout,
                )
                if built:
                    executable = find_executable(self.build_dir)
                    if executable:
                        return executable
        return None

    def prepare_case(self, commit: Commit, suffix: str = "") -> str:
        """A copy of the case for a run, next to the meshes it refers to."""
        directory = self.run_dir(commit, suffix)
        shutil.rmtree(directory, ignore_errors=True)
        os.makedirs(directory)
        name = os.path.basename(self.case.rstrip("/"))
        target = os.path.join(directory, name)
        shutil.copytree(
            self.case, target, ignore=shutil.ignore_patterns("precomputed", "output")
        )
        os.makedirs(
            os.path.dirname(os.path.join(target, self.output_prefix)), exist_ok=True
        )
        meshes = os.path.join(os.path.dirname(os.path.abspath(self.case)), "meshes")
        if os.path.isdir(meshes):
            os.symlink(meshes, os.path.join(directory, "meshes"))
        return target

    def run(self, commit: Commit, executable: str, suffix: str = "") -> bool:
        case = self.prepare_case(commit, suffix)
        env = dict(os.environ, OMP_NUM_THREADS=str(self.threads))
        command = shlex.split(self.launcher.format(ranks=self.ranks)) + [
            executable,
            "parameters.par",
        ]
        log = self.log(commit, "run" + suffix)
        ran = self._call(command, log, self.run_timeout, case, env)
        output = os.path.dirname(os.path.join(case, self.output_prefix))
        return ran and bool(os.listdir(output))

    def evaluate(self, commit: Commit) -> bool:
        record = self.state.get(commit.sha, {})
        if "ok" in record:
            return record["ok"]
        started = time.time()
        print(f"[{time.strftime('%H:%M:%S')}] {commit.short} {commit.subject[:60]}")
        executable = self.build(commit)
        record = {"built": executable is not None}
        if commit.sha in self.patched:
            record["patched"] = self.patched[commit.sha]
        record["ok"] = executable is not None and self.run(commit, executable)
        record["seconds"] = round(time.time() - started)
        print(f"    {'ok' if record['ok'] else 'failed'} after {record['seconds']} s")
        self.state[commit.sha] = record
        self.save()
        return record["ok"]


# --- report ---------------------------------------------------------------------------


def _number(value, digits=3) -> str:
    if value is None:
        return "–"
    return "0" if value == 0 else f"{value:.{digits}g}"


def write_report(
    work,
    commits,
    result: SearchResult,
    measured: dict,
    repeat,
    left_out: dict,
    patched: dict | None = None,
) -> None:
    """changes.json and changes.md; left_out counts, per change, the commits between
    its two candidates that the search left out (as irrelevant), any of which may be
    the change instead. patched holds the paths that the fixups changed, by the sha of
    the commit; a change between two commits that the fixups changed differently may
    come from the fixups."""
    patched = patched or {}

    def fixups_differ(a, b):
        return set(patched.get(commits[a].sha, ())) != set(
            patched.get(commits[b].sha, ())
        )

    changes = []
    for a, b in result.changes:
        values = measured[(a, b)]
        changes.append(
            {
                "before": commits[a].sha,
                "at": commits[b].sha,
                "left_out_before": left_out.get((a, b), 0),
                "fixups_differ": fixups_differ(a, b),
                **values,
            }
        )
    unresolved = [
        {
            "after": commits[a].sha,
            "until": commits[b].sha,
            "failed": [commits[k].sha for k in range(a + 1, b)],
        }
        for a, b in result.unresolved
    ]
    with open(os.path.join(work, "changes.json"), "w") as handle:
        json.dump(
            {
                "changes": changes,
                "unresolved": unresolved,
                "repeat": repeat,
                "patched": {
                    commits[i].sha: patched[commits[i].sha]
                    for i in sorted(result.evaluated)
                    if commits[i].sha in patched
                },
            },
            handle,
            indent=1,
        )

    runs = sum(1 for ok in result.evaluated.values() if ok)
    failed = [commits[i] for i, ok in sorted(result.evaluated.items()) if not ok]
    lines = [
        "# Commits that change the result",
        "",
        f"{len(commits)} candidates, {runs} runs, {len(failed)} failed builds or runs, "
        f"{len(changes)} changes located, {len(unresolved)} only up to failed commits.",
        "",
        "Against the run of the candidate before; station and map nrms in percent, the "
        "rupture time in ms, M0 relative. A change can also come from one of the "
        "commits left out before it, which the search took for irrelevant.",
        "",
    ]
    if repeat is not None:
        lines += [
            "Two runs of the first build differ by: "
            + ", ".join(f"{key} {_number(value)}" for key, value in repeat.items())
            + "",
            "",
        ]
    head = [
        "commit",
        "date",
        "subject",
        "kind",
        "left out before",
        "station nrms",
        "map nrms",
        "RT max",
        "RT bias",
        "M0",
    ]
    lines += ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for a, b in result.changes:
        values, commit = measured[(a, b)], commits[b]
        lines.append(
            "| "
            + " | ".join(
                [
                    f"`{commit.short}`",
                    commit.date,
                    commit.subject.replace("|", "/"),
                    commit.relevance
                    + (", fixups differ" if fixups_differ(a, b) else ""),
                    str(left_out.get((a, b), 0)),
                    _number(values["station_nrms"]),
                    _number(values["map_nrms"]),
                    _number(values["rt_max_ms"]),
                    _number(values["rt_bias_ms"]),
                    _number(values["m0_rel"]),
                ]
            )
            + " |"
        )
    if unresolved:
        lines += ["", "## Changes among commits that failed", ""]
        for entry in unresolved:
            lines.append(
                f"- after `{entry['after'][:9]}`, up to `{entry['until'][:9]}`; "
                f"failed: {', '.join(f'`{sha[:9]}`' for sha in entry['failed'])}"
            )
    if failed:
        lines += ["", "## Failed builds or runs", ""]
        lines += [f"- `{c.short}` {c.date} {c.subject}" for c in failed]
    fixed = [commits[i] for i in sorted(result.evaluated) if commits[i].sha in patched]
    if fixed:
        lines += ["", "## Commits that the fixups changed", ""]
        lines += [
            f"- `{c.short}` {c.date} {c.subject}: {', '.join(patched[c.sha])}"
            for c in fixed
        ]
    with open(os.path.join(work, "changes.md"), "w") as handle:
        handle.write("\n".join(lines) + "\n")


def write_manifest(work, campaign: Campaign, commits, result, region) -> str:
    """A manifest of compare-history.py with every run in the order of the commits."""
    runs = []
    for index in sorted(result.evaluated):
        if not result.evaluated[index]:
            continue
        commit = commits[index]
        number = re.search(r"#(\d+)", commit.subject)
        runs.append(
            {
                "id": commit.short,
                "label": f"#{number.group(1)}" if number else commit.short,
                "output": campaign.output(commit),
                "commit": commit.short,
                "date": commit.date,
                "note": commit.subject,
            }
        )
    manifest = {
        "title": f"{os.path.basename(campaign.case.rstrip('/'))}: commits that change it",
        "region": list(region),
        "runs": runs,
    }
    path = os.path.join(work, "manifest.json")
    with open(path, "w") as handle:
        json.dump(manifest, handle, indent=1)
    return path


# --- main -----------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        epilog="Details at the top of this script.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--source", required=True, help="clone of SeisSol to build")
    parser.add_argument("--case", required=True, help="directory of the case")
    parser.add_argument("--work", required=True, help="directory for everything else")
    parser.add_argument("--range", help="commits to search, e.g. v1.3.1..origin/master")
    parser.add_argument("--commits", help="file listing the commits to search instead")
    parser.add_argument(
        "--kinds",
        default="rupture,core",
        help="relevance of the commits of a range that are searched: "
        "rupture, core, unlikely, none (default: rupture,core)",
    )
    parser.add_argument("--stride", type=int, default=8, help="run every n-th first")
    parser.add_argument("--cmake-arg", action="append", default=[], dest="cmake_args")
    parser.add_argument(
        "--fixup",
        action="append",
        default=[],
        dest="fixups",
        help="shell command run in the source after each checkout (with the commit "
        "in BISECT_COMMIT), e.g. a script that patches what today's tools reject",
    )
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    parser.add_argument(
        "--launcher",
        default="mpirun -n {ranks} --bind-to none",
        help="command that starts SeisSol, with {ranks}; without binding, the "
        "threads of a single rank are not all pinned to one core",
    )
    parser.add_argument(
        "--output-prefix",
        default="output/tpv",
        help="OutputFile of the case, relative to its directory",
    )
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--build-timeout", type=float, default=4 * 3600.0)
    parser.add_argument("--run-timeout", type=float, default=4 * 3600.0)
    parser.add_argument(
        "--region",
        type=float,
        nargs=4,
        default=(-15e3, 15e3, 0.0, 15e3),
        metavar=("ALONG_MIN", "ALONG_MAX", "DOWN_MIN", "DOWN_MAX"),
        help="region of the rupture time comparison [m]",
    )
    parser.add_argument("--tolerance", type=float, default=Tolerances.nrms)
    parser.add_argument(
        "--rupture-tolerance", type=float, default=Tolerances.rupture_ms
    )
    parser.add_argument("--moment-tolerance", type=float, default=Tolerances.moment)
    parser.add_argument(
        "--repeat", action="store_true", help="run the first commit twice"
    )
    parser.add_argument(
        "--list", action="store_true", help="only write and print the candidates"
    )
    parser.add_argument(
        "--no-history", action="store_true", help="skip compare-history.py at the end"
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    source = os.path.abspath(args.source)
    if os.path.realpath(source) == os.path.realpath(os.path.join(HERE, "..", "..")):
        print("error: --source must be another clone than this one", file=sys.stderr)
        return 2
    if bool(args.range) == bool(args.commits):
        print("error: give either --range or --commits", file=sys.stderr)
        return 2
    os.makedirs(args.work, exist_ok=True)
    work = os.path.abspath(args.work)

    if args.range:
        kinds = {kind.strip() for kind in args.kinds.split(",")}
        commits = candidates(source, args.range, kinds)
    else:
        commits = read_commits(source, args.commits)
    with open(os.path.join(work, "commits.tsv"), "w") as handle:
        for commit in commits:
            handle.write(
                f"{commit.sha}\t{commit.date}\t{commit.relevance}\t{commit.subject}\n"
            )
    print(f"{len(commits)} candidates, written to {os.path.join(work, 'commits.tsv')}")
    if args.list or not commits:
        return 0

    state_path = os.path.join(work, "state.json")
    state = {}
    if os.path.isfile(state_path):
        with open(state_path) as handle:
            state = json.load(handle)
    campaign = Campaign(
        source=source,
        case=os.path.abspath(args.case),
        work=work,
        cmake_args=args.cmake_args,
        jobs=args.jobs,
        launcher=args.launcher,
        ranks=args.ranks,
        threads=args.threads,
        build_timeout=args.build_timeout,
        run_timeout=args.run_timeout,
        output_prefix=args.output_prefix,
        commits=commits,
        state=state,
        fixups=args.fixups,
    )
    tolerances = Tolerances(
        args.tolerance, args.rupture_tolerance, args.moment_tolerance
    )
    axes = history.Axes.from_spec(None)
    solutions = {}

    def solution(index):
        if index not in solutions:
            commit = commits[index]
            solutions[index] = history.load_run(
                {"id": commit.short, "output": campaign.output(commit)}, axes, work
            )
        return solutions[index]

    measured, cache = {}, {}

    def differs(a, b):
        if (a, b) not in measured:
            measured[(a, b)] = distance(
                solution(a), solution(b), args.region, cache=cache
            )
        return differ(measured[(a, b)], tolerances)

    repeat = None
    if args.repeat and campaign.evaluate(commits[0]):
        # a second run of a build of the first commit, beside the first run
        executable = campaign.build(commits[0])
        if executable and campaign.run(commits[0], executable, "-repeat"):
            again = history.load_run(
                {"id": "repeat", "output": campaign.output(commits[0], "-repeat")},
                axes,
                work,
            )
            repeat = distance(solution(0), again, args.region)

    result = search(
        len(commits), lambda i: campaign.evaluate(commits[i]), differs, args.stride
    )
    left_out = {}
    for a, b in result.changes:
        count = git(
            source,
            "rev-list",
            "--first-parent",
            "--count",
            f"{commits[a].sha}..{commits[b].sha}",
            check=False,
        ).strip()
        left_out[(a, b)] = max(int(count or "1") - 1, 0)
    patched = {
        sha: record["patched"]
        for sha, record in campaign.state.items()
        if record.get("patched")
    }
    write_report(work, commits, result, measured, repeat, left_out, patched)
    manifest = write_manifest(work, campaign, commits, result, args.region)
    print(f"changes: {os.path.join(work, 'changes.md')}")
    if not args.no_history:
        history.main([manifest, "--out", os.path.join(work, "history")])
    return 0


if __name__ == "__main__":
    sys.exit(main())
