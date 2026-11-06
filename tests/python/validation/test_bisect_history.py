# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Tests for scripts/validate/bisect-history.py

The search has to find every change of a sequence with few evaluations, and say where
failed builds hide one. The whole script runs on a small repository whose "SeisSol" is a
script that writes receivers, with a value that changes at known commits.
"""

import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

SEISSOL_ROOT = Path(__file__).resolve().parents[3]
_spec = importlib.util.spec_from_file_location(
    "bisect_history",
    SEISSOL_ROOT / "scripts" / "validate" / "bisect-history.py",
)
bh = importlib.util.module_from_spec(_spec)
sys.modules["bisect_history"] = bh
_spec.loader.exec_module(bh)


# ============================================================================
# relevance of the paths a commit changes
# ============================================================================


class TestRelevance:
    def test_documentation_and_tests_are_not_built(self):
        assert bh.relevance(["Documentation/x.rst", "tests/python/a.py"]) == "none"
        assert (
            bh.relevance(["src/Kernels/Time.t.h", ".github/workflows/a.yml"]) == "none"
        )

    def test_code_the_case_does_not_run(self):
        paths = [
            "src/DynamicRupture/FrictionLaws/GpuImpl/RateAndState.h",
            "src/Parallel/Runtime/StreamSycl.cpp",
            "src/Equations/poroelastic/Model/Datastructures.h",
        ]
        assert bh.relevance(paths) == "unlikely"

    def test_the_rupture_and_the_rest(self):
        assert (
            bh.relevance(["src/DynamicRupture/FrictionLaws/RateAndState.h"])
            == "rupture"
        )
        assert bh.relevance(["src/Kernels/Local.cpp", "docs/a.md"]) == "core"
        assert bh.relevance(["submodules/yateto"]) == "core"


# ============================================================================
# the search
# ============================================================================


def run_search(values, stride, failing=()):
    evaluations = []

    def evaluate(index):
        evaluations.append(index)
        return index not in failing

    result = bh.search(
        len(values), evaluate, lambda a, b: values[a] != values[b], stride
    )
    return result, evaluations


class TestSearch:
    def test_finds_every_change(self):
        values = [0] * 10 + [1] * 7 + [2] + [3] * 20 + [4] * 3
        truth = [i for i in range(1, len(values)) if values[i] != values[i - 1]]
        for stride in (1, 3, 8, 64):
            result, evaluations = run_search(values, stride)
            assert [b for _, b in result.changes] == truth
            assert all(b - a == 1 for a, b in result.changes)
            assert not result.unresolved
            assert len(evaluations) == len(set(evaluations))

    def test_evaluates_little_where_nothing_changes(self):
        result, evaluations = run_search([0] * 100, 10)
        assert not result.changes
        assert sorted(evaluations) == list(range(0, 100, 10)) + [99]

    def test_bisects_rather_than_walks(self):
        values = [0] * 60 + [1] * 60
        _, evaluations = run_search(values, 64)
        assert len(evaluations) < 12

    def test_a_failed_commit_is_stepped_around(self):
        values = [0] * 10 + [1] * 10
        result, _ = run_search(values, 4, failing={9})
        # the change is at 9 or 10; 9 does not build
        assert result.changes == []
        assert result.unresolved == [(8, 10)]

    def test_a_failure_away_from_the_change_does_no_harm(self):
        values = [0] * 10 + [1] * 10
        result, _ = run_search(values, 4, failing={4, 12})
        assert result.changes == [(9, 10)]
        assert not result.unresolved

    def test_trivial_sequences(self):
        assert bh.search(0, lambda i: True, lambda a, b: True, 4).changes == []
        result = bh.search(1, lambda i: True, lambda a, b: True, 4)
        assert result.evaluated == {0: True} and result.changes == []


# ============================================================================
# comparing two runs
# ============================================================================


def solution(name, scale=1.0):
    time = np.linspace(0.0, 2.0, 401)
    rate = scale * np.exp(-(((time - 1.0) / 0.1) ** 2))
    result = bh.history.Solution(id=name, label=name)
    result.stations["fault"]["faultst000dp075"] = bh.history.Station(
        "faultst000dp075",
        time,
        {"h-slip-rate": rate, "v-slip-rate": 0 * rate, "psi": 0.5 + rate},
    )
    return result


class TestDistance:
    region = (-1e3, 1e3, 0.0, 1e3)

    def test_equal_runs(self):
        measured = bh.distance(solution("a"), solution("b"), self.region)
        assert measured["station_nrms"] == 0.0
        assert not bh.differ(measured, bh.Tolerances())

    def test_a_tiny_difference_counts(self):
        measured = bh.distance(solution("a"), solution("b", 1 + 1e-9), self.region)
        assert 0.0 < measured["station_nrms"] < 1e-6
        assert bh.differ(measured, bh.Tolerances(nrms=1e-12))
        assert not bh.differ(measured, bh.Tolerances(nrms=1e-6))

    def test_nothing_to_compare_is_a_difference(self):
        empty = bh.history.Solution(id="empty", label="")
        measured = bh.distance(solution("a"), empty, self.region)
        assert measured["compared"] == 0 and bh.differ(measured, bh.Tolerances())


def test_the_executable_is_found_by_its_name(tmp_path):
    def make(name, executable=True):
        path = tmp_path / "app" / name
        path.parent.mkdir(exist_ok=True)
        path.write_text("x")
        if executable:
            path.chmod(0o755)
        return path

    make("libseissol.a")
    make("seissollib-elastic-o6-f64.so")
    make("seissol-proxy-exe")
    make("seissol-elastic-o6-f64-notes", executable=False)
    assert bh.find_executable(str(tmp_path)) is None
    binary = make("SeisSol_Release_dhsw_6_elastic")
    assert bh.find_executable(str(tmp_path)) == str(binary)


# ============================================================================
# the whole script, on a repository of its own
# ============================================================================

FAKE_SEISSOL = """#!{python}
import os
value = float("@VALUE@")
os.makedirs("output", exist_ok=True)
with open("output/tpv-faultreceiver-00001-00000.dat", "w") as handle:
    handle.write('TITLE = "Temporal Signal for fault receiver number 1"\\n')
    handle.write('VARIABLES = "Time","SRs","SRd","StV"\\n')
    handle.write("# x1 0.0\\n# x2 0.0\\n# x3 -7500.0\\n")
    for step in range(201):
        t = step * 0.01
        rate = value * 2.718281828 ** (-((t - 1.0) / 0.1) ** 2)
        handle.write(f"{{t:.6e}} {{rate:.15e}} 0.0 {{0.5 + rate:.15e}}\\n")
"""

CMAKE_LISTS = """cmake_minimum_required(VERSION 3.10)
project(fake NONE)
file(READ ${CMAKE_SOURCE_DIR}/seissol.in content)
file(READ ${CMAKE_SOURCE_DIR}/value.txt value)
string(REPLACE "@VALUE@" "${value}" content "${content}")
file(WRITE ${CMAKE_BINARY_DIR}/seissol-staging "${content}")
add_custom_target(seissol-bin ALL
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/seissol-staging
          ${CMAKE_BINARY_DIR}/seissol-elastic-o6-f64
  COMMAND chmod +x ${CMAKE_BINARY_DIR}/seissol-elastic-o6-f64)
"""


def git(repo, *args):
    subprocess.run(
        ["git", "-C", str(repo), *args],
        check=True,
        capture_output=True,
        env={
            **os.environ,
            "GIT_AUTHOR_NAME": "t",
            "GIT_AUTHOR_EMAIL": "t@t",
            "GIT_COMMITTER_NAME": "t",
            "GIT_COMMITTER_EMAIL": "t@t",
        },
    )


@pytest.mark.skipif(
    shutil.which("cmake") is None or shutil.which("git") is None,
    reason="needs cmake and git",
)
def test_end_to_end(tmp_path):
    repo = tmp_path / "source"
    repo.mkdir()
    git(repo, "init", "-q", "-b", "main")
    (repo / "CMakeLists.txt").write_text(CMAKE_LISTS)
    (repo / "seissol.in").write_text(FAKE_SEISSOL.format(python=sys.executable))
    (repo / "src").mkdir()
    values = [1.0, 1.0, 1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.0]
    for index, value in enumerate(values):
        (repo / "value.txt").write_text(str(value))
        (repo / "src" / "code.cpp").write_text(f"// revision {index}\n")
        git(repo, "add", "-A")
        git(repo, "commit", "-q", "-m", f"revision {index} (#{100 + index})")
        git(repo, "tag", f"r{index}")
    # a commit that builds nothing new is left out of the search
    (repo / "README.md").write_text("notes\n")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "notes")

    case = tmp_path / "precomputed" / "fake"
    case.mkdir(parents=True)
    (case / "parameters.par").write_text("&output\n/\n")
    work = tmp_path / "work"
    code = bh.main(
        [
            "--source",
            str(repo),
            "--case",
            str(case),
            "--work",
            str(work),
            "--range",
            "r0..main",
            "--stride",
            "3",
            "--launcher",
            "",
            "--threads",
            "1",
            "--jobs",
            "1",
            "--no-history",
            "--repeat",
        ]
    )
    assert code == 0
    commits = (work / "commits.tsv").read_text().splitlines()
    assert len(commits) == len(values)
    changes = json.loads((work / "changes.json").read_text())
    tags = {
        subprocess.run(
            ["git", "-C", str(repo), "rev-parse", f"r{i}"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip(): i
        for i in range(len(values))
    }
    assert [tags[change["at"]] for change in changes["changes"]] == [3, 7]
    assert all(change["station_nrms"] > 1.0 for change in changes["changes"])
    state = json.loads((work / "state.json").read_text())
    assert len(state) < len(values)
    assert "| `" in (work / "changes.md").read_text()
    manifest = json.loads((work / "manifest.json").read_text())
    assert [run["label"] for run in manifest["runs"]][0] == "#100"
    # the first commit ran twice, alike
    assert changes["repeat"]["station_nrms"] == 0.0
    first = commits[0].split()[0][:9]
    repeat_log = work / "logs" / f"{first}.run-repeat.log"
    assert repeat_log.read_text().count("$ ") == 1

    # a second call resumes: nothing is built or run again
    before = (work / "state.json").read_text()
    assert (
        bh.main(
            [
                "--source",
                str(repo),
                "--case",
                str(case),
                "--work",
                str(work),
                "--range",
                "r0..main",
                "--stride",
                "3",
                "--launcher",
                "",
                "--no-history",
                "--repeat",
            ]
        )
        == 0
    )
    assert (work / "state.json").read_text() == before
    assert repeat_log.read_text().count("$ ") == 1


# a build that fails while broken.txt says so, as an old commit may with today's tools
CMAKE_LISTS_CHECKED = CMAKE_LISTS.replace(
    "project(fake NONE)\n",
    "project(fake NONE)\n"
    "file(READ ${CMAKE_SOURCE_DIR}/broken.txt broken)\n"
    'if (broken MATCHES "yes")\n'
    '  message(FATAL_ERROR "this revision does not build")\n'
    "endif()\n",
)


@pytest.mark.skipif(
    shutil.which("cmake") is None or shutil.which("git") is None,
    reason="needs cmake and git",
)
def test_fixups_make_old_commits_build(tmp_path):
    repo = tmp_path / "source"
    repo.mkdir()
    git(repo, "init", "-q", "-b", "main")
    (repo / "CMakeLists.txt").write_text(CMAKE_LISTS_CHECKED)
    (repo / "seissol.in").write_text(FAKE_SEISSOL.format(python=sys.executable))
    (repo / "src").mkdir()
    revisions = [(1.0, "no"), (1.0, "yes"), (2.0, "yes"), (2.0, "no")]
    for index, (value, broken) in enumerate(revisions):
        (repo / "value.txt").write_text(str(value))
        (repo / "broken.txt").write_text(broken)
        (repo / "src" / "code.cpp").write_text(f"// revision {index}\n")
        git(repo, "add", "-A")
        git(repo, "commit", "-q", "-m", f"revision {index}")
        git(repo, "tag", f"r{index}")
    tags = {
        subprocess.run(
            ["git", "-C", str(repo), "rev-parse", f"r{i}"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip(): i
        for i in range(len(revisions))
    }
    case = tmp_path / "precomputed" / "fake"
    case.mkdir(parents=True)
    (case / "parameters.par").write_text("&output\n/\n")

    def bisect(work, fixups):
        arguments = [
            "--source",
            str(repo),
            "--case",
            str(case),
            "--work",
            str(work),
            "--range",
            "r0..main",
            "--stride",
            "1",
            "--launcher",
            "",
            "--threads",
            "1",
            "--jobs",
            "1",
            "--no-history",
        ]
        assert bh.main(arguments + [f"--fixup={fixup}" for fixup in fixups]) == 0
        return json.loads((work / "changes.json").read_text())

    # without a fixup, the change lies among the two commits that do not build
    plain = bisect(tmp_path / "plain", [])
    assert plain["changes"] == []
    assert [tags[sha] for sha in plain["unresolved"][0]["failed"]] == [1, 2]
    assert plain["patched"] == {}

    # with one, every commit builds, and the commits it changed are named
    seen = tmp_path / "seen.txt"
    work = tmp_path / "fixed"
    fixed = bisect(
        work, ["sed -i s/yes/no/ broken.txt", f'echo "$BISECT_COMMIT" >> {seen}']
    )
    assert [tags[change["at"]] for change in fixed["changes"]] == [2]
    assert not fixed["changes"][0]["fixups_differ"]
    assert sorted(tags[sha] for sha in fixed["patched"]) == [1, 2]
    assert all(paths == ["broken.txt"] for paths in fixed["patched"].values())
    assert set(seen.read_text().split()) == set(tags)
    report = (work / "changes.md").read_text()
    assert "## Commits that the fixups changed" in report
    assert "fixups differ" not in report
    log = (work / "logs" / f"{next(iter(fixed['patched']))[:9]}.build.log").read_text()
    assert "fixups changed: broken.txt" in log
