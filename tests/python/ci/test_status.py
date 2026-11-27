# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""scripts/ci/status.py tells why the jobs of the CI have failed, from what they upload: which
artifacts belong to which job, and what the build reports, the results of the tests (JUnit of ctest
and pytest), the findings of clang-tidy and the logs of the builds in them say."""

import importlib.util
import io
import json
import sys
import zipfile
from pathlib import Path

SEISSOL_ROOT = Path(__file__).resolve().parents[3]
_spec = importlib.util.spec_from_file_location(
    "ci_status", SEISSOL_ROOT / "scripts" / "ci" / "status.py"
)
status = importlib.util.module_from_spec(_spec)
sys.modules["ci_status"] = status
_spec.loader.exec_module(status)


def test_the_artifacts_of_a_build_are_those_with_the_parts_of_its_matrix():
    job = "seissol-build-test<gcc-13, hsw, Release>"
    assert status.artifact_matches(job, "build-report-Release-hsw-gcc-13")
    assert status.artifact_matches(job, "logs-Release-hsw-gcc-13")
    # the installed SeisSol does not tell why, and another compiler is another job
    assert not status.artifact_matches(job, "seissol-Release-hsw-gcc-13")
    assert not status.artifact_matches(job, "logs-Release-hsw-clang-19")


def test_the_artifacts_of_the_other_jobs():
    assert status.artifact_matches("clang-tidy (2/4)", "tidy-cpu-2")
    assert not status.artifact_matches("clang-tidy (2/4)", "tidy-gpu-2")
    assert status.artifact_matches("clang-tidy-gpu (3/4)", "tidy-gpu-3")
    assert status.artifact_matches("pytest-3.12", "pytest-3.12")
    assert status.artifact_matches("sanitizer<thread>", "sanitizer-thread-logs")
    assert status.artifact_matches("emulated<rvv128>", "emulated-rvv128-logs")
    assert status.artifact_matches("build-windows", "windows-logs")
    assert status.artifact_matches(
        "seissol-run-test<gcc-13, hsw, Release><r12×t1>", "report-gcc-13-r12t1"
    )
    assert not status.artifact_matches("flake8", "flake8-logs")


CTEST_JUNIT = """<?xml version="1.0" encoding="UTF-8"?>
<testsuite name="ci" tests="3" failures="1">
  <testcase name="io/IO/BinaryFile: a write" classname="io" status="run">
    <system-out>passed</system-out>
  </testcase>
  <testcase name="verification.null-mat48&lt;elastic-o3-f64-s8&gt;" classname="verification"
      status="fail">
    <failure message="Failed"/>
    <system-out>info 0 : Check the mesh for geometric errors.
the receiver 3 deviates by 1e-3 from zero
done</system-out>
  </testcase>
  <testcase name="skipped one" classname="x" status="notrun"/>
</testsuite>
"""

PYTEST_JUNIT = """<?xml version="1.0" encoding="utf-8"?>
<testsuites><testsuite name="pytest" errors="0" failures="1" tests="2">
  <testcase classname="codegen.test_configboundary" name="test_alone" time="0.01">
    <failure message="AttributeError: no attribute 'precision'">Traceback
E   AttributeError: 'SimpleNamespace' object has no attribute 'precision'</failure>
  </testcase>
  <testcase classname="codegen.test_configboundary" name="test_other" time="0.01"/>
</testsuite></testsuites>
"""


def test_the_failed_tests_of_ctest_with_what_went_wrong():
    ((name, lines),) = status.junit_failures(CTEST_JUNIT, 5)
    assert name.startswith("verification.null-mat48<elastic-o3-f64-s8>")
    # the information of the log of SeisSol is no error, the deviation is
    assert lines == ["the receiver 3 deviates by 1e-3 from zero"]


def test_the_failed_tests_of_pytest():
    ((name, lines),) = status.junit_failures(PYTEST_JUNIT, 5)
    assert name == "codegen.test_configboundary::test_alone"
    assert any("has no attribute 'precision'" in line for line in lines)


def test_without_a_line_of_an_error_the_end_of_the_output():
    assert status.excerpt("one\ntwo\nthree\nFailed", 2) == ["three", "Failed"]


def test_the_failed_variants_of_a_build_report():
    report = {
        "variants": [
            {"id": "ok", "status": "ok"},
            {
                "id": "ci-cpu-hsw-gcc-13",
                "status": "failed",
                "phase_failed": "test",
                "testing": {"total": 4473, "failed": 2, "failed_tests": ["a", "b"]},
            },
        ]
    }
    assert status.build_report_failures(json.dumps(report)) == [
        ("ci-cpu-hsw-gcc-13 failed in test: 2 of 4473 tests", ["a", "b"])
    ]


def test_the_errors_of_a_build_log():
    log = "[1/2] Building x.o\nsrc/x.cpp:3:1: error: expected ';'\nFAILED: x.o\nninja: build stopped"
    assert status.build_log_errors(log, 2) == [
        "src/x.cpp:3:1: error: expected ';'",
        "FAILED: x.o",
    ]


def zip_of(files):
    data = io.BytesIO()
    with zipfile.ZipFile(data, "w") as archive:
        for name, content in files.items():
            archive.writestr(name, content)
    return data.getvalue()


def test_an_artifact_with_the_logs_of_a_build():
    found = status.read_artifact(
        "logs-Release-hsw-gcc-13",
        zip_of(
            {
                "ci-cpu-hsw-gcc-13/ctest-results.xml": CTEST_JUNIT,
                "ci-cpu-hsw-gcc-13/build.log": "all fine",
                "ci-cpu-hsw-gcc-13/Testing/Temporary/LastTestsFailed.log": "7:x\n",
            }
        ),
        3,
    )
    assert found[0] == "ci-cpu-hsw-gcc-13/ctest-results.xml: 1 failed"
    assert found[1].startswith("  - verification.null-mat48")
    # the list of the failed tests of ctest only where there are no results as JUnit
    assert not any("LastTestsFailed" in line for line in found)


def test_an_artifact_with_the_findings_of_clang_tidy():
    found = status.read_artifact(
        "tidy-cpu-0",
        zip_of(
            {
                "tidy-findings.txt": "a.h:1:2: error: redundant 'typename' [x]\n",
                "tidy-output.txt": "lots",
            }
        ),
        3,
    )
    assert found == [
        "1 findings of clang-tidy",
        "  a.h:1:2: error: redundant 'typename' [x]",
    ]


class FakeGitHub:
    def __init__(self, answers):
        self.answers = answers

    def get(self, path):
        return self.answers[path]


def test_a_pull_request_a_branch_and_a_commit():
    github = FakeGitHub(
        {
            "pulls/1615": {"head": {"sha": "abc123", "ref": "davschneller/simcopy"}},
            "commits/davschneller%2Fsimcopy": {"sha": "abc123"},
            "commits/abc1": {"sha": "abc123"},
        }
    )
    assert status.resolve(github, "1615") == (
        "abc123",
        "pull request #1615 (davschneller/simcopy)",
    )
    assert status.resolve(github, "davschneller/simcopy") == (
        "abc123",
        "branch davschneller/simcopy",
    )
    assert status.resolve(github, "abc1") == ("abc123", "commit")


def test_the_log_of_a_job_up_to_its_error():
    log = "\n".join(
        ["setup", "build 1", "build 2", "##[error]Process completed with exit code 1."]
    )
    assert status.log_excerpt(log, 2) == [
        "build 1",
        "build 2",
        "##[error]Process completed with exit code 1.",
    ]
    assert (
        status.run_of({"html_url": "https://github.com/o/r/actions/runs/42/job/7"})
        == 42
    )
