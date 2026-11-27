#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Run the verification suite of an installed SeisSol, without its build.

A build describes its verification tests to suite.json when it is installed
(see CMakeLists.txt): what CTest knows of each of them (ctest
--show-only=json-v1), and where the build had the suite, its own files, the
executable and the MPI launcher. This script runs the same tests with the
suite that the installation holds, the executable given, and a working
directory of its own, in the order that their fixtures demand, as CTest would:
a test whose setup failed is not run, a test that ends with its
SKIP_RETURN_CODE is skipped, and the cleanup of a fixture runs after the tests
that require it. The output of each test goes to logs/ of the working
directory, and that of a failed one to the terminal as well.

For example, on a machine with a GPU, against the installation of a build for
it, leaving out the spectral comparisons, which need no device, and the
snapshots, which hold the references of a CPU build:

    run_suite.py --binary /opt/seissol/bin/seissol --work /tmp/verification \\
                 --exclude-label 'spectral|snapshot' -j 4
"""

import argparse
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path


def load(path):
    """The origin of the suite, and its tests, from suite.json."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    tests = []
    for entry in data["ctest"]["tests"]:
        properties = {
            item["name"]: item["value"] for item in entry.get("properties", [])
        }
        tests.append({"name": entry["name"], "command": entry["command"], **properties})
    return data["origin"], tests


def relocation(origin, suite, work, binary, python, mpiexec):
    """The replacements that move a test from the build to this machine.

    The executable, the interpreter and the MPI launcher are replaced as a
    whole, also as the value of an option (--mpiexec=...); the source tree is
    not installed, so the commit a test would record is unknown. The paths of
    the suite, of its references and of the build directory of the suite are
    replaced by their prefixes, the longest first.
    """
    whole = {
        origin["binary"]: str(binary),
        origin["python"]: str(python),
        origin["project"]: "",
    }
    if origin.get("mpiexec"):
        whole[origin["mpiexec"]] = str(mpiexec)
    prefixes = sorted(
        [
            (origin["source"], str(suite)),
            (origin["references"], str(suite / "references")),
            (origin["build"], str(work)),
        ],
        key=lambda pair: -len(pair[0]),
    )

    def move(value):
        if isinstance(value, list):
            return [move(item) for item in value]
        if not isinstance(value, str):
            return value
        if value in whole:
            return whole[value]
        option, equals, argument = value.partition("=")
        if equals and argument in whole:
            return option + equals + whole[argument]
        for old, new in prefixes:
            value = value.replace(old, new)
        return value

    return move


def select(tests, args):
    """The tests the options ask for, with the setups of their fixtures."""

    def wanted(test):
        labels = test.get("LABELS", [])
        if args.regex and not re.search(args.regex, test["name"]):
            return False
        if args.exclude and re.search(args.exclude, test["name"]):
            return False
        if args.label and not any(re.search(args.label, label) for label in labels):
            return False
        if args.exclude_label and any(
            re.search(args.exclude_label, label) for label in labels
        ):
            return False
        return True

    chosen = {test["name"] for test in tests if wanted(test)}
    named = {test["name"]: test for test in tests}
    providers = {}
    for test in tests:
        for fixture in test.get("FIXTURES_SETUP", []) + test.get(
            "FIXTURES_CLEANUP", []
        ):
            providers.setdefault(fixture, []).append(test["name"])
    # as CTest does, the setups and cleanups of the fixtures of a chosen test
    # are chosen as well, whatever their names and labels
    pending = list(chosen)
    while pending:
        for fixture in named[pending.pop()].get("FIXTURES_REQUIRED", []):
            for provider in providers.get(fixture, []):
                if provider not in chosen:
                    chosen.add(provider)
                    pending.append(provider)
    return [test for test in tests if test["name"] in chosen]


def run(tests, move, jobs, timeout_factor, logs, log):
    """Run the tests, `jobs` processors at a time, in the order of their
    fixtures, each writing its output to a file of its own in `logs`; returns
    the outcome per test."""
    setups, cleanups = {}, {}
    for test in tests:
        for fixture in test.get("FIXTURES_SETUP", []):
            setups.setdefault(fixture, set()).add(test["name"])
        for fixture in test.get("FIXTURES_CLEANUP", []):
            cleanups.setdefault(fixture, set()).add(test["name"])
    requiring = {}
    for test in tests:
        for fixture in test.get("FIXTURES_REQUIRED", []):
            requiring.setdefault(fixture, set()).add(test["name"])

    def blockers(test):
        """The tests that have to end before this one starts."""
        names = set(test.get("DEPENDS", []))
        for fixture in test.get("FIXTURES_REQUIRED", []):
            names |= setups.get(fixture, set())
        for fixture in test.get("FIXTURES_CLEANUP", []):
            names |= requiring.get(fixture, set())
        return names - {test["name"]}

    known = {test["name"] for test in tests}
    waiting = {test["name"]: test for test in tests}
    outcome, running, locks = {}, {}, set()
    while waiting or running:
        started = False
        for name, test in list(waiting.items()):
            before = blockers(test) & known
            if not before <= set(outcome):
                continue
            failed = [
                other
                for fixture in test.get("FIXTURES_REQUIRED", [])
                for other in setups.get(fixture, set())
                if outcome.get(other) in ("failed", "not run")
            ]
            if failed:
                outcome[name] = "not run"
                log(f"not run  {name} (failed setup {failed[0]})")
                del waiting[name]
                started = True
                continue
            processors = max(1, int(test.get("PROCESSORS", 1)))
            used = sum(item["processors"] for item in running.values())
            held = set(test.get("RESOURCE_LOCK", []))
            if (running and used + processors > jobs) or held & locks:
                continue
            directory = move(test.get("WORKING_DIRECTORY", os.getcwd()))
            Path(directory).mkdir(parents=True, exist_ok=True)
            environment = dict(os.environ)
            for setting in move(test.get("ENVIRONMENT", [])):
                key, _, value = setting.partition("=")
                environment[key] = value
            output = open(logs / (re.sub(r"[^\w.-]+", "_", name) + ".log"), "w+")
            # a group of its own, so that a timeout ends the launcher and the ranks as well
            process = subprocess.Popen(
                move(test["command"]),
                cwd=directory,
                env=environment,
                stdout=output,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
            running[name] = {
                "process": process,
                "output": output,
                "processors": processors,
                "locks": held,
                "start": time.monotonic(),
                "timeout": float(test.get("TIMEOUT", 1500.0)) * timeout_factor,
                "skip": test.get("SKIP_RETURN_CODE"),
            }
            locks |= held
            del waiting[name]
            started = True
        for name, item in list(running.items()):
            process = item["process"]
            elapsed = time.monotonic() - item["start"]
            if process.poll() is None:
                if elapsed <= item["timeout"]:
                    continue
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
                state = "failed"
                note = f" (timeout after {elapsed:.0f} s)"
            else:
                state = (
                    "skipped"
                    if item["skip"] is not None
                    and process.returncode == int(item["skip"])
                    else "passed" if process.returncode == 0 else "failed"
                )
                note = f" ({elapsed:.1f} s)"
            item["output"].seek(0)
            output = item["output"].read()
            item["output"].close()
            outcome[name] = state
            locks -= item["locks"]
            del running[name]
            log(f"{state:8s} {name}{note}")
            if state == "failed":
                log(output.rstrip())
        if not started:
            time.sleep(0.05)
    return outcome


def main(argv=None):
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--suite", default=here / "suite.json", type=Path)
    parser.add_argument(
        "--binary", required=True, type=Path, help="the SeisSol executable"
    )
    parser.add_argument(
        "--work", required=True, type=Path, help="where the tests write"
    )
    parser.add_argument(
        "--mpiexec",
        help="the MPI launcher; by default the one of the build, if it is there, else mpiexec",
    )
    parser.add_argument("-R", "--regex", help="only tests whose name matches")
    parser.add_argument("-E", "--exclude", help="no tests whose name matches")
    parser.add_argument("-L", "--label", help="only tests with a label that matches")
    parser.add_argument(
        "-LE", "--exclude-label", help="no tests with a label that matches"
    )
    parser.add_argument("-j", "--jobs", type=int, default=1, help="processors to fill")
    parser.add_argument("--timeout-factor", type=float, default=1.0)
    parser.add_argument("--list", action="store_true", help="name the tests, run none")
    args = parser.parse_args(argv)

    origin, tests = load(args.suite)
    work = args.work.resolve()
    mpiexec = args.mpiexec
    if mpiexec is None:
        built = origin.get("mpiexec")
        mpiexec = (
            built
            if built and Path(built).exists()
            else shutil.which("mpiexec") or built
        )
    move = relocation(
        origin, here, work, args.binary.resolve(), sys.executable, mpiexec
    )
    tests = select(tests, args)
    if args.list:
        for test in tests:
            print(test["name"])
        return 0

    # what the configuration of the build left in its directory of the suite: the
    # capabilities of each configuration, the cases it declares, its environment
    for relative in origin["files"]:
        target = work / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(here / "build" / relative, target)

    logs = work / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    outcome = run(
        tests,
        move,
        args.jobs,
        args.timeout_factor,
        logs,
        lambda line: print(line, flush=True),
    )
    counts = {
        state: list(outcome.values()).count(state) for state in set(outcome.values())
    }
    print(", ".join(f"{count} {state}" for state, count in sorted(counts.items())))
    return 1 if counts.get("failed") or counts.get("not run") else 0


if __name__ == "__main__":
    sys.exit(main())
