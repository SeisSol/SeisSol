# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The executables of a campaign: built with compile_variants.py, or given.

A build of a variant-set file is built anew whenever the campaign submits
runs, so that a series measures the source as it is; CMake rebuilds only what
changed. The source it was built from is recorded with it.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess  # nosec B404
import sys
from pathlib import Path
from typing import Any

from .campaign import Layout as Directory
from .campaign import read_json, utc_now
from .config import Build, Campaign

SOURCE_ROOT = Path(__file__).resolve().parents[3]
COMPILE_VARIANTS = SOURCE_ROOT / "scripts" / "compile" / "compile_variants.py"


class BuildError(Exception):
    """An executable that could not be built or found."""


def source_revision(root: Path = SOURCE_ROOT) -> dict[str, Any]:
    """The commit of the source tree, and whether it has changes on top."""
    if shutil.which("git") is None:
        return {"commit": None, "dirty": None}
    commit = subprocess.run(  # nosec B603 B607
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        check=False,
        capture_output=True,
        text=True,
    )
    status = subprocess.run(  # nosec B603 B607
        ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"],
        check=False,
        capture_output=True,
        text=True,
    )
    return {
        "commit": commit.stdout.strip() or None,
        "dirty": bool(status.stdout.strip()) if status.returncode == 0 else None,
    }


def fingerprint(paths: list[Path]) -> str:
    """A hash of the names, sizes and times of the files in @p paths (files or directories).

    An installation that CMake brought up to date without a change keeps its files as they
    were, and its fingerprint.
    """
    entries = []
    for path in paths:
        files = (
            [path]
            if path.is_file()
            else sorted(p for p in path.rglob("*") if p.is_file())
        )
        for file in files:
            status = file.stat()
            entries.append(f"{file}|{status.st_size}|{status.st_mtime_ns}")
    return hashlib.sha256("\n".join(entries).encode("utf-8")).hexdigest()[:16]


def compile_command(build: Build, directory: Directory) -> list[str]:
    build_root = directory.builds / build.name
    command = [
        sys.executable,
        str(COMPILE_VARIANTS),
        "--config",
        str(build.variants_file),
        "--build-root",
        str(build_root),
        "--install-prefix",
        str(directory.installs / build.name),
        "--report",
        str(build_root / "build-report.json"),
        "--no-test",
    ]
    if build.preset:
        command += ["--preset", build.preset]
    if build.toolchain:
        command += ["--toolchain", build.toolchain]
    for override in build.overrides:
        command += ["--set", override]
    return command


def binary_of_report(build: Build, report: dict[str, Any]) -> Path:
    """The executable a build report of compile_variants.py names for @p build."""
    variants = report.get("variants", [])
    failed = [v for v in variants if v.get("status") != "ok"]
    if failed:
        logs = ", ".join(str(v.get("log")) for v in failed)
        raise BuildError(
            f"build {build.name}: {len(failed)} variant(s) failed; see {logs}"
        )
    candidates = variants
    if build.select:
        candidates = [v for v in variants if v.get("binary_name") == build.select]
    if len(candidates) != 1:
        names = ", ".join(str(v.get("binary_name")) for v in variants)
        raise BuildError(
            f"build {build.name}: the variant set built {names}; name the executable "
            "to take with select"
        )
    path = Path(candidates[0]["binary_path"])
    if not path.exists():
        raise BuildError(f"build {build.name}: {path} was not installed")
    return path


def build_all(
    campaign: Campaign, directory: Directory, state_builds: dict[str, Any], log=print
) -> dict[str, Path]:
    """Builds what the campaign builds, and finds what it is given; the executables by build."""
    binaries: dict[str, Path] = {}
    revision = source_revision()
    for build in campaign.builds:
        if build.binary is not None:
            if not build.binary.exists():
                raise BuildError(f"build {build.name}: {build.binary} does not exist")
            if not os.access(build.binary, os.X_OK):
                raise BuildError(
                    f"build {build.name}: {build.binary} is not executable"
                )
            binaries[build.name] = build.binary
            state_builds[build.name] = {
                "binary": str(build.binary),
                "given": True,
                "fingerprint": fingerprint([build.binary]),
            }
            continue
        if not build.variants_file.exists():
            raise BuildError(
                f"build {build.name}: {build.variants_file} does not exist"
            )
        command = compile_command(build, directory)
        (directory.builds / build.name).mkdir(parents=True, exist_ok=True)
        log_path = directory.builds / build.name / "compile-variants.log"
        log(f"building {build.name} (log: {log_path})")
        with log_path.open("w", encoding="utf-8") as stream:
            result = subprocess.run(  # nosec B603
                command, check=False, stdout=stream, stderr=subprocess.STDOUT
            )
        report_path = directory.builds / build.name / "build-report.json"
        if result.returncode != 0 and not report_path.exists():
            raise BuildError(f"build {build.name} failed; see {log_path}")
        report = read_json(report_path)
        binaries[build.name] = binary_of_report(build, report)
        prefix = directory.installs / build.name
        state_builds[build.name] = {
            "binary": str(binaries[build.name]),
            "given": False,
            "built": utc_now(),
            "source": revision,
            "report": str(report_path),
            "log": str(log_path),
            "fingerprint": fingerprint(
                [
                    p
                    for p in (prefix / "bin", prefix / "lib", prefix / "lib64")
                    if p.exists()
                ]
            ),
        }
    return binaries


def known_fingerprints(
    campaign: Campaign, state_builds: dict[str, Any]
) -> dict[str, str]:
    """The fingerprints of the executables: of those given as they are, of the others as built."""
    fingerprints = {}
    for build in campaign.builds:
        if build.binary is not None and build.binary.exists():
            fingerprints[build.name] = fingerprint([build.binary])
        elif state_builds.get(build.name, {}).get("fingerprint"):
            fingerprints[build.name] = state_builds[build.name]["fingerprint"]
    return fingerprints


def known_binaries(campaign: Campaign, state_builds: dict[str, Any]) -> dict[str, Path]:
    """The executables of the builds as given or built before, without building."""
    binaries = {}
    for build in campaign.builds:
        if build.binary is not None:
            binaries[build.name] = build.binary
        elif build.name in state_builds:
            binaries[build.name] = Path(state_builds[build.name]["binary"])
    return binaries
