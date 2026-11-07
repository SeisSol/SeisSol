# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Builds with compile_variants.py, and the executables they install."""

import os
import time

import pytest
from strongscaling import build, config
from strongscaling.campaign import Layout


def test_the_build_command(campaign_file, tmp_path):
    variants = tmp_path / "variants.yaml"
    variants.write_text("version: 1\n")
    campaign = config.load_campaign(
        campaign_file(
            builds=[
                {
                    "name": "current",
                    "variants": "variants.yaml",
                    "preset": "quick",
                    "toolchain": "gcc",
                    "set": {"host_arch": "skx", "defines.ASAGI": "OFF"},
                }
            ]
        )
    )
    directory = Layout(campaign.directory)
    command = build.compile_command(campaign.builds[0], directory)
    assert command[1] == str(build.COMPILE_VARIANTS)
    assert command[command.index("--config") + 1] == str(variants)
    assert command[command.index("--install-prefix") + 1] == str(
        directory.installs / "current"
    )
    assert command[command.index("--preset") + 1] == "quick"
    assert command[command.index("--toolchain") + 1] == "gcc"
    sets = [command[i + 1] for i, word in enumerate(command) if word == "--set"]
    assert sets == ["host_arch=skx", "defines.ASAGI=OFF"]
    assert "--no-test" in command


def report(*variants):
    return {"variants": list(variants)}


def variant(tmp_path, name, status="ok"):
    binary = tmp_path / name
    binary.write_text("")
    return {
        "status": status,
        "binary_name": name,
        "binary_path": str(binary),
        "log": "log",
    }


def test_the_executable_of_a_build(tmp_path):
    one = config.Build(name="b")
    assert build.binary_of_report(one, report(variant(tmp_path, "seissol-a"))) == (
        tmp_path / "seissol-a"
    )
    two = report(variant(tmp_path, "seissol-a"), variant(tmp_path, "seissol-b"))
    with pytest.raises(build.BuildError, match="select"):
        build.binary_of_report(one, two)
    chosen = config.Build(name="b", select="seissol-b")
    assert build.binary_of_report(chosen, two) == tmp_path / "seissol-b"
    with pytest.raises(build.BuildError, match="failed"):
        build.binary_of_report(one, report(variant(tmp_path, "seissol-a", "failed")))


def test_fingerprints(tmp_path):
    folder = tmp_path / "install" / "bin"
    folder.mkdir(parents=True)
    binary = folder / "seissol"
    binary.write_text("one")
    first = build.fingerprint([folder])
    assert build.fingerprint([folder]) == first
    os.utime(binary, ns=(time.time_ns(), time.time_ns() + 10**9))
    assert build.fingerprint([folder]) != first


def test_given_executables(campaign_file, fakes):
    campaign = config.load_campaign(campaign_file())
    state_builds = {}
    directory = Layout(campaign.directory)
    binaries = build.build_all(campaign, directory, state_builds, log=lambda _: None)
    assert binaries == {"fake": fakes["seissol"]}
    assert state_builds["fake"]["given"] is True
    assert build.known_fingerprints(campaign, state_builds)["fake"] == (
        state_builds["fake"]["fingerprint"]
    )


def test_a_missing_executable(campaign_file, tmp_path):
    campaign = config.load_campaign(
        campaign_file(builds=[{"name": "gone", "binary": str(tmp_path / "nothere")}])
    )
    with pytest.raises(build.BuildError, match="does not exist"):
        build.build_all(campaign, Layout(campaign.directory), {}, log=lambda _: None)
