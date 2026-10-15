# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""verify.py runs a case in the executable of its configuration, or else in an executable of
several configurations, which it tells the configuration by its canonical name."""

import sys
from pathlib import Path

import pytest

SEISSOL_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(SEISSOL_ROOT / "scripts" / "validate"))

import verify  # noqa: E402


@pytest.mark.parametrize(
    "config, canonical",
    [
        ("elastic-o6-f64", "elastic-linearck-o6-f64-stroud"),
        ("acoustic-o6-f32", "acoustic-linearck-o6-f32-stroud"),
        ("anisotropic-o6-f64", "anisotropic-linearck-o6-f64-stroud"),
        ("poroelastic-o6-f32", "poroelastic-stp-o6-f32-stroud"),
        ("viscoelastic-3-o6-f64", "viscoelastic-linearckanelastic-m3-o6-f64-stroud"),
        ("viscoacoustic-3-o6-f32", "viscoacoustic-linearckanelastic-m3-o6-f32-stroud"),
        ("elastic-o3-f64-s8", "elastic-linearck-o3-f64-stroud-s8"),
        ("elastic-linearck-o4-f32-dunavant", "elastic-linearck-o4-f32-dunavant"),
    ],
)
def test_canonical_config(config, canonical):
    assert verify.canonical_config(config) == canonical


def test_unknown_config_is_rejected():
    with pytest.raises(ValueError):
        verify.canonical_config("plastic-o6-f64")


def executable(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\n")
    path.chmod(0o755)
    return path


def test_own_executable_needs_no_configuration(tmp_path):
    own = executable(tmp_path / "bin" / "seissol-elastic-o6-f64")
    executable(tmp_path / "bin" / "seissol")
    assert verify.resolve_binary(tmp_path, "elastic-o6-f64") == (own.absolute(), None)


def test_shared_executable_is_told_the_configuration(tmp_path):
    shared = executable(tmp_path / "bin" / "seissol")
    assert verify.resolve_binary(tmp_path, "viscoelastic-3-o6-f32") == (
        shared.absolute(),
        "viscoelastic-linearckanelastic-m3-o6-f32-stroud",
    )


def test_shared_executable_without_the_configuration_is_passed_over(tmp_path):
    # as SeisSol's help lists the configurations of the executable
    shared = tmp_path / "bin" / "seissol"
    shared.parent.mkdir(parents=True)
    shared.write_text(
        "#!/bin/sh\n"
        "echo elastic-linearck-o6-f64-stroud\n"
        "echo '  Material: elastic'\n"
        "echo acoustic-linearck-o6-f32-stroud\n"
        "exit 1\n"
    )
    shared.chmod(0o755)
    assert verify.resolve_binary(tmp_path, "elastic-o6-f64") == (
        shared.absolute(),
        "elastic-linearck-o6-f64-stroud",
    )
    with pytest.raises(FileNotFoundError):
        verify.resolve_binary(tmp_path, "viscoelastic-3-o6-f64")


def test_a_configuration_under_another_name_is_no_reason_to_skip(tmp_path):
    # the same material, order and precision, but with a quadrature rule verify.py does not
    # assume: renamed in SeisSol, every case would be skipped and the run pass unrun
    shared = tmp_path / "bin" / "seissol"
    shared.parent.mkdir(parents=True)
    shared.write_text(
        "#!/bin/sh\n"
        "echo elastic-linearck-o6-f64-dunavant\n"
        "echo acoustic-linearck-o6-f32-stroud\n"
        "exit 1\n"
    )
    shared.chmod(0o755)
    with pytest.raises(verify.ConfigurationMismatch):
        verify.resolve_binary(tmp_path, "elastic-o6-f64")
    with pytest.raises(FileNotFoundError):
        verify.resolve_binary(tmp_path, "elastic-o6-f32")


def test_no_executable(tmp_path):
    with pytest.raises(FileNotFoundError):
        verify.resolve_binary(tmp_path, "elastic-o6-f64")


def test_a_run_that_skips_every_case_fails(tmp_path):
    # as one given an artifact without the executable: each case is skipped for want of one, and
    # skipped cases pass
    precomputed = tmp_path / "precomputed"
    (precomputed / "tpv5").mkdir(parents=True)
    (precomputed / "cases.json").write_text(
        '{"tpv5/double": {"config": "elastic-o6-f64", "parameters": "tpv5/parameters.par", '
        '"output": {"prefix": "tpv", "directory": "tpv5/precomputed/double"}}}'
    )
    (tmp_path / "install" / "bin").mkdir(parents=True)
    argv = [
        "--precomputed-dir",
        str(precomputed),
        "--binary-dir",
        str(tmp_path / "install"),
        "--all",
        "--color",
        "never",
        "--step-summary",
        "never",
    ]
    assert verify.main(argv) == 1
