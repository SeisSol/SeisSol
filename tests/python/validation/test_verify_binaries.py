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


def test_no_executable(tmp_path):
    with pytest.raises(FileNotFoundError):
        verify.resolve_binary(tmp_path, "elastic-o6-f64")
