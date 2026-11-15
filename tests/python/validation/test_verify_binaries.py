# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""verify.py runs a case in the executable of its configuration, or else in an executable of
several configurations, which it tells the configuration by its name; and it compares a case with
dynamic rupture only with a reference of the same quadrature rule on the faults."""

import sys
from pathlib import Path

import pytest

SEISSOL_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(SEISSOL_ROOT / "scripts" / "validate"))

import verify  # noqa: E402


@pytest.mark.parametrize(
    "config",
    [
        "elastic-o6-f64",
        "acoustic-o6-f32",
        "anisotropic-o6-f64",
        "poroelastic-o6-f32",
        "viscoelastic-m3-o6-f64",
        "viscoacoustic-m3-o6-f32",
        "elastic-o3-f64-s8",
    ],
)
def test_canonical_config(config):
    # cases.json selects a configuration by the name SeisSol gives it
    assert verify.canonical_config(config) == config


@pytest.mark.parametrize(
    "config",
    [
        "plastic",
        "viscoelastic-3-o6-f64",
        "elastic-linearck-o6-f64-stroud",
        "elastic-o6-double",
    ],
)
def test_unknown_config_is_rejected(config):
    with pytest.raises(ValueError):
        verify.canonical_config(config)


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
    assert verify.resolve_binary(tmp_path, "viscoelastic-m3-o6-f32") == (
        shared.absolute(),
        "viscoelastic-m3-o6-f32",
    )


def help_script(path, lines):
    """An executable that prints `lines` as SeisSol prints its help."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "#!/bin/sh\n" + "".join(f"echo '{line}'\n" for line in lines) + "exit 1\n"
    )
    path.chmod(0o755)
    return path


def test_shared_executable_without_the_configuration_is_passed_over(tmp_path):
    shared = help_script(
        tmp_path / "bin" / "seissol",
        [
            "elastic-o6-f64",
            "  Material: elastic",
            "  Solver: linearck",
            "acoustic-o6-f32",
            "Dynamic rupture quadrature rule (all configurations): dunavant",
        ],
    )
    assert verify.resolve_binary(tmp_path, "elastic-o6-f64") == (
        shared.absolute(),
        "elastic-o6-f64",
    )
    with pytest.raises(FileNotFoundError):
        verify.resolve_binary(tmp_path, "viscoelastic-m3-o6-f64")


def test_the_quadrature_rule_is_read_from_the_help(tmp_path):
    binary = help_script(
        tmp_path / "bin" / "seissol",
        [
            "elastic-o6-f64",
            "Dynamic rupture quadrature rule (all configurations): stroud",
        ],
    )
    assert verify.built_dr_quad_rule(binary) == "stroud"
    silent = help_script(tmp_path / "other" / "seissol", ["elastic-o6-f64"])
    assert verify.built_dr_quad_rule(silent) is None


def rupture_case(tmp_path, rule):
    """A case with fault output whose reference states `rule`, if any."""
    precomputed = tmp_path / "precomputed"
    output = precomputed / "tpv5" / "precomputed" / "double"
    output.mkdir(parents=True)
    (output / "tpv-fault.xdmf").write_text("")
    if rule is not None:
        (precomputed / "references.json").write_text(
            '{"tpv5/double": {"drQuadRule": "' + rule + '"}}'
        )
    case = verify.CaseSpec(
        key="tpv5/double",
        config="elastic-o6-f64",
        params_path=precomputed / "tpv5" / "parameters.par",
        case_dir=precomputed / "tpv5",
        prefix="tpv",
        output_dir=output,
        reference_dir=output,
    )
    return precomputed, case


def test_a_reference_of_another_quadrature_rule_is_not_compared(tmp_path):
    precomputed, case = rupture_case(tmp_path, "stroud")
    binary = help_script(
        tmp_path / "bin" / "seissol",
        [
            "elastic-o6-f64",
            "Dynamic rupture quadrature rule (all configurations): dunavant",
        ],
    )
    with pytest.raises(verify.RuleMismatch):
        verify.check_dr_quad_rule(precomputed, case, binary)


def test_a_reference_of_the_same_or_no_stated_rule_is_compared(tmp_path):
    binary = help_script(
        tmp_path / "bin" / "seissol",
        [
            "elastic-o6-f64",
            "Dynamic rupture quadrature rule (all configurations): dunavant",
        ],
    )
    precomputed, case = rupture_case(tmp_path / "same", "dunavant")
    verify.check_dr_quad_rule(precomputed, case, binary)
    precomputed, case = rupture_case(tmp_path / "unstated", None)
    verify.check_dr_quad_rule(precomputed, case, binary)


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
