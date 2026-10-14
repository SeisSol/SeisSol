# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""compile_variants.py builds variants of one configuration, which equations, order and
precision describe, and of several, which `configs` names as CMake's CONFIGS takes them.
"""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

SEISSOL_ROOT = Path(__file__).resolve().parents[3]
_spec = importlib.util.spec_from_file_location(
    "compile_variants", SEISSOL_ROOT / "scripts" / "compile" / "compile_variants.py"
)
cv = importlib.util.module_from_spec(_spec)
# the dataclasses of the module look themselves up there
sys.modules["compile_variants"] = cv
_spec.loader.exec_module(cv)


def variants(tmp_path, preset):
    path = tmp_path / "variants.json"
    path.write_text(json.dumps({"version": 1, "presets": {"p": preset}}))
    return cv.build_variants(cv.load_config(path), "p")


def defines(variant):
    return [arg for arg in variant.cmake_configure_cmd() if arg.startswith("-D")]


def test_one_configuration(tmp_path):
    (variant,) = variants(
        tmp_path,
        {"variants": [{"equations": "elastic", "order": 6, "precision": "double"}]},
    )
    assert variant.binary_name == "seissol-elastic-o6-f64"
    assert "-DEQUATIONS=elastic" in defines(variant)
    assert not any(arg.startswith("-DCONFIGS=") for arg in defines(variant))


def test_several_configurations(tmp_path):
    (variant,) = variants(tmp_path, {"variants": [{"configs": "ci-cpu"}]})
    # the executable of several configurations is named after none of them
    assert variant.binary_name == "seissol"
    assert variant.slug.startswith("ci-cpu")
    assert "-DCONFIGS=ci-cpu" in defines(variant)
    assert not any(arg.startswith("-DEQUATIONS=") for arg in defines(variant))


def test_several_configurations_as_a_list(tmp_path):
    (variant,) = variants(
        tmp_path,
        {
            "variants": [
                {
                    "configs": [
                        "elastic-linearck-o4-f64-stroud",
                        "elastic-linearck-o4-f32-stroud",
                    ],
                    "custom_suffix": "o4",
                }
            ]
        },
    )
    assert variant.binary_name == "seissol-o4"
    assert (
        "-DCONFIGS=elastic-linearck-o4-f64-stroud;elastic-linearck-o4-f32-stroud"
        in defines(variant)
    )


def test_suffix_keeps_executables_of_several_configurations_apart(tmp_path):
    # as the GPU CI installs its two executables into one prefix
    first, second = variants(
        tmp_path,
        {
            "variants": [
                {"configs": "ci-gpu", "device_codegen": "tensorforge"},
                {
                    "configs": "elastic-linearck-o6-f32-stroud",
                    "device_codegen": "gemmforge-chainforge",
                    "custom_suffix": "gemmforge",
                },
            ]
        },
    )
    assert (first.binary_name, second.binary_name) == ("seissol", "seissol-gemmforge")
    assert "-DCUSTOM_BINARY_SUFFIX=gemmforge" in defines(second)


def test_configurations_given_twice_are_refused(tmp_path):
    with pytest.raises(SystemExit):
        variants(
            tmp_path,
            {"variants": [{"configs": "ci-cpu", "equations": "elastic"}]},
        )


def test_a_variant_needs_its_configurations(tmp_path):
    with pytest.raises(SystemExit):
        variants(tmp_path, {"variants": [{"order": 6, "precision": "double"}]})


def test_ci_presets():
    pytest.importorskip("yaml")
    config = cv.load_config(
        SEISSOL_ROOT / "scripts" / "compile" / "compile-variants.yaml"
    )
    (cpu,) = cv.build_variants(config, "ci-cpu")
    assert cpu.binary_name == "seissol"
    assert "-DCONFIGS=ci-cpu" in defines(cpu)
    gpu = cv.build_variants(config, "ci-gpu")
    assert [v.binary_name for v in gpu] == ["seissol", "seissol-gemmforge"]
