# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""Tests for scripts/validate/components.py"""

import sys
from pathlib import Path

import pytest

SEISSOL_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(SEISSOL_ROOT / "scripts" / "validate"))

import components  # noqa: E402


class TestComponentGroup:
    @pytest.mark.parametrize(
        "name, group",
        [
            ("v1", "velocity"),
            ("s_xy", "stress"),
            ("u3", "displacement"),
            ("P_n", "traction"),
            ("Sls", "slip"),
            # a fluid velocity is not a solid velocity with a suffix
            ("v1_f", "fluid velocity"),
            # the suffixes of the fused simulations in the three outputs
            ("v13", "velocity3"),
            ("SRs-4", "slip rate-4"),
            ("v2-1", "velocity-1"),
            # no component
            ("Vr", "Vr"),
            ("eta", "eta"),
            ("u_n", "u_n"),
        ],
    )
    def test_group(self, name, group):
        assert components.component_group(name) == group

    def test_the_simulations_of_a_fused_run_are_apart(self):
        assert components.component_group("v1-1") != components.component_group("v2-2")


class TestGroupScales:
    def test_largest_norm_of_each_group(self):
        scales = components.group_scales(
            {"v1": 4.0, "v2": 1e-12, "v3": 2.0, "u1": 0.5, "Vr": 3.0}
        )
        assert scales == {"velocity": 4.0, "displacement": 0.5, "Vr": 3.0}
