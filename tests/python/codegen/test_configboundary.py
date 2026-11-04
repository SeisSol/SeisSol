# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The kernels that the configurations of a build get for their faces towards
other configurations (kernels.configboundary.plans)."""

from types import SimpleNamespace

import kernels  # noqa: F401 — establish sys.path via conftest
from kernels.configboundary import Plan, plans


def config(equations, order, simulations=1):
    return SimpleNamespace(
        equations=equations, order=order, multipleSimulations=simulations
    )


def test_a_configuration_alone_converts_nothing():
    assert plans([config("elastic", 6)]) == [Plan(None, False, None)]


def test_anisotropic_is_of_the_family_of_the_isotropic_solids():
    result = plans(
        [config("elastic", 6), config("viscoelastic", 5), config("anisotropic", 4)]
    )
    assert result == [Plan(6, True, None)] * 3


def test_an_anisotropic_solid_does_not_couple_to_a_fluid():
    elastic, anisotropic, acoustic = plans(
        [config("elastic", 5), config("anisotropic", 6), config("acoustic", 4)]
    )
    # the canonical form of the solids has the order of the anisotropic one, which the fluid reads
    assert elastic == Plan(6, True, 4)
    assert anisotropic == Plan(6, True, None)
    assert acoustic == Plan(4, False, 6)


def test_a_fluid_with_anisotropic_solids_only_couples_to_nothing():
    anisotropic, acoustic = plans([config("anisotropic", 6), config("acoustic", 4)])
    assert anisotropic == Plan(None, False, None)
    assert acoustic == Plan(None, False, None)


def test_families_are_apart_by_their_fused_simulations():
    result = plans([config("elastic", 3, 8), config("anisotropic", 3)])
    assert result == [Plan(None, False, None)] * 2
