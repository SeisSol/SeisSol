# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Tests for pure math and constants in codegen/kernels.

These cover:
 - closed-form basis-function / quadrature-point counts
 - DoF-index sparsity patterns (extractVelocities, extractTractions)
 - Sparsity patterns for transformation / flux / godunov matrices

All of these are SILENT-FAILURE domains: a wrong index or off-by-one
produces a compiling, running simulation whose output is wrong.

We bypass the ADERDGBase constructor (which parses XML matrix files)
by constructing concrete subclasses with object.__new__ + direct
attribute injection, so we can test the static methods and the
spp helpers without needing a full yateto environment.
"""

import importlib.util  # noqa: F401 ; yateto expects importlib.util pre-loaded
from math import comb

import numpy as np
import pytest


# Lazy imports so the module loads even if yateto is broken
def _mk(cls, **kwargs):
    """Construct a codegen equation class bypassing its heavy __init__.

    The __init__ parses XML/JSON matrix files; we only want to test the
    pure-computation methods. Replace __init__ with a no-op that records
    the user-provided attributes, then restore it.
    """
    instance = object.__new__(cls)
    for k, v in kwargs.items():
        setattr(instance, k, v)
    return instance


# =============================================================================
# Basis-function counts — closed-form ADER-DG formulas
# =============================================================================


class TestBasisFunctionCounts:
    """2D/3D basis function counts and 3D quadrature point counts follow
    closed-form formulas. These are called everywhere the codegen allocates
    tensors, so off-by-one here means wrong memory allocation everywhere.

    In these tests, `order` means the convergence order in all cases.
    """

    @pytest.fixture
    def elastic(self):
        from kernels.equations.elastic import ElasticADERDG

        # Orders in the actual CMake build matrix: 2..8
        return lambda order: _mk(ElasticADERDG, order=order)

    @pytest.mark.parametrize(
        "order,expected",
        [
            (1, 1),  # triangular numbers T_{order}
            (2, 3),
            (3, 6),
            (4, 10),
            (5, 15),
            (6, 21),
            (7, 28),
            (8, 36),
        ],
    )
    def test_numberOf2DBasisFunctions(self, elastic, order, expected):
        adg = elastic(order)
        assert adg.numberOf2DBasisFunctions() == expected

    @pytest.mark.parametrize(
        "order,expected",
        [
            (1, 1),  # tetrahedral numbers Te_{order}
            (2, 4),
            (3, 10),
            (4, 20),
            (5, 35),
            (6, 56),
            (7, 84),
            (8, 120),
        ],
    )
    def test_numberOf3DBasisFunctions(self, elastic, order, expected):
        adg = elastic(order)
        assert adg.numberOf3DBasisFunctions() == expected

    @pytest.mark.parametrize("order", [1, 2, 3, 4, 5, 6, 7, 8])
    def test_numberOf3DQuadraturePoints_is_cubic(self, elastic, order):
        adg = elastic(order)
        assert adg.numberOf3DQuadraturePoints() == (order + 1) ** 3

    @pytest.mark.parametrize("order", [1, 3, 5, 7])
    def test_2d_count_matches_pascal_identity(self, elastic, order):
        """T_n = C(n+1, 2) — another way to state the formula."""
        adg = elastic(order)
        assert adg.numberOf2DBasisFunctions() == comb(order + 1, 2)

    @pytest.mark.parametrize("order", [1, 3, 5, 7])
    def test_3d_count_matches_pascal_identity(self, elastic, order):
        """Te_n = C(n+2, 3)."""
        adg = elastic(order)
        assert adg.numberOf3DBasisFunctions() == comb(order + 2, 3)

    def test_counts_are_monotonically_increasing(self, elastic):
        for order in range(1, 8):
            assert (
                elastic(order + 1).numberOf2DBasisFunctions()
                > elastic(order).numberOf2DBasisFunctions()
            )
            assert (
                elastic(order + 1).numberOf3DBasisFunctions()
                > elastic(order).numberOf3DBasisFunctions()
            )

    def test_counts_are_integers(self, elastic):
        # Integer division // must not yield floats
        for order in [2, 3, 4, 5, 6]:
            adg = elastic(order)
            assert isinstance(adg.numberOf2DBasisFunctions(), int)
            assert isinstance(adg.numberOf3DBasisFunctions(), int)


# =============================================================================
# DoF-index sparsity patterns — the silent-failure domain
# =============================================================================


class TestElasticDofIndices:
    """ElasticADERDG encodes:
      - 9 quantities: 6 stress components + 3 velocity components
      - Velocities at indices 6, 7, 8
      - Tractions at indices 0, 3, 5 (s_xx, s_yy, s_zz in Voigt-like order)

    If a refactor permutes these, simulations silently output nonsense.
    """

    @pytest.fixture
    def adg(self):
        from kernels.equations.elastic import ElasticADERDG

        return _mk(ElasticADERDG, order=4)

    def test_extractVelocities_shape(self, adg):
        v = adg.extractVelocities()
        assert v.shape == (3, 9)  # 3 vel components × 9 elastic quantities

    def test_extractVelocities_is_sparse(self, adg):
        v = adg.extractVelocities()
        assert np.count_nonzero(v) == 3  # exactly one 1 per row

    def test_extractVelocities_indices(self, adg):
        v = adg.extractVelocities()
        assert v[0, 6] == 1
        assert v[1, 7] == 1
        assert v[2, 8] == 1

    def test_extractTractions_shape(self, adg):
        t = adg.extractTractions()
        assert t.shape == (3, 9)

    def test_extractTractions_indices(self, adg):
        t = adg.extractTractions()
        assert t[0, 0] == 1
        assert t[1, 3] == 1
        assert t[2, 5] == 1
        assert np.count_nonzero(t) == 3

    def test_mapTo_is_transpose_of_extract(self, adg):
        """mapToVelocities/Tractions is defined as .T — invariant."""
        np.testing.assert_array_equal(adg.mapToVelocities(), adg.extractVelocities().T)
        np.testing.assert_array_equal(adg.mapToTractions(), adg.extractTractions().T)

    def test_velocity_and_traction_indices_disjoint(self, adg):
        """A DoF can't be both a velocity and a traction component."""
        v_indices = np.nonzero(adg.extractVelocities().any(axis=0))[0]
        t_indices = np.nonzero(adg.extractTractions().any(axis=0))[0]
        assert not set(v_indices) & set(
            t_indices
        ), "Velocity and traction index sets must be disjoint"


class TestAcousticDofIndices:
    """AcousticADERDG: 4 quantities = 1 (negative) pressure + 3 velocities.
    - Pressure (= single traction component) at index 0
    - Velocities at 1, 2, 3
    """

    @pytest.fixture
    def adg(self):
        from kernels.equations.acoustic import AcousticADERDG

        return _mk(AcousticADERDG, order=4)

    def test_numberOfQuantities(self, adg):
        assert adg.numberOfQuantities() == 4

    def test_extractVelocities(self, adg):
        v = adg.extractVelocities()
        assert v.shape == (3, 4)
        assert v[0, 1] == 1 and v[1, 2] == 1 and v[2, 3] == 1
        assert np.count_nonzero(v) == 3

    def test_extractTractions_only_one_component(self, adg):
        """Acoustic has a single scalar pressure, not three tractions.
        The pattern has 3 rows (for interface shape consistency) but only
        row 0 is non-zero.
        """
        t = adg.extractTractions()
        assert t.shape == (3, 4)
        assert t[0, 0] == 1
        assert np.count_nonzero(t) == 1  # only the pressure

    def test_pressure_not_in_velocity_indices(self, adg):
        v = adg.extractVelocities()
        t = adg.extractTractions()
        p_idx = np.nonzero(t[0])[0]
        v_idxs = set()
        for i in range(3):
            v_idxs.update(np.nonzero(v[i])[0].tolist())
        assert set(p_idx) & v_idxs == set()


class TestPoroelasticDofIndices:
    """PoroelasticADERDG: 13 quantities.
    - Solid velocities at 6, 7, 8
    - Fluid velocity component at 10
    - Tractions at 0, 3, 5, 9
    """

    @pytest.fixture
    def adg(self):
        from kernels.equations.poroelastic import PoroelasticADERDG

        return _mk(PoroelasticADERDG, order=4)

    def test_numberOfQuantities(self, adg):
        assert adg.numberOfQuantities() == 13

    def test_extractVelocities_has_4_components(self, adg):
        v = adg.extractVelocities()
        assert v.shape == (4, 13)
        assert np.count_nonzero(v) == 4
        assert v[0, 6] == 1
        assert v[1, 7] == 1
        assert v[2, 8] == 1
        assert v[3, 10] == 1  # the fluid-velocity-z component

    def test_extractTractions(self, adg):
        t = adg.extractTractions()
        assert t.shape == (4, 13)
        assert np.count_nonzero(t) == 4
        for row, col in [(0, 0), (1, 3), (2, 5), (3, 9)]:
            assert t[row, col] == 1


# =============================================================================
# Sparsity patterns for flux/godunov/transformation
# =============================================================================


class TestLinearADERDGSpp:
    """Base class sparsity patterns — fully dense for elastic/acoustic."""

    @pytest.fixture
    def adg(self):
        from kernels.equations.elastic import ElasticADERDG

        return _mk(ElasticADERDG, order=4)

    def test_godunov_is_all_ones(self, adg):
        spp = adg.godunov_spp()
        assert spp.dtype == bool
        assert spp.shape == (9, 9)
        assert np.all(spp)

    def test_flux_solver_is_all_ones(self, adg):
        spp = adg.flux_solver_spp()
        assert np.all(spp)

    def test_transformation_spp_is_square(self, adg):
        spp = adg.transformation_spp()
        assert spp.shape[0] == spp.shape[1]

    def test_transformation_inv_spp_matches_godunov(self, adg):
        """Elastic base implementation: transformation_inv_spp == godunov_spp."""
        np.testing.assert_array_equal(adg.transformation_inv_spp(), adg.godunov_spp())


class TestViscoelasticSpp:
    """Viscoelastic has BLOCK-structured sparsity patterns that depend on
    numberOfMechanisms. Miscounting blocks here breaks multi-mechanism
    simulations silently.
    """

    @pytest.fixture
    def adg(self):
        # viscoelastic needs numberOfElasticQuantities attribute set
        from kernels.equations.viscoelastic import ViscoelasticADERDG

        return _mk(
            ViscoelasticADERDG,
            order=4,
            numberOfMechanisms=3,
            numberOfElasticQuantities=9,
        )

    def test_numberOfQuantities_scales_with_mechanisms(self):
        """n_q = 9 + 6 * numberOfMechanisms — off-by-one here breaks memory."""
        from kernels.equations.viscoelastic import ViscoelasticADERDG

        for m in [0, 1, 3, 5, 10]:
            adg = _mk(
                ViscoelasticADERDG,
                order=4,
                numberOfMechanisms=m,
                numberOfElasticQuantities=9,
            )
            assert adg.numberOfQuantities() == 9 + 6 * m

    def test_godunov_spp_block_structure(self, adg):
        """First 9 rows are true, rest are false — coupling is only via
        the elastic rows (Kaeser & Dumbser 2006)."""
        spp = adg.godunov_spp()
        # Total quantities: 9 + 6*3 = 27
        assert spp.shape == (27, 27)
        assert np.all(spp[0:9, :])  # elastic rows: all ones
        assert not np.any(spp[9:, :])  # mechanism rows: all zeros

    def test_transformation_spp_diagonal_blocks(self, adg):
        """Structure:
        [6x6 | 0 | 0 | 0 | 0]
        [0 | 3x3 | 0 | 0 | 0]
        [0 | 0 | 6x6 | 0 | 0]   <- mechanism 0
        [0 | 0 | 0 | 6x6 | 0]   <- mechanism 1
        [0 | 0 | 0 | 0 | 6x6]   <- mechanism 2
        """
        spp = adg.transformation_spp()
        assert spp.shape == (27, 27)
        # Block (0..6, 0..6) is all ones, but (0..6, 6..9) is zero
        assert np.all(spp[0:6, 0:6])
        assert not np.any(spp[0:6, 6:9])
        # Block (6..9, 6..9) is all ones
        assert np.all(spp[6:9, 6:9])
        # Mechanism blocks along the diagonal starting at 9, 15, 21
        for m in range(3):
            off = 9 + m * 6
            assert np.all(spp[off : off + 6, off : off + 6])
            # Off-diagonal mechanism blocks are zero
            for m2 in range(3):
                if m != m2:
                    off2 = 9 + m2 * 6
                    assert not np.any(spp[off : off + 6, off2 : off2 + 6])

    def test_transformation_inv_equals_transformation(self, adg):
        """For viscoelastic: transformation_inv_spp == transformation_spp."""
        np.testing.assert_array_equal(
            adg.transformation_inv_spp(), adg.transformation_spp()
        )


# =============================================================================
# Poroelastic choose() — binomial coefficient helper
# =============================================================================


class TestPoroelasticChoose:
    """poroelastic.choose(n, k) — binomial, used for STP basis functions."""

    @pytest.fixture
    def choose(self):
        from kernels.equations.poroelastic import choose

        return choose

    def test_small_values(self, choose):
        assert choose(5, 2) == 10
        assert choose(10, 3) == 120
        assert choose(8, 0) == 1
        assert choose(8, 8) == 1

    @pytest.mark.parametrize(
        "n,k",
        [
            (5, 2),
            (8, 3),
            (10, 5),
            (12, 4),
            (20, 7),
        ],
    )
    def test_matches_math_comb(self, choose, n, k):
        """Local implementation must agree with stdlib math.comb."""
        assert choose(n, k) == comb(n, k)

    @pytest.mark.parametrize("n", [5, 10, 15])
    def test_pascal_symmetry(self, choose, n):
        """C(n, k) == C(n, n-k)"""
        for k in range(0, n + 1):
            assert choose(n, k) == choose(n, n - k)

    def test_pascal_triangle_identity(self, choose):
        """C(n+1, k+1) == C(n, k) + C(n, k+1)"""
        for n in range(3, 10):
            for k in range(0, n):
                assert choose(n + 1, k + 1) == choose(n, k) + choose(n, k + 1)
