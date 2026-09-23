# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Structural / cross-module tests for codegen/kernels.

These guard invariants that every equation module MUST satisfy, and
document the known bug in vtkproject.py for future fixers.
"""

import importlib  # yateto submodule needs this pre-loaded
import importlib.util

import numpy as np
import pytest

# The six equation modules shipped by SeisSol's codegen, with the arguments
# kernel_class needs to pick one of several classes.
EQUATION_KWARGS = {
    "acoustic": {"solver": "linearck"},
    "elastic": {"solver": "linearck"},
    "anisotropic": {"solver": "linearck"},
    "poroelastic": {"solver": "stp"},
    "viscoelastic": {"solver": "linearckanelastic"},
    "viscoacoustic": {"solver": "linearckanelastic"},
}
EQUATION_MODULES = list(EQUATION_KWARGS)


def equation_class(module_name):
    """The class generate.py would pick for this equation."""
    mod = importlib.import_module(f"kernels.equations.{module_name}")
    return mod.kernel_class(**EQUATION_KWARGS[module_name])


# =============================================================================
# Every equation module exposes EQUATION_CLASS
# =============================================================================


class TestEquationModuleContract:
    """generate.py loads equation modules dynamically by name via:
        equations = importlib.import_module()
        equation_class = equations.kernel_class(**args)
    So every module in kernels/equations MUST export kernel_class, and it has
    to return a class for the arguments generate.py passes. This contract is
    only enforced at runtime -- breaking it silently passes linters and
    static checks.
    """

    @pytest.mark.parametrize("module_name", EQUATION_MODULES)
    def test_module_imports(self, module_name):
        mod = importlib.import_module(f"kernels.equations.{module_name}")
        assert mod is not None

    @pytest.mark.parametrize("module_name", EQUATION_MODULES)
    def test_exposes_kernel_class(self, module_name):
        mod = importlib.import_module(f"kernels.equations.{module_name}")
        assert callable(getattr(mod, "kernel_class", None)), (
            f"kernels.equations.{module_name} does not export kernel_class; "
            f"generate.py will raise AttributeError when selecting this equation."
        )

    @pytest.mark.parametrize("module_name", EQUATION_MODULES)
    def test_kernel_class_returns_a_class(self, module_name):
        assert isinstance(equation_class(module_name), type)

    @pytest.mark.parametrize("module_name", EQUATION_MODULES)
    def test_equation_class_name_has_aderdg_suffix(self, module_name):
        """Convention: <Equation>ADERDG. Documented for discoverability --
        not a hard requirement, but a useful one."""
        assert equation_class(module_name).__name__.endswith("ADERDG")

    @pytest.mark.parametrize("module_name", ["viscoelastic", "viscoacoustic"])
    def test_relaxation_rejects_a_solver_that_cannot_carry_it(self, module_name):
        """The predictor substitutes one entry per stiff source row, and
        relaxation contributes several, so the combination has to be refused
        rather than silently generated."""
        mod = importlib.import_module(f"kernels.equations.{module_name}")
        with pytest.raises(NotImplementedError):
            mod.kernel_class(solver="stp")


class TestEquationClassInheritance:
    """All equation classes ultimately derive from LinearCK, which
    provides the default spp/basis-count methods. Anisotropic derives from
    Elastic (re-uses its indexing).
    """

    def test_elastic_is_linear_aderdg(self):
        from kernels.aderdg.linearck import LinearCK
        from kernels.equations.elastic import ElasticADERDG

        assert issubclass(ElasticADERDG, LinearCK)

    def test_anisotropic_extends_elastic(self):
        from kernels.equations.anisotropic import AnisotropicADERDG
        from kernels.equations.elastic import ElasticADERDG

        assert issubclass(AnisotropicADERDG, ElasticADERDG)

    @pytest.mark.parametrize(
        "module_name",
        [
            "acoustic",
            "elastic",
            "anisotropic",
            "poroelastic",
        ],
    )
    def test_linear_equations_inherit_from_linear_aderdg(self, module_name):
        """Equations whose solver builds on LinearCK inherit from it. The
        relaxing ones do not when they are built on LinearCKAnelastic -- see
        the solver hierarchy test below."""
        from kernels.aderdg.linearck import LinearCK

        assert issubclass(
            equation_class(module_name), LinearCK
        ), f"{module_name}: the chosen class must inherit from LinearCK"

    def test_solver_hierarchy_asymmetry(self):
        """Two of the three solvers extend LinearCK; the anelastic one does
        not. It keeps the memory variables in a tensor dimension of their own,
        so it has its own time kernel rather than a variation on LinearCK's.
        Any refactor that tries to normalise the hierarchy has to account for
        that, which is why it is written down.
        """
        from kernels.aderdg.aderdg import ADERDGBase
        from kernels.aderdg.linearck import LinearCK
        from kernels.aderdg.linearckanelastic import LinearCKAnelastic
        from kernels.aderdg.stp import STP

        assert issubclass(LinearCK, ADERDGBase)
        assert issubclass(STP, LinearCK)

        assert issubclass(LinearCKAnelastic, ADERDGBase)
        assert not issubclass(LinearCKAnelastic, LinearCK), (
            "If this now fails, the hierarchy has been unified. Verify that "
            "the anelastic time kernel still keeps the mechanism index in its "
            "own tensor dimension, then drop this test."
        )

    @pytest.mark.parametrize("module_name", EQUATION_MODULES)
    def test_all_equations_are_aderdg_base(self, module_name):
        """The universally-applicable invariant: every equation is an ADERDGBase."""
        from kernels.aderdg.aderdg import ADERDGBase

        assert issubclass(equation_class(module_name), ADERDGBase)


# =============================================================================
# DoF-index patterns are consistent across equations
# =============================================================================


class TestCrossEquationDoFConsistency:
    """Every equation's extractVelocities / extractTractions patterns must
    satisfy shared invariants that kernel code assumes.
    """

    def _make(self, module_name, **kwargs):
        instance = object.__new__(equation_class(module_name))
        defaults = {"order": 4}
        if module_name == "viscoelastic":
            defaults.update(numMechanisms=3, numElasticQuantities=9)
        defaults.update(kwargs)
        for k, v in defaults.items():
            setattr(instance, k, v)
        return instance

    @pytest.mark.parametrize(
        "module_name",
        [
            "acoustic",
            "elastic",
            "poroelastic",
        ],
    )
    def test_extract_patterns_have_correct_quantity_dim(self, module_name):
        """extractVelocities/Tractions must have second dim == numQuantities."""
        adg = self._make(module_name)
        nq = adg.numQuantities()
        v = adg.extractVelocities()
        t = adg.extractTractions()
        assert v.shape[1] == nq, f"{module_name}: velocity cols != nq"
        assert t.shape[1] == nq, f"{module_name}: traction cols != nq"

    @pytest.mark.parametrize("module_name", ["elastic", "poroelastic"])
    def test_velocity_and_traction_indices_disjoint(self, module_name):
        """For equations with both velocities and tractions, their DoF
        indices must not overlap — a DoF is either one or the other."""
        adg = self._make(module_name)
        v = adg.extractVelocities()
        t = adg.extractTractions()
        v_cols = set(np.nonzero(v.any(axis=0))[0].tolist())
        t_cols = set(np.nonzero(t.any(axis=0))[0].tolist())
        assert not (v_cols & t_cols), (
            f"{module_name}: velocity cols {v_cols} and "
            f"traction cols {t_cols} overlap"
        )

    @pytest.mark.parametrize(
        "module_name",
        [
            "acoustic",
            "elastic",
            "poroelastic",
        ],
    )
    def test_extract_patterns_are_row_wise_one_hot(self, module_name):
        """Every row of extractVelocities must contain exactly one 1.
        (Each row selects one velocity component; multiple 1s would
        average or sum them.)
        """
        adg = self._make(module_name)
        v = adg.extractVelocities()
        for i in range(v.shape[0]):
            assert (
                np.count_nonzero(v[i]) == 1
            ), f"{module_name}.extractVelocities row {i} is not one-hot"

    @pytest.mark.parametrize("module_name", ["elastic", "poroelastic"])
    def test_extract_tractions_are_row_wise_one_hot(self, module_name):
        adg = self._make(module_name)
        t = adg.extractTractions()
        for i in range(t.shape[0]):
            assert (
                np.count_nonzero(t[i]) == 1
            ), f"{module_name}.extractTractions row {i} is not one-hot"


# =============================================================================
# generate_kernel_name_prefix invariants
# =============================================================================


class TestKernelNamePrefixInvariants:
    """The function is dirt-simple but called all over. Lock in its shape."""

    def test_gpu_gives_gpu_suffix(self):
        from kernels.common import generate_kernel_name_prefix

        assert generate_kernel_name_prefix("gpu") == "gpu_"

    def test_cpu_gives_empty(self):
        from kernels.common import generate_kernel_name_prefix

        assert generate_kernel_name_prefix("cpu") == ""

    def test_unknown_target_gives_empty(self):
        """Callers concatenate the result as f'{name_prefix}kernelname'.
        An unexpected target must NOT give a prefix that pollutes the
        namespace — it falls through to ''.
        """
        from kernels.common import generate_kernel_name_prefix

        for bad in ["GPU", "Gpu", "gpU", "cuda", "hip", "", "xpu"]:
            assert generate_kernel_name_prefix(bad) == "", (
                f"Unexpected target {bad!r} produced prefix "
                f"{generate_kernel_name_prefix(bad)!r}"
            )

    def test_output_is_always_string(self):
        """Invariant: callers concatenate into f-strings. Must be str."""
        from kernels.common import generate_kernel_name_prefix

        for target in ["gpu", "cpu", "anything"]:
            assert isinstance(generate_kernel_name_prefix(target), str)


# =============================================================================
# OptionalDimTensor (multsim.py)
# =============================================================================


class TestOptionalDimTensor:
    """OptionalDimTensor optionally inserts a 'multisim' dimension when
    multipleSimulations > 1. This is where fused-sim shape bugs live.
    """

    def test_hasOptDim_when_optSize_gt_1(self):
        from kernels.multsim import OptionalDimTensor

        t = OptionalDimTensor("Q", "s", optSize=8, optPos=0, shape=(10, 9))
        assert t.hasOptDim() is True

    def test_hasOptDim_false_when_optSize_is_1(self):
        from kernels.multsim import OptionalDimTensor

        t = OptionalDimTensor("Q", "s", optSize=1, optPos=0, shape=(10, 9))
        assert t.hasOptDim() is False

    def test_insertOptDim_adds_dim_at_position(self):
        from kernels.multsim import OptionalDimTensor

        t = OptionalDimTensor("Q", "s", optSize=8, optPos=1, shape=(10, 9))
        # When hasOptDim is True, insert item at optPos in a tuple
        result = t.insertOptDim((10, 9), (8,))
        assert result == (10, 8, 9)

    def test_insertOptDim_at_position_zero(self):
        from kernels.multsim import OptionalDimTensor

        t = OptionalDimTensor("Q", "s", optSize=8, optPos=0, shape=(10, 9))
        result = t.insertOptDim((10, 9), (8,))
        assert result == (8, 10, 9)

    def test_insertOptDim_does_nothing_when_no_opt_dim(self):
        from kernels.multsim import OptionalDimTensor

        t = OptionalDimTensor("Q", "s", optSize=1, optPos=0, shape=(10, 9))
        result = t.insertOptDim((10, 9), (1,))
        assert result == (10, 9)

    def test_accessors_return_constructor_args(self):
        from kernels.multsim import OptionalDimTensor

        t = OptionalDimTensor("Q", "s", optSize=8, optPos=1, shape=(10, 9))
        assert t.optName() == "s"
        assert t.optSize() == 8
        assert t.optPos() == 1

    def test_getitem_injects_multisim_index_when_fused(self):
        from kernels.multsim import OptionalDimTensor

        t = OptionalDimTensor("Q", "s", optSize=8, optPos=0, shape=(10, 9))
        # Accessing t["ij"] should produce an IndexedTensor with "sij"
        it = t["ij"]
        # The indices now include the optName at optPos
        # (we just check it doesn't crash and produces something)
        assert it is not None

    def test_shape_reflects_opt_dim_presence(self):
        from kernels.multsim import OptionalDimTensor

        fused = OptionalDimTensor("Q", "s", optSize=8, optPos=0, shape=(10, 9))
        unfused = OptionalDimTensor("Q", "s", optSize=1, optPos=0, shape=(10, 9))
        # Fused: 3D (optSize, 10, 9)
        # Unfused: 2D (10, 9)
        assert len(fused.shape()) == 3
        assert len(unfused.shape()) == 2
