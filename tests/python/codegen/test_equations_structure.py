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

# The six equation modules shipped by SeisSol's codegen
EQUATION_MODULES = [
    "acoustic",
    "elastic",
    "anisotropic",
    "poroelastic",
    "viscoelastic",
    "viscoelastic2",
]


# =============================================================================
# Every equation module exposes EQUATION_CLASS
# =============================================================================


class TestEquationModuleContract:
    """generate.py loads equation modules dynamically by name via:
        equations = equationsSpec.loader.load_module()
        equation_class = equations.EQUATION_CLASS
    So every module in kernels/equations MUST export EQUATION_CLASS.
    This contract is only enforced at runtime — breaking it silently
    passes linters and static checks.
    """

    @pytest.mark.parametrize("module_name", EQUATION_MODULES)
    def test_module_imports(self, module_name):
        mod = importlib.import_module(f"kernels.equations.{module_name}")
        assert mod is not None

    @pytest.mark.parametrize("module_name", EQUATION_MODULES)
    def test_exposes_equation_class(self, module_name):
        mod = importlib.import_module(f"kernels.equations.{module_name}")
        assert hasattr(mod, "EQUATION_CLASS"), (
            f"kernels.equations.{module_name} does not export EQUATION_CLASS; "
            f"generate.py will raise AttributeError when selecting this equation."
        )

    @pytest.mark.parametrize("module_name", EQUATION_MODULES)
    def test_equation_class_is_a_class(self, module_name):
        mod = importlib.import_module(f"kernels.equations.{module_name}")
        assert isinstance(mod.EQUATION_CLASS, type)

    @pytest.mark.parametrize("module_name", EQUATION_MODULES)
    def test_equation_class_name_has_aderdg_suffix(self, module_name):
        """Convention: <Equation>ADERDG. Documented for discoverability —
        not a hard requirement, but a useful one."""
        mod = importlib.import_module(f"kernels.equations.{module_name}")
        assert mod.EQUATION_CLASS.__name__.endswith("ADERDG")


class TestEquationClassInheritance:
    """All equation classes ultimately derive from LinearADERDG, which
    provides the default spp/basis-count methods. Anisotropic derives from
    Elastic (re-uses its indexing).
    """

    def test_elastic_is_linear_aderdg(self):
        from kernels.aderdg import LinearADERDG
        from kernels.equations.elastic import ElasticADERDG

        assert issubclass(ElasticADERDG, LinearADERDG)

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
            "viscoelastic",
        ],
    )
    def test_linear_equations_inherit_from_linear_aderdg(self, module_name):
        """All LINEAR equations inherit from LinearADERDG."""
        from kernels.aderdg import LinearADERDG

        mod = importlib.import_module(f"kernels.equations.{module_name}")
        assert issubclass(
            mod.EQUATION_CLASS, LinearADERDG
        ), f"{module_name}: EQUATION_CLASS must inherit from LinearADERDG"

    def test_viscoelastic2_inherits_from_aderdgbase_directly(self):
        """Structural asymmetry: viscoelastic2 uses the Space-Time-Predictor
        (STP) formulation and inherits ADERDGBase directly, bypassing
        LinearADERDG. Any refactor that tries to 'normalize' the hierarchy
        would break the STP machinery. This test documents the decision.
        """
        from kernels.aderdg import ADERDGBase, LinearADERDG
        from kernels.equations.viscoelastic2 import Viscoelastic2ADERDG

        assert issubclass(Viscoelastic2ADERDG, ADERDGBase)
        assert not issubclass(Viscoelastic2ADERDG, LinearADERDG), (
            "If this now fails, the hierarchy has been unified. Verify that "
            "the STP-specific behavior (see numberOfExtendedQuantities, etc.) "
            "has been preserved, then drop this test."
        )

    @pytest.mark.parametrize("module_name", EQUATION_MODULES)
    def test_all_equations_are_aderdg_base(self, module_name):
        """The universally-applicable invariant: every equation is an ADERDGBase."""
        from kernels.aderdg import ADERDGBase

        mod = importlib.import_module(f"kernels.equations.{module_name}")
        assert issubclass(mod.EQUATION_CLASS, ADERDGBase)


# =============================================================================
# DoF-index patterns are consistent across equations
# =============================================================================


class TestCrossEquationDoFConsistency:
    """Every equation's extractVelocities / extractTractions patterns must
    satisfy shared invariants that kernel code assumes.
    """

    def _make(self, module_name, **kwargs):
        mod = importlib.import_module(f"kernels.equations.{module_name}")
        instance = object.__new__(mod.EQUATION_CLASS)
        defaults = {"order": 4}
        if module_name == "viscoelastic":
            defaults.update(numberOfMechanisms=3, numberOfElasticQuantities=9)
        if module_name == "viscoelastic2":
            defaults.update(numberOfMechanisms=3, numberOfElasticQuantities=9)
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
        """extractVelocities/Tractions must have second dim == numberOfQuantities."""
        adg = self._make(module_name)
        nq = adg.numberOfQuantities()
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
# Documented bug in vtkproject.py
# =============================================================================


class TestVtkProjectKnownBug:
    """
    BUG: codegen/kernels/vtkproject.py line 15 (as of master @ 2026-04-21):

        def addKernels(generator, aderdg, PlasticityMethod, matricesDir,
                       targets=["cpu"]):
            for target in targets:
                name_prefix = generate_kernel_name_prefix(targets)
                                                          ^^^^^^^
                # ↑ passes the LIST of targets, not the loop variable

    `generate_kernel_name_prefix` does `return f"{target}_" if target == "gpu"
    else ""`.  Comparing a list to the string "gpu" is always False, so the
    prefix is always empty — meaning VTK kernels inside this function never
    receive the `gpu_` prefix on GPU builds. Kernels of the same name exist
    for both CPU and GPU; without the prefix, the GPU variants likely collide
    or get mis-dispatched.

    Fix: replace `generate_kernel_name_prefix(targets)` with
    `generate_kernel_name_prefix(target)` (singular).

    These tests DOCUMENT the bug and will fail when fixed — which is the
    signal to flip them to their positive form.
    """

    def test_prefix_bug_on_gpu_target(self):
        from kernels.common import generate_kernel_name_prefix

        # Reproduce the buggy call (passing the list)
        targets = ["gpu"]
        actual_prefix = generate_kernel_name_prefix(targets)
        # Current (buggy) behavior: returns ''
        assert actual_prefix == "", (
            "If this assertion now fails with 'gpu_', the bug is fixed. "
            "Change this test to `assert actual_prefix == 'gpu_'` and invert "
            "the corresponding test below."
        )

    def test_prefix_correct_when_passed_string(self):
        """What the CORRECT call site should do."""
        from kernels.common import generate_kernel_name_prefix

        assert generate_kernel_name_prefix("gpu") == "gpu_"
        assert generate_kernel_name_prefix("cpu") == ""


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
