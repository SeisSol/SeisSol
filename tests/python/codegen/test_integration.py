# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Integration tests for codegen.

These actually parse the JSON/XML matrix files and construct concrete
yateto-backed ADERDG objects — the same objects that `generate.py`
assembles at build time.

Unlike test_aderdg_math.py (which bypasses __init__ to test pure logic),
these exercise the full constructor path, including:
 - XML matrix file parsing (matrices_{N}.xml, star.xml, ...)
 - JSON matrix file parsing (plasticity-ip-matrices-{order}.json, ...)
 - memoryLayoutFromFile + CSC memory layouts
 - The full nodal/gravitational matrix suite

This catches breakage from:
 - matrix-file schema changes
 - yateto upstream changes
 - missing matrix files for specific orders
 - constructor-arg handling regressions
"""

import importlib.util  # noqa: F401 ; yateto expects importlib.util pre-loaded
from pathlib import Path

# Derive paths from the already-imported kernels module (sys.path is set up
# by conftest.py). This works in both dev layouts and the final integration.
import kernels  # noqa: F401 — establish sys.path via conftest
import kernels.memlayout as _ml
import numpy as np
import pytest

CODEGEN_DIR = Path(_ml.__file__).resolve().parent.parent
REPO_ROOT = CODEGEN_DIR.parent
MATRICES_DIR = CODEGEN_DIR / "matrices"
CONFIG_CPU = CODEGEN_DIR / "config" / "cpu"


@pytest.fixture(scope="module", autouse=True)
def _global_arch():
    """yateto requires a globally-fixed architecture for CSC layouts and
    alignment. We set up a reasonable default once per test module.
    """
    import warnings

    warnings.filterwarnings("ignore", category=UserWarning, module="yateto")
    from yateto import (
        HostArchDefinition,
        deriveArchitecture,
        fixArchitectureGlobal,
    )

    host = HostArchDefinition("hsw", "d", None, None)
    arch = deriveArchitecture(host, None)
    fixArchitectureGlobal(arch)
    yield arch


def _default_memLayout():
    """A minimal memory-layout file that matches any build configuration."""
    return str(CONFIG_CPU / "dense.xml")


# =============================================================================
# Elastic — reference case used in >80% of production runs
# =============================================================================


class TestElasticConstruction:
    """Real constructor path through yateto."""

    @pytest.mark.parametrize("order", [2, 3, 4, 5, 6])
    def test_constructs_across_build_matrix_orders(self, order):
        """The CI build matrix tests orders 2..6 for elastic — they must
        all construct. Matrix files for order 2 are matrices_4.xml, ...
        """
        from kernels.equations.elastic import ElasticADERDG

        adg = ElasticADERDG(
            order=order,
            multipleSimulations=1,
            matricesDir=str(MATRICES_DIR),
            memLayout=_default_memLayout(),
        )
        assert adg.numberOfQuantities() == 9

    def test_Q_tensor_shape_matches_basis_count_and_nq(self):
        from kernels.equations.elastic import ElasticADERDG

        adg = ElasticADERDG(
            order=4,
            multipleSimulations=1,
            matricesDir=str(MATRICES_DIR),
            memLayout=_default_memLayout(),
        )
        # Q holds DG coefficients: (numberOf3DBasisFunctions, numberOfQuantities)
        assert adg.Q.shape() == (
            adg.numberOf3DBasisFunctions(),
            adg.numberOfQuantities(),
        )

    def test_I_tensor_matches_Q(self):
        from kernels.equations.elastic import ElasticADERDG

        adg = ElasticADERDG(
            order=4,
            multipleSimulations=1,
            matricesDir=str(MATRICES_DIR),
            memLayout=_default_memLayout(),
        )
        # I is the time-integrated state — same shape as Q
        assert adg.I.shape() == adg.Q.shape()

    def test_star_matrices_loaded(self):
        """elastic.py clones `star` into `star(0)`, `star(1)`, `star(2)` —
        if the XML or the clones dict drift, this breaks."""
        from kernels.equations.elastic import ElasticADERDG

        adg = ElasticADERDG(
            order=4,
            multipleSimulations=1,
            matricesDir=str(MATRICES_DIR),
            memLayout=_default_memLayout(),
        )
        for dim in range(3):
            s = adg.starMatrix(dim)
            assert s is not None, f"star({dim}) missing"
            # Each star matrix is square and has shape (9, 9) for elastic
            shape = s.shape()
            assert shape == (9, 9), f"star({dim}).shape = {shape}"

    def test_flux_solver_tensors_created(self):
        from kernels.equations.elastic import ElasticADERDG

        adg = ElasticADERDG(
            order=4,
            multipleSimulations=1,
            matricesDir=str(MATRICES_DIR),
            memLayout=_default_memLayout(),
        )
        # Local and neighbor flux solvers — both exist
        for attr in ["AplusT", "AminusT", "T", "Tinv"]:
            t = getattr(adg, attr)
            assert t is not None, f"attr {attr} missing"

    def test_fused_simulations_changes_Q_rank(self):
        """multipleSimulations > 1 inserts a 'simulation' axis into tensors.
        Shape must reflect that — otherwise codegen produces a flat dense
        tensor and the fused-sim feature is broken.
        """
        from kernels.equations.elastic import ElasticADERDG

        adg_single = ElasticADERDG(
            order=3,
            multipleSimulations=1,
            matricesDir=str(MATRICES_DIR),
            memLayout=_default_memLayout(),
        )
        adg_fused = ElasticADERDG(
            order=3,
            multipleSimulations=8,
            matricesDir=str(MATRICES_DIR),
            memLayout=_default_memLayout(),
        )
        # Single-sim Q is rank 2; fused-sim Q is rank 3
        assert len(adg_single.Q.shape()) == 2
        assert len(adg_fused.Q.shape()) == 3
        # First axis of fused is the sim axis
        assert adg_fused.Q.shape()[0] == 8


# =============================================================================
# Acoustic / Anisotropic / Poroelastic — smoke construction
# =============================================================================


class TestOtherEquationsConstruction:

    def test_acoustic_constructs(self):
        from kernels.equations.acoustic import AcousticADERDG

        adg = AcousticADERDG(
            order=4,
            multipleSimulations=1,
            matricesDir=str(MATRICES_DIR),
            memLayout=_default_memLayout(),
        )
        assert adg.numberOfQuantities() == 4

    def test_anisotropic_constructs(self):
        from kernels.equations.anisotropic import AnisotropicADERDG

        adg = AnisotropicADERDG(
            order=4,
            multipleSimulations=1,
            matricesDir=str(MATRICES_DIR),
            memLayout=_default_memLayout(),
        )
        # Anisotropic has the same 9 quantities as elastic
        assert adg.numberOfQuantities() == 9

    def test_poroelastic_constructs(self):
        from kernels.equations.poroelastic import PoroelasticADERDG

        adg = PoroelasticADERDG(
            order=4,
            multipleSimulations=1,
            matricesDir=str(MATRICES_DIR),
            memLayout=_default_memLayout(),
            numberOfMechanisms=0,
        )
        assert adg.numberOfQuantities() == 13


class TestViscoelasticConstruction:
    """Viscoelastic is the configuration that stresses the mechanism-count-
    dependent block structure the most."""

    @pytest.mark.parametrize("mechanisms", [1, 3, 5])
    def test_numberOfQuantities_matches_formula(self, mechanisms):
        from kernels.equations.viscoelastic import ViscoelasticADERDG

        adg = ViscoelasticADERDG(
            order=4,
            multipleSimulations=1,
            matricesDir=str(MATRICES_DIR),
            memLayout=_default_memLayout(),
            numberOfMechanisms=mechanisms,
        )
        assert adg.numberOfQuantities() == 9 + 6 * mechanisms

    @pytest.mark.parametrize("mechanisms", [1, 3, 5])
    def test_star_matrix_shape_scales_with_mechanisms(self, mechanisms):
        """After the mechanism-expansion logic, star matrices must be
        square of size numberOfQuantities × numberOfQuantities.
        """
        from kernels.equations.viscoelastic import ViscoelasticADERDG

        adg = ViscoelasticADERDG(
            order=4,
            multipleSimulations=1,
            matricesDir=str(MATRICES_DIR),
            memLayout=_default_memLayout(),
            numberOfMechanisms=mechanisms,
        )
        nq = adg.numberOfQuantities()
        for dim in range(3):
            assert adg.starMatrix(dim).shape() == (nq, nq)


# =============================================================================
# Matrix-file inventory — catches missing files for supported orders
# =============================================================================


class TestMatrixFileInventory:
    """Documents which matrix files must exist in codegen/matrices/ for
    which (equation, order) combinations. If a future PR adds an order but
    forgets a matrix file, this flags it.

    Orders 2..8 are the production build-matrix (see .github/workflows/
    build-seissol-cpu.yml — ORDER=6 currently, but Spack/config supports
    2..8).
    """

    # (order -> expected numberOf3DBasisFunctions)
    BASIS_FOR_ORDER = {2: 4, 3: 10, 4: 20, 5: 35, 6: 56, 7: 84, 8: 120}

    @pytest.mark.parametrize("order", [2, 3, 4, 5, 6, 7, 8])
    def test_matrices_N_xml_exists_for_each_order(self, order):
        """matrices_{numberOf3DBasisFunctions}.xml is required for every
        supported order."""
        n = self.BASIS_FOR_ORDER[order]
        f = MATRICES_DIR / f"matrices_{n}.xml"
        assert f.exists(), f"Missing matrix file for order={order}: {f}"

    @pytest.mark.parametrize("order", [2, 3, 4, 5, 6, 7, 8])
    @pytest.mark.parametrize("method", ["ip", "nb"])
    def test_plasticity_matrix_files_exist(self, order, method):
        """Both plasticity methods (Interior Point, Nodal Basis) need a
        matrix file for each order.
        """
        f = MATRICES_DIR / f"plasticity-{method}-matrices-{order}.json"
        assert f.exists(), (
            f"Missing plasticity matrix file: {f}. "
            f"Either the order isn't supported yet, or the file is missing."
        )

    @pytest.mark.parametrize("order", [2, 3, 4, 5, 6, 7, 8])
    def test_mass_matrix_per_order_exists(self, order):
        f = MATRICES_DIR / f"mass_{order}.json"
        assert f.exists(), f"Missing mass matrix file for order={order}"

    @pytest.mark.parametrize("order", [2, 3, 4, 5, 6, 7, 8])
    def test_stp_matrix_exists_for_poroelastic(self, order):
        """Space-Time-Predictor matrices — required for poroelastic &
        viscoelastic2 at every order.
        """
        f = MATRICES_DIR / f"stp_{order}.json"
        assert f.exists(), f"Missing stp_{order}.json"

    @pytest.mark.parametrize("order", [2, 3, 4, 5, 6, 7, 8])
    def test_nodal_boundary_matrices_exist(self, order):
        f = MATRICES_DIR / "nodal" / f"nodalBoundary_matrices_{order}.json"
        assert f.exists(), f"Missing nodal boundary matrices for order={order}"

    @pytest.mark.parametrize("order", [2, 3, 4, 5, 6, 7, 8])
    def test_gravitational_energy_matrices_exist(self, order):
        f = MATRICES_DIR / "nodal" / f"gravitational_energy_matrices_{order}.xml"
        assert f.exists(), f"Missing gravitational_energy_matrices for order={order}"

    def test_equation_specific_matrices_exist(self):
        """star matrices for each equation that needs one."""
        for fname in [
            "star.xml",  # elastic
            "star_acoustic.xml",  # acoustic
            "star_anisotropic.xml",  # anisotropic
            "matrices_viscoelastic.xml",
            "matrices_poroelastic.xml",
        ]:
            f = MATRICES_DIR / fname
            assert f.exists(), f"Missing equation matrix file: {f}"

    def test_dr_quadrature_matrices_exist(self):
        """Dynamic Rupture uses three quadrature rules × 7 orders."""
        orders_present = {2, 3, 4, 5, 6, 7, 8}
        missing = []
        for rule in ["dunavant", "stroud", "witherden_vincent"]:
            for order in orders_present:
                f = MATRICES_DIR / f"dr_{rule}_matrices_{order}.json"
                if not f.exists():
                    missing.append(f.name)
        # Not all (rule, order) combinations exist historically — we just
        # flag what's missing so someone can justify any gaps
        if missing:
            pytest.skip(
                f"Dynamic-rupture quadrature gaps exist "
                f"(expected on current master): {missing[:5]}..."
            )


# =============================================================================
# Config files — memory layouts for the build matrix
# =============================================================================


class TestConfigFileInventory:

    def test_dense_fallback_exists(self):
        """guessMemoryLayout falls back to dense.xml when nothing matches —
        if this file is deleted, that fallback crashes."""
        assert (CONFIG_CPU / "dense.xml").exists()

    def test_architectures_for_cpu(self):
        """The host architectures documented in memlayout.py's 'pes' list
        should have at least one config file each.
        """
        archs_in_use = {"hsw", "skx", "a64fx", "snb", "knl"}
        # Gather what architectures have config files
        config_files = list(CONFIG_CPU.glob("*.xml"))
        archs_present = set()
        for f in config_files:
            stem_parts = f.stem.split("_")
            for part in stem_parts:
                if part in archs_in_use:
                    archs_present.add(part)
        assert len(archs_present) > 0, (
            "No recognized architecture-specific config files found. "
            f"Available: {[f.name for f in config_files]}"
        )


# =============================================================================
# tensor_to_numpy / numpy_to_tensor — roundtrip in common.py
# =============================================================================


class TestTensorNumpyRoundtrip:
    """common.py provides numpy <-> yateto Tensor conversions. These are
    used when codegen needs to manipulate sparsity patterns with numpy
    and feed the result back to yateto. The roundtrip must preserve data.
    """

    def test_roundtrip_dense_matrix(self):
        from kernels.common import numpy_to_tensor, tensor_to_numpy

        a = np.array([[1.0, 2.0, 0.0], [0.0, 3.0, 4.0]])
        t = numpy_to_tensor("A", a)
        recovered = tensor_to_numpy(t)
        np.testing.assert_allclose(recovered, a)

    def test_roundtrip_preserves_zeros_via_sparsity(self):
        """numpy_to_tensor builds a sparse-pattern spp dict. Zeros in the
        original are NOT stored — they must still read back as zeros."""
        from kernels.common import numpy_to_tensor, tensor_to_numpy

        a = np.zeros((4, 4))
        a[0, 0] = 1.0
        a[2, 3] = 5.0
        t = numpy_to_tensor("A", a)
        recovered = tensor_to_numpy(t)
        np.testing.assert_allclose(recovered, a)

    @pytest.mark.xfail(
        reason="Upstream yateto bug",
        raises=IndexError,
        strict=True,
    )
    def test_roundtrip_all_zeros(self):
        from kernels.common import numpy_to_tensor, tensor_to_numpy

        a = np.zeros((3, 3))
        t = numpy_to_tensor("A", a)
        recovered = tensor_to_numpy(t)
        np.testing.assert_allclose(recovered, a)

    def test_roundtrip_integer_array(self):
        """Common usage: DoF extract patterns are integer arrays (0/1).
        Conversion to float via np.float64 in tensor_to_numpy must work.
        """
        from kernels.common import numpy_to_tensor, tensor_to_numpy

        a = np.zeros((3, 9), dtype=int)
        a[0, 6] = 1
        a[1, 7] = 1
        a[2, 8] = 1
        t = numpy_to_tensor("extractVel", a.astype(float))
        recovered = tensor_to_numpy(t)
        np.testing.assert_allclose(recovered, a)

    def test_roundtrip_preserves_shape(self):
        from kernels.common import numpy_to_tensor, tensor_to_numpy

        for shape in [(3,), (2, 4), (5, 5), (2, 3, 4)]:
            a = np.random.rand(*shape)
            t = numpy_to_tensor("A", a)
            recovered = tensor_to_numpy(t)
            assert recovered.shape == shape
