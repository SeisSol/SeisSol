# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""End-to-end smoke + golden-file tests for generate.py.

These invoke the REAL generator with minimal configurations and verify:
 - All the forward-files listed in generate.py are produced
 - The per-equation subfolder is created with its expected name pattern
 - The generated header files contain the kernel-family names the rest
   of SeisSol's C++ code expects

Design notes:
 - Uses `subprocess` so the test exercises the CLI path CMake uses,
   not a hot-wired Python function call.
 - Uses `tmp_path` per test, so golden-content tests don't interfere.
 - Module-scoped fixture caches one generator run across tests, since
   the full generation takes ~1–2 seconds.
"""

import importlib.util  # noqa: F401
import subprocess
import sys
from pathlib import Path

# Resolve codegen dir via the sys.path set up in conftest.py
import kernels.memlayout as _ml
import pytest

CODEGEN_DIR = Path(_ml.__file__).resolve().parent.parent
GENERATE = CODEGEN_DIR / "generate.py"


def _invoke_generate(
    outdir, equation="elastic", order=3, precision="d", multi_sims=1, mechanisms=0
):
    """Run generate.py with the given config. Returns CompletedProcess."""
    return subprocess.run(
        [
            sys.executable,
            str(GENERATE),
            "--equations",
            equation,
            "--matricesDir",
            str(CODEGEN_DIR / "matrices"),
            "--outputDir",
            str(outdir),
            "--host_arch",
            "hsw",
            "--order",
            str(order),
            "--precision",
            precision,
            "--numberOfMechanisms",
            str(mechanisms),
            "--memLayout",
            "auto",
            "--multipleSimulations",
            str(multi_sims),
            "--PlasticityMethod",
            "ip",
            "--gemm_tools",
            "none",
            "--drQuadRule",
            "dunavant",
            "--device_backend",
            "none",
        ],
        cwd=str(CODEGEN_DIR),
        capture_output=True,
        text=True,
        timeout=120,
    )


# =============================================================================
# Single cached run of the minimal config, reused across tests
# =============================================================================


@pytest.fixture(scope="module")
def generated_elastic_o3(tmp_path_factory):
    """Generate once per module: elastic, order=3, double precision."""
    outdir = tmp_path_factory.mktemp("gen_elastic_o3")
    result = _invoke_generate(outdir)
    if result.returncode != 0:
        pytest.fail(
            f"generate.py failed (code {result.returncode}):\n"
            f"stdout:\n{result.stdout[-2000:]}\n"
            f"stderr:\n{result.stderr[-2000:]}"
        )
    return outdir, result


# =============================================================================
# Existence of output files
# =============================================================================


class TestGeneratedFilesExist:
    """generate.py's `forward_files` list — we assert each one appears."""

    # Must match the literal list at the bottom of generate.py's main().
    # If the list drifts, this catches it.
    EXPECTED_FORWARD_FILES = {
        "init.h",
        "kernel.h",
        "tensor.h",
    }

    def test_forward_files_produced(self, generated_elastic_o3):
        outdir, _ = generated_elastic_o3
        produced = {p.name for p in outdir.iterdir() if p.is_file()}
        missing = self.EXPECTED_FORWARD_FILES - produced
        assert not missing, f"Missing forward files: {missing}"

    def test_equation_subfolder_name_pattern(self, generated_elastic_o3):
        """The per-equation subfolder is named
        `equation-<name>-<order>-<precision>[-f<multi>]`.
        Test covers elastic/3/double."""
        outdir, _ = generated_elastic_o3
        subdirs = [p for p in outdir.iterdir() if p.is_dir()]
        names = {p.name for p in subdirs}
        assert (
            "equation-elastic-3-double" in names
        ), f"Expected equation-elastic-3-double subfolder; got {names}"
        assert "general" in names

    def test_equation_subfolder_contains_per_equation_files(self, generated_elastic_o3):
        outdir, _ = generated_elastic_o3
        sub = outdir / "equation-elastic-3-double"
        files = {p.name for p in sub.iterdir() if p.is_file()}
        # The per-equation subfolder must contain its own copies of these
        for expected in [
            "init.h",
            "kernel.h",
            "tensor.h",
            "init.cpp",
            "kernel.cpp",
            "tensor.cpp",
            "test-kernel.cpp",
        ]:
            assert (
                expected in files
            ), f"Expected {expected} in equation subfolder, got {files}"

    def test_general_subfolder_exists(self, generated_elastic_o3):
        """The `general` subfolder holds kernels used only in initialization
        (see generate_general()). Double-precision only."""
        outdir, _ = generated_elastic_o3
        assert (outdir / "general").is_dir()

    def test_top_level_subroutine_generated(self, generated_elastic_o3):
        """subroutine.cpp/.h — GEMM routine cache produced at top level."""
        outdir, _ = generated_elastic_o3
        assert (outdir / "subroutine.cpp").exists()
        assert (outdir / "subroutine.h").exists()


# =============================================================================
# Content invariants: headers reference the correct kernel families
# =============================================================================


class TestGeneratedContent:
    """We don't snapshot-hash the full output (that would flip-flop on
    yateto upgrades). Instead we assert invariants about what MUST be
    present, which is what the C++ caller relies on.
    """

    def test_init_h_uses_iwyu_pragma(self, generated_elastic_o3):
        outdir, _ = generated_elastic_o3
        content = (outdir / "init.h").read_text()
        assert "IWYU pragma: begin_exports" in content
        assert "IWYU pragma: end_exports" in content

    def test_forward_files_include_equation_subdir(self, generated_elastic_o3):
        """The forward init.h should #include the equation-subdir's init.h."""
        outdir, _ = generated_elastic_o3
        content = (outdir / "init.h").read_text()
        assert "equation-elastic-3-double/init.h" in content
        assert "general/init.h" in content

    def test_kernel_h_declares_aderdg_kernels(self, generated_elastic_o3):
        """ADER-DG pipeline kernels must appear in kernel.h. If generate.py
        or the kernel-name contract changes, the C++ side breaks."""
        outdir, _ = generated_elastic_o3
        content = (outdir / "equation-elastic-3-double" / "kernel.h").read_text()
        # A handful of canonical kernel names
        for name in [
            "computeFluxSolverLocal",
            "computeFluxSolverNeighbor",
        ]:
            assert name in content, f"Expected kernel '{name}' in kernel.h"

    def test_test_kernel_cpp_not_empty(self, generated_elastic_o3):
        """test-kernel.cpp feeds the TESTING_GENERATED compile path."""
        outdir, _ = generated_elastic_o3
        content = (outdir / "equation-elastic-3-double" / "test-kernel.cpp").read_text()
        assert len(content.strip()) > 0
        # A handful of canonical kernel names
        for name in [
            "computeFluxSolverLocal",
            "computeFluxSolverNeighbor",
        ]:
            assert name in content, f"Expected kernel '{name}' in test-kernel.cpp"


# =============================================================================
# Acoustic smoke — catches equation-specific regressions
# =============================================================================


class TestAcousticSmoke:
    """Acoustic has fewer quantities (4 vs 9) — a separate pass verifies
    no equation-specific path is ELBOW-DEPENDENT on numberOfQuantities=9."""

    def test_acoustic_generates(self, tmp_path):
        result = _invoke_generate(tmp_path, equation="acoustic", order=3)
        assert result.returncode == 0, (
            f"generate.py acoustic failed:\n"
            f"stdout:\n{result.stdout[-1000:]}\n"
            f"stderr:\n{result.stderr[-1000:]}"
        )
        assert (tmp_path / "equation-acoustic-3-double").is_dir()


# =============================================================================
# Poroelastic smoke — catches equation-specific regressions
# =============================================================================


class TestPoroelasticSmoke:
    """Poroelastic has more quantities (13 vs 9) — a separate pass verifies
    no equation-specific path is ELBOW-DEPENDENT on numberOfQuantities=13."""

    def test_acoustic_generates(self, tmp_path):
        result = _invoke_generate(tmp_path, equation="poroelastic", order=3)
        assert result.returncode == 0, (
            f"generate.py poroelastic failed:\n"
            f"stdout:\n{result.stdout[-1000:]}\n"
            f"stderr:\n{result.stderr[-1000:]}"
        )
        assert (tmp_path / "equation-poroelastic-3-double").is_dir()


# =============================================================================
# Single precision — catches f32-specific alignment regressions
# =============================================================================


class TestSinglePrecisionSmoke:

    def test_elastic_single_precision_generates(self, tmp_path):
        result = _invoke_generate(tmp_path, equation="elastic", order=3, precision="s")
        assert result.returncode == 0, result.stderr[-1000:]
        # Subfolder name uses "single" not "s"
        assert (tmp_path / "equation-elastic-3-single").is_dir()

    def test_subfolder_name_for_f32(self, tmp_path):
        """generate.py maps --precision=f32 to subfolder suffix 'single'."""
        result = _invoke_generate(
            tmp_path, equation="elastic", order=3, precision="f32"
        )
        assert result.returncode == 0
        assert (tmp_path / "equation-elastic-3-single").is_dir()


# =============================================================================
# Fused simulations — catches multi-sim layout regressions
# =============================================================================


class TestFusedSimsSmoke:

    def test_fused_8_sims_generates(self, tmp_path):
        result = _invoke_generate(tmp_path, equation="elastic", order=3, multi_sims=8)
        assert result.returncode == 0, result.stderr[-1000:]
        # Subfolder name suffixes the multi-sim count with -f<N>
        assert (tmp_path / "equation-elastic-3-double-f8").is_dir()


# =============================================================================
# Error paths
# =============================================================================


class TestGenerateErrorPaths:

    def test_unknown_equation_fails_cleanly(self, tmp_path):
        """A typo in --equations should fail with a clear RuntimeError,
        not a silent miscompile."""
        result = _invoke_generate(tmp_path, equation="fantasy_equation", order=3)
        assert result.returncode != 0
        # Message from generate.py's own raise
        assert (
            "Could not find kernels for fantasy_equation" in result.stderr
            or "fantasy_equation" in result.stderr
        )

    def test_unknown_gemm_tool_fails(self, tmp_path):
        """--gemm_tools=fakegemm triggers the 'Unknown GEMM tool' path."""
        result = subprocess.run(
            [
                sys.executable,
                str(GENERATE),
                "--equations",
                "elastic",
                "--matricesDir",
                str(CODEGEN_DIR / "matrices"),
                "--outputDir",
                str(tmp_path),
                "--host_arch",
                "hsw",
                "--order",
                "3",
                "--precision",
                "d",
                "--numberOfMechanisms",
                "0",
                "--memLayout",
                "auto",
                "--multipleSimulations",
                "1",
                "--PlasticityMethod",
                "ip",
                "--gemm_tools",
                "fakegemm",  # <-- bogus
                "--drQuadRule",
                "dunavant",
                "--device_backend",
                "none",
            ],
            cwd=str(CODEGEN_DIR),
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode != 0
        # Error message from generate.py explicitly
        combined = result.stdout + result.stderr
        assert "Unknown GEMM tool" in combined or "fakegemm" in combined
