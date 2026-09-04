# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""CLI smoke tests for codegen/generate.py

generate.py is invoked from CMake with a large set of build-parameter-dependent
arguments. An argument-schema regression there is expensive because it only
surfaces late in the build pipeline. These tests isolate the CLI surface.

Scope:
 - argparse validation (choices enforcement on --precision)
 - exit-code contract
 - pure helper functions (generate_kernel_name_prefix)

Out of scope:
 - Running full generation — that requires yateto + matrices + a target
   compiler architecture. That belongs in a golden-file test (see plan).
"""

import subprocess
import sys
from pathlib import Path

import pytest

SEISSOL_ROOT = Path(__file__).resolve().parents[3]
GENERATE = SEISSOL_ROOT / "codegen" / "generate.py"


def run_generate(*args):
    """Run generate.py with the given arguments; return CompletedProcess."""
    return subprocess.run(
        [sys.executable, str(GENERATE), *args],
        capture_output=True,
        text=True,
        cwd=str(SEISSOL_ROOT / "codegen"),  # needs relative `kernels.` imports
        timeout=30,
    )


# =============================================================================
# argparse surface
# =============================================================================


class TestGenerateArgparse:

    def test_help_exits_zero(self):
        result = run_generate("--help")
        assert result.returncode == 0, result.stderr

    def test_help_lists_all_documented_flags(self):
        result = run_generate("--help")
        out = result.stdout
        # All flags currently exposed by the argparse schema
        for flag in [
            "--equations",
            "--matricesDir",
            "--outputDir",
            "--host_arch",
            "--device_backend",
            "--device_arch",
            "--device_vendor",
            "--order",
            "--precision",
            "--numberOfMechanisms",
            "--vectorsize",
            "--memLayout",
            "--multipleSimulations",
            "--PlasticityMethod",
            "--gemm_tools",
            "--device_codegen",
            "--drQuadRule",
            "--enable_premultiply_flux",
            "--disable_premultiply_flux",
            "--executable_libxsmm",
            "--executable_pspamm",
        ]:
            assert flag in out, f"Flag {flag} missing from --help"

    def test_invalid_precision_is_rejected(self):
        result = run_generate("--precision", "float")
        # argparse writes to stderr and exits 2 for invalid choice
        assert result.returncode == 2
        assert "invalid choice" in result.stderr
        assert "'float'" in result.stderr

    @pytest.mark.parametrize("precision", ["s", "d", "f32", "f64"])
    def test_valid_precision_choices_accepted(self, precision):
        # --help still exits 0 with any valid choice
        # This tests that the schema hasn't been narrowed accidentally
        result = run_generate("--precision", precision, "--help")
        assert result.returncode == 0

    def test_non_integer_order_rejected(self):
        result = run_generate("--order", "abc")
        # type=int rejection → argparse exits 2
        assert result.returncode == 2
        assert "invalid int value" in result.stderr

    def test_non_integer_vectorsize_rejected(self):
        result = run_generate("--vectorsize", "quarter")
        assert result.returncode == 2

    def test_mutually_inverse_flux_flags(self):
        # --enable_premultiply_flux and --disable_premultiply_flux are paired
        # via store_true/store_false on the same dest. The LAST one wins.
        # We test just that neither crashes the parser.
        r1 = run_generate("--enable_premultiply_flux", "--help")
        r2 = run_generate("--disable_premultiply_flux", "--help")
        assert r1.returncode == 0
        assert r2.returncode == 0


# =============================================================================
# Running without required args — surfaces the "unprotected None" pattern
# =============================================================================


class TestGenerateWithoutArgs:
    """None of the arguments are actually marked required=True in argparse,
    yet downstream code dereferences them (e.g. cmdLineArgs.vectorsize == 0,
    cmdLineArgs.gemm_tools.replace(" ", "")). Running with no args surfaces
    this: either you get a TypeError/AttributeError, or the script silently
    produces nonsense.

    Documents current behavior without prescribing a fix.
    """

    def test_missing_args_fails_cleanly_or_with_traceback(self):
        # Expected: process crashes, but with a recognizable failure mode
        # NOT: silent success and bogus output
        result = run_generate()
        assert result.returncode != 0


# =============================================================================
# Pure helpers that are importable
# =============================================================================


class TestKernelNamePrefix:
    """generate_kernel_name_prefix is called all over the codegen."""

    def test_gpu_prefix(self):
        # Lazy import — requires yateto; skip if unavailable
        try:
            import sys as _s

            _s.path.insert(0, str(SEISSOL_ROOT / "codegen"))
            from kernels.common import generate_kernel_name_prefix
        except ImportError:
            pytest.skip("yateto submodule not initialized")
        assert generate_kernel_name_prefix("gpu") == "gpu_"

    def test_cpu_no_prefix(self):
        try:
            import sys as _s

            _s.path.insert(0, str(SEISSOL_ROOT / "codegen"))
            from kernels.common import generate_kernel_name_prefix
        except ImportError:
            pytest.skip("yateto submodule not initialized")
        assert generate_kernel_name_prefix("cpu") == ""

    def test_other_targets_return_empty(self):
        """Any non-'gpu' value should return empty — this is the invariant
        callers rely on when they concatenate the prefix."""
        try:
            import sys as _s

            _s.path.insert(0, str(SEISSOL_ROOT / "codegen"))
            from kernels.common import generate_kernel_name_prefix
        except ImportError:
            pytest.skip("yateto submodule not initialized")
        assert generate_kernel_name_prefix("cpu") == ""
        assert generate_kernel_name_prefix("anything_else") == ""
        assert generate_kernel_name_prefix("") == ""


# =============================================================================
# Platform derivation — inline in generate.py, currently untestable without
# refactoring. Documents the logic for future extraction.
# =============================================================================


class TestPlatformDerivation:
    """generate.py contains this logic inline:
        gpu_platforms = ["cuda", "hip", "hipsycl", "acpp", "oneapi"]
        targets = ["gpu", "cpu"] if device_backend in gpu_platforms else ["cpu"]

    If this logic were extracted to a function `derive_targets(device_backend)`,
    it would be trivially testable. Recommended refactor.
    """

    @pytest.mark.parametrize(
        "backend,expected",
        [
            ("cuda", ["gpu", "cpu"]),
            ("hip", ["gpu", "cpu"]),
            ("hipsycl", ["gpu", "cpu"]),
            ("acpp", ["gpu", "cpu"]),
            ("oneapi", ["gpu", "cpu"]),
            ("none", ["cpu"]),
            (None, ["cpu"]),
            ("", ["cpu"]),
            ("fantasy", ["cpu"]),
        ],
    )
    def test_platform_derivation_oracle(self, backend, expected):
        """Oracle for the derive-targets logic.

        If someone extracts this logic, drop the import-the-function call
        in place of the re-implementation here — the test stays valid.
        """
        gpu_platforms = ["cuda", "hip", "hipsycl", "acpp", "oneapi"]
        targets = ["gpu", "cpu"] if backend in gpu_platforms else ["cpu"]
        assert targets == expected
