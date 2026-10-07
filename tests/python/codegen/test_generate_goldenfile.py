# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""End-to-end smoke + golden-file tests for generate.py.

These invoke the REAL generator with minimal configurations and verify:
 - All the headers at the top level are produced
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
import json
import os
import re
import subprocess
import sys
from pathlib import Path

# Resolve codegen dir via the sys.path set up in conftest.py
import kernels.memlayout as _ml
import pytest

CODEGEN_DIR = Path(_ml.__file__).resolve().parent.parent
GENERATE = CODEGEN_DIR / "generate.py"


def _invoke_generate(
    outdir,
    equation="elastic",
    order=3,
    precision="d",
    multi_sims=1,
    mechanisms=0,
    mode="codegen",
    solver=None,
    target=None,
):
    """Run generate.py with the given config. Returns CompletedProcess.

    The hash seed is pinned because the generator iterates over sets in a few
    places, which permutes the order temporaries are declared in. Without it
    two runs of the same input differ, and no diff of the generated code means
    anything.
    """
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
            "--numMechanisms",
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
            "--mode",
            mode,
            *(["--solver", solver] if solver is not None else []),
            *(["--codegen_target", target] if target is not None else []),
        ],
        env={**os.environ, "PYTHONHASHSEED": "0"},
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
    """The headers at the top level: those of the metagen, which name the code
    of the equation by key, and generate.py's `forward_files` list."""

    EXPECTED_TOP_LEVEL_HEADERS = {
        "init.h",
        "kernel.h",
        "pool.h",
        "quantities.h",
        "tensor.h",
    }

    def test_top_level_headers_produced(self, generated_elastic_o3):
        outdir, _ = generated_elastic_o3
        produced = {p.name for p in outdir.iterdir() if p.is_file()}
        missing = self.EXPECTED_TOP_LEVEL_HEADERS - produced
        assert not missing, f"Missing top-level headers: {missing}"

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

    def test_quantities_h_uses_iwyu_pragma(self, generated_elastic_o3):
        outdir, _ = generated_elastic_o3
        content = (outdir / "quantities.h").read_text()
        assert "IWYU pragma: begin_exports" in content
        assert "IWYU pragma: end_exports" in content

    def test_top_level_headers_include_the_equation_only(self, generated_elastic_o3):
        """init.h includes the equation's init.h; the code of general/ belongs
        to no configuration and is included from there."""
        outdir, _ = generated_elastic_o3
        content = (outdir / "init.h").read_text()
        assert "equation-elastic-3-double/init.h" in content
        assert "general/init.h" not in content

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
# runtime.h — the kernels reached by the variant of their configuration
# =============================================================================


class TestRuntime:
    """The equation is generated through a yateto metagen, which adds runtime.h:
    kernels reached by the variant of the configuration they compute for, with
    operands as views. CMake builds its units from what `--mode collect` lists.
    """

    EQUATION = "equation-elastic-3-double"

    def test_runtime_files_produced(self, generated_elastic_o3):
        outdir, _ = generated_elastic_o3
        for name in ["runtime.h", "runtime.cpp", "variant.h"]:
            assert (outdir / name).is_file(), f"Expected {name} at the top level"
        # the unit that binds views to the kernels of the equation
        assert (outdir / self.EQUATION / "runtime.cpp").is_file()

    def test_variant_h_keys_the_configuration(self, generated_elastic_o3):
        """The id of a configuration is its variant, runtime::variantOf<Config0>()."""
        outdir, _ = generated_elastic_o3
        content = (outdir / "variant.h").read_text()
        assert '#include "Config.h"' in content
        assert "VariantOf<seissol::Config0>" in content

    def test_code_is_named_by_the_key_of_its_configuration(self, generated_elastic_o3):
        """The code of the equation is in a namespace of its own, and the headers
        at the top level name it by the key of its configuration:
        seissol::kernel::X<seissol::Config0>, seissol::Pool<seissol::Config0>."""
        outdir, _ = generated_elastic_o3
        space = "yatetometagen_" + self.EQUATION.replace("-", "_")
        content = (outdir / self.EQUATION / "kernel.h").read_text()
        assert f"namespace seissol {{\n  namespace {space} {{" in content
        for name, prefix in [
            ("init", "init::"),
            ("kernel", "kernel::"),
            ("tensor", "tensor::"),
        ]:
            typed = (outdir / f"{name}.h").read_text()
            assert f'#include "{self.EQUATION}/{name}.h"' in typed
            assert f"using Type = ::seissol::{space}::{prefix}" in typed
        pool = (outdir / "pool.h").read_text()
        assert (
            f"struct Internal_Pool<seissol::Config0> {{ using Type = ::seissol::{space}::Pool; }};"
            in pool
        )

    def test_optional_tensors_are_named_in_every_configuration(
        self, generated_elastic_o3
    ):
        """Qane is no tensor of the elastic equation, but its name exists, as
        `void`, which kernels::size counts as empty."""
        outdir, _ = generated_elastic_o3
        tensor = (outdir / "tensor.h").read_text()
        assert "template<typename Arg0> using Qane = " in tensor
        assert "Internal_Qane<seissol::Config0>" not in tensor

    def test_collect_lists_what_codegen_writes(self, generated_elastic_o3, tmp_path):
        outdir, _ = generated_elastic_o3
        result = _invoke_generate(tmp_path, mode="collect")
        assert result.returncode == 0, result.stderr[-1000:]
        targets = json.loads((tmp_path / "targets.json").read_text())

        assert targets["runtime"]["kernels"] == ["runtime.cpp"]
        assert sorted(targets["runtime"]["headers"]) == [
            "init.h",
            "kernel.h",
            "pool.h",
            "runtime.h",
            "tensor.h",
            "variant.h",
        ]
        assert f"{self.EQUATION}/runtime.cpp" in targets[self.EQUATION]["kernels"]

        listed = [
            path
            for target in targets.values()
            for kind in ("kernels", "device", "tests", "headers")
            for path in target.get(kind, [])
        ]
        missing = [path for path in listed if not (outdir / path).is_file()]
        assert (
            not missing
        ), f"collect lists files that codegen does not write: {missing}"

    def test_steps_list_what_they_write(self, generated_elastic_o3, tmp_path):
        """A build runs every step of the code generation as a command of its
        own, and has to know what each one writes: every file that codegen
        writes is an output of exactly one step (see STEP_ALL in
        generate.py). alignment.h is written when CMake runs as well."""
        outdir, _ = generated_elastic_o3
        result = _invoke_generate(tmp_path, mode="collect")
        assert result.returncode == 0, result.stderr[-1000:]
        steps = json.loads((tmp_path / "steps.json").read_text())

        listed = [path for step in steps for path in step["outputs"]]
        twice = {path for path in listed if listed.count(path) > 1}
        assert not twice, f"files listed by more than one step: {twice}"
        written = {
            str(path.relative_to(outdir))
            for path in outdir.rglob("*")
            if path.is_file()
        } - {"alignment.h"}
        assert set(listed) == written, (
            f"written, but listed by no step: {written - set(listed)}; "
            f"listed, but not written: {set(listed) - written}"
        )

    def test_steps_write_what_one_run_writes(self, generated_elastic_o3, tmp_path):
        """Run one after the other, as a build may, the steps write the same
        code as one run that takes all of them."""
        outdir, _ = generated_elastic_o3
        result = _invoke_generate(tmp_path, mode="collect")
        assert result.returncode == 0, result.stderr[-1000:]
        steps = json.loads((tmp_path / "steps.json").read_text())
        for step in steps:
            result = _invoke_generate(tmp_path, target=step["target"])
            assert result.returncode == 0, result.stderr[-1000:]
        for step in steps:
            for path in step["outputs"]:
                assert (tmp_path / path).read_bytes() == (
                    outdir / path
                ).read_bytes(), f"{path} differs"

    @staticmethod
    def _runtime_kernels(outdir):
        """The kernels runtime.h declares, in any namespace."""
        content = (outdir / "runtime.h").read_text()
        blocks = re.findall(
            r"namespace kernel \{(.*?)\} // namespace kernel", content, re.S
        )
        return {
            name
            for block in blocks
            for name in re.findall(r"^\s*struct (\w+) \{", block, re.M)
        }

    def test_kernels_of_setup_and_output_take_views(self, generated_elastic_o3):
        """The kernels only setup and output run are reached through runtime.h,
        the ones of the time step are not: they keep their operands as pointers
        (see kernels.common.cold_kernel_attrs)."""
        outdir, _ = generated_elastic_o3
        kernels = self._runtime_kernels(outdir)
        for name in [
            "computeFluxSolverLocal",
            "foldDirichlet",
            "projectIniCond",
            "transformNRF",
            "evaluateDOFSAtPoint",
            "evalAtQP",
            "momentQQCompute",
            "plProject",
            "projectNodalToVtkFace",
            "rotateFluxMatrix",
            "evaluateFaceAlignedDOFSAtPoint",
            "accumulateStaticFrictionalWork",
        ]:
            assert name in kernels, f"Expected kernel '{name}' in runtime.h"
        for name in [
            "volume",
            "localFlux",
            "neighboringFlux",
            "derivative",
            "evaluateAndRotateQAtInterpolationPoints",
        ]:
            assert name not in kernels, f"Kernel '{name}' of the time step in runtime.h"

    def test_anelastic_moments_take_views(self, tmp_path):
        """The energy output of the viscoelastic material runs the moments of
        the anelastic unknowns through runtime.h."""
        result = _invoke_generate(
            tmp_path, equation="viscoelastic", mechanisms=3, solver="linearckanelastic"
        )
        assert result.returncode == 0, result.stderr[-1000:]
        kernels = self._runtime_kernels(tmp_path)
        assert {"momentQaneQaneCompute", "momentQQaneCompute"} <= kernels


# =============================================================================
# Acoustic smoke — catches equation-specific regressions
# =============================================================================


class TestAcousticSmoke:
    """Acoustic has fewer quantities (4 vs 9) — a separate pass verifies
    no equation-specific path is ELBOW-DEPENDENT on numQuantities=9."""

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
    no equation-specific path is ELBOW-DEPENDENT on numQuantities=13."""

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
                "--numMechanisms",
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


class TestReproducibility:
    """The generated code has to be a function of the inputs alone.

    Every refactor of the generator in this tree has leaned on the same check:
    generate before, generate after, diff. That is only evidence if two runs of
    the *same* input agree -- and they do not by default, because the generator
    iterates over sets and the declaration order of temporaries follows the
    hash seed.
    """

    @staticmethod
    def _snapshot(root):
        return {
            path.relative_to(root).as_posix(): path.read_bytes()
            for path in sorted(root.rglob("*"))
            if path.is_file()
        }

    def test_two_runs_of_the_same_input_agree(self, tmp_path):
        first, second = tmp_path / "first", tmp_path / "second"
        for outdir in (first, second):
            outdir.mkdir()
            result = _invoke_generate(outdir)
            assert result.returncode == 0, result.stderr

        left, right = self._snapshot(first), self._snapshot(second)
        assert sorted(left) == sorted(right), "the two runs produced different files"
        differing = [name for name in left if left[name] != right[name]]
        assert not differing, (
            f"generated code is not reproducible; {len(differing)} file(s) differ, "
            f"e.g. {differing[:3]}"
        )
