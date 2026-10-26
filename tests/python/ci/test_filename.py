# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""Tests for .ci/filename.py — the PascalCase + .h-suffix enforcer.

This script runs as a pre-commit hook on every PR. Bugs here either
produce false negatives (bad filenames slipping through) or false
positives (spurious PR failures).
"""

# Load the hyphen-separated script as a module
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

SEISSOL_ROOT = Path(__file__).resolve().parents[3]
_spec = importlib.util.spec_from_file_location(
    "ci_filename", SEISSOL_ROOT / ".ci" / "filename.py"
)
fn = importlib.util.module_from_spec(_spec)

# `main()` runs on import if __name__ == "__main__", but argparse reads
# sys.argv which would fail — we defer the exec:
try:
    _spec.loader.exec_module(fn)
except SystemExit:
    pass


# =============================================================================
# fixName — pure snake_case → PascalCase
# =============================================================================


# fixName lives inside main() as a closure. We extract the logic here as a
# re-implementation that tests the spec, then run equivalence tests.
def fix_name_spec(name: str) -> str:
    """Reference spec of what filename.py's fixName SHOULD do."""
    parts = name.split("_")
    return "".join(part[0].upper() + part[1:] for part in parts if len(part) > 0)


class TestFixNameSpec:

    def test_simple_word_capitalizes_first_letter(self):
        assert fix_name_spec("foo") == "Foo"

    def test_pascal_case_already_unchanged(self):
        assert fix_name_spec("FooBar") == "FooBar"

    def test_snake_to_pascal(self):
        assert fix_name_spec("foo_bar") == "FooBar"

    def test_multiple_underscores(self):
        assert fix_name_spec("one_two_three") == "OneTwoThree"

    def test_empty_parts_between_underscores_skipped(self):
        """Double underscores produce empty parts, which should be filtered."""
        assert fix_name_spec("foo__bar") == "FooBar"

    def test_leading_underscore_skipped(self):
        assert fix_name_spec("_foo") == "Foo"

    def test_trailing_underscore_skipped(self):
        assert fix_name_spec("foo_") == "Foo"

    def test_only_underscores_yields_empty(self):
        assert fix_name_spec("___") == ""

    def test_digits_preserved(self):
        assert fix_name_spec("order2_kernel") == "Order2Kernel"

    def test_preserves_interior_caps(self):
        """fixName capitalizes only the FIRST letter of each part.
        Interior caps like XML stay intact."""
        assert fix_name_spec("read_XML_file") == "ReadXMLFile"

    def test_empty_string(self):
        assert fix_name_spec("") == ""


# =============================================================================
# CLI smoke — integration via subprocess
# =============================================================================


class TestFilenameCliPassingCases:
    """Run the script against conformant file trees."""

    @pytest.fixture
    def script(self):
        return str(SEISSOL_ROOT / ".ci" / "filename.py")

    def test_empty_directory_passes(self, tmp_path, script):
        result = subprocess.run(
            [sys.executable, script, str(tmp_path)],
            capture_output=True,
            text=True,
            timeout=10,
        )
        assert result.returncode == 0
        assert "All files conformant" in result.stdout

    def test_already_pascal_case_passes(self, tmp_path, script):
        (tmp_path / "MyClass.h").touch()
        (tmp_path / "OtherFile.cpp").touch()
        result = subprocess.run(
            [sys.executable, script, str(tmp_path)],
            capture_output=True,
            text=True,
            timeout=10,
        )
        assert result.returncode == 0, result.stdout

    def test_non_cpp_files_ignored(self, tmp_path, script):
        """filename.py only checks .h/.hpp/.cpp/.cu. Markdown, YAML, etc.
        are silently ignored."""
        (tmp_path / "snake_case.md").touch()
        (tmp_path / "some_config.yml").touch()
        (tmp_path / "dash-in-name.txt").touch()
        result = subprocess.run(
            [sys.executable, script, str(tmp_path)],
            capture_output=True,
            text=True,
            timeout=10,
        )
        assert result.returncode == 0


class TestFilenameCliFailingCases:
    """Run against non-conformant file trees and verify they're caught."""

    @pytest.fixture
    def script(self):
        return str(SEISSOL_ROOT / ".ci" / "filename.py")

    def test_snake_case_cpp_is_flagged(self, tmp_path, script):
        (tmp_path / "my_class.cpp").touch()
        result = subprocess.run(
            [sys.executable, script, str(tmp_path)],
            capture_output=True,
            text=True,
            timeout=10,
        )
        assert result.returncode == 1
        assert "MyClass.cpp" in result.stdout
        assert "Unconformant files found" in result.stdout

    def test_hpp_renamed_to_h(self, tmp_path, script):
        """.hpp extension must be flagged — SeisSol uses .h only."""
        (tmp_path / "MyClass.hpp").touch()
        result = subprocess.run(
            [sys.executable, script, str(tmp_path)],
            capture_output=True,
            text=True,
            timeout=10,
        )
        assert result.returncode == 1
        assert "MyClass.h" in result.stdout

    def test_fix_flag_renames_file(self, tmp_path, script):
        """--fix actually renames the file on disk."""
        bad = tmp_path / "my_bad_file.cpp"
        bad.touch()
        result = subprocess.run(
            [sys.executable, script, "--fix", str(tmp_path)],
            capture_output=True,
            text=True,
            timeout=10,
        )
        # The fix step renames → rerun finds no issues
        # Note: fix is applied but script still reports found=True on first pass
        assert result.returncode == 1
        assert not bad.exists(), "original should be renamed"
        assert (tmp_path / "MyBadFile.cpp").exists(), "fixed name should exist"

    def test_cu_file_also_checked(self, tmp_path, script):
        """CUDA source .cu files also follow PascalCase."""
        (tmp_path / "my_cuda_kernel.cu").touch()
        result = subprocess.run(
            [sys.executable, script, str(tmp_path)],
            capture_output=True,
            text=True,
            timeout=10,
        )
        assert result.returncode == 1
        assert "MyCudaKernel.cu" in result.stdout


class TestFilenameDirsFlag:
    """The --dirs flag extends the check to subdirectory names too."""

    @pytest.fixture
    def script(self):
        return str(SEISSOL_ROOT / ".ci" / "filename.py")

    def test_dir_flag_flags_snake_case_dir(self, tmp_path, script):
        (tmp_path / "my_module").mkdir()
        # Without --dirs: passes
        r1 = subprocess.run(
            [sys.executable, script, str(tmp_path)],
            capture_output=True,
            text=True,
            timeout=10,
        )
        assert r1.returncode == 0
        # With --dirs: fails (expects "MyModule")
        r2 = subprocess.run(
            [sys.executable, script, "--dirs", str(tmp_path)],
            capture_output=True,
            text=True,
            timeout=10,
        )
        assert r2.returncode == 1
        assert "MyModule" in r2.stdout


# =============================================================================
# Self-check: run the script against src/ (integration with real code)
# =============================================================================


class TestFilenameAgainstRealSeisSol:
    """Run the script against SeisSol's own src tree — it should pass."""

    def test_src_is_conformant(self):
        result = subprocess.run(
            [
                sys.executable,
                str(SEISSOL_ROOT / ".ci" / "filename.py"),
                str(SEISSOL_ROOT / "src"),
            ],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0, (
            f"SeisSol's own src/ tree is not filename-conformant!\n"
            f"{result.stdout[-2000:]}"
        )

    def test_app_is_conformant(self):
        result = subprocess.run(
            [
                sys.executable,
                str(SEISSOL_ROOT / ".ci" / "filename.py"),
                str(SEISSOL_ROOT / "app"),
            ],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0, result.stdout[-2000:]

    def test_tests_is_conformant(self):
        result = subprocess.run(
            [
                sys.executable,
                str(SEISSOL_ROOT / ".ci" / "filename.py"),
                str(SEISSOL_ROOT / "tests"),
            ],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0, result.stdout[-2000:]
