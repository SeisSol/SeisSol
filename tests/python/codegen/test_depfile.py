# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The inputs that a step of the code generator writes into its depfile
(kernels.depfile)."""

import json
import os

import kernels  # noqa: F401 — establish sys.path via conftest
from kernels.depfile import InputRecorder, imported_files, write_depfile


def test_the_recorder_keeps_the_files_read_under_its_roots(tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    read = root / "matrix.json"
    read.write_text("{}")
    outside = tmp_path / "outside.json"
    outside.write_text("{}")

    recorder = InputRecorder([root])
    with open(read) as file:
        json.load(file)
    with open(outside) as file:
        json.load(file)
    # a file written is an output, not an input
    with open(root / "written.json", "w") as file:
        file.write("{}")
    # nor is a compiled module, which the imports write themselves
    compiled = root / "__pycache__" / "module.cpython.pyc"
    compiled.parent.mkdir()
    compiled.write_bytes(b"")
    with open(compiled, "rb") as file:
        file.read()

    assert recorder.files == {str(read)}


def test_the_imported_modules_leave_out_the_standard_library():
    files = imported_files()
    assert os.path.abspath(kernels.depfile.__file__) in files
    assert os.path.abspath(json.__file__) not in files


def test_the_depfile_names_the_first_output_and_escapes_the_paths(tmp_path):
    depfile = tmp_path / "step.d"
    write_depfile(
        str(depfile), str(tmp_path / "first.h"), {"/b/with space.py", "/a/plain.py"}
    )
    target, dependencies = depfile.read_text().split(":", 1)
    assert target == str(tmp_path / "first.h")
    assert dependencies == " \\\n  /a/plain.py \\\n  /b/with\\ space.py\n"
