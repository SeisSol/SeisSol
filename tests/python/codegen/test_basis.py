# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Unit tests for the basis of a build and where its matrices lie
(codegen/kernels/basis.py)."""

import json
import os

import pytest
from kernels.basis import NAMES, Basis, count, directory, shared


class TestDirectory:
    def test_dubiner_lies_in_the_matrix_directory(self):
        assert directory("matrices", "dubiner") == "matrices"

    @pytest.mark.parametrize("name", NAMES[1:])
    def test_another_basis_lies_in_a_directory_of_its_name(self, name):
        assert directory("matrices", name) == os.path.join("matrices", name)

    def test_unknown_basis_is_rejected(self):
        with pytest.raises(ValueError):
            directory("matrices", "legendre")


class TestShared:
    def test_shared_files_of_the_dubiner_basis_lie_next_to_its_matrices(self, tmp_path):
        assert shared(str(tmp_path), "vtkbase.json") == str(tmp_path / "vtkbase.json")

    def test_shared_files_of_another_basis_lie_one_directory_up(self, tmp_path):
        own = tmp_path / "bernstein"
        own.mkdir()
        (own / "basis.json").write_text("{}")
        assert shared(str(own), "vtkbase.json") == str(tmp_path / "vtkbase.json")


class TestBasis:
    def test_without_a_description_it_is_the_dubiner_basis(self, tmp_path):
        basis = Basis(str(tmp_path))
        assert basis.name == "dubiner"
        assert basis.hierarchical and basis.orthogonal and not basis.orthonormal
        assert basis.mass == "diagonal"

    def test_reads_the_description(self, tmp_path):
        description = {
            "family": "bernstein",
            "hierarchical": False,
            "orthogonal": False,
            "orthonormal": False,
            "mass": "dense",
        }
        (tmp_path / "basis.json").write_text(json.dumps(description))
        basis = Basis(str(tmp_path))
        assert basis.name == "bernstein"
        assert not basis.hierarchical and not basis.orthogonal
        assert basis.mass == "dense"


def test_count():
    # the polynomials up to a degree in three variables, and in two
    assert [count(degree) for degree in range(-1, 5)] == [0, 1, 4, 10, 20, 35]
    assert [count(degree, 2) for degree in range(4)] == [1, 3, 6, 10]
