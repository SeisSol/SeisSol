# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The basis of the discretisation, and where its matrices lie.

The matrix files of the Dubiner basis lie directly in the matrix directory
(codegen/matrices), next to the files that do not depend on the basis: the
equations, the time basis of the space-time predictor, the sampling directions
and the points of the visualisation. Every other basis has a directory of its
own in there, with the files that depend on it and a basis.json that states
what may be relied on, both written by SeisSol/matrices (generate.py seissol).
"""

import json
import os

import numpy as np

#: The bases a build can have; the first one is the default.
NAMES = ("dubiner", "dubiner-normalized", "bernstein")

#: What the Dubiner basis states of itself, in the form of basis.json.
_DUBINER = {
    "family": "dubiner",
    "hierarchical": True,
    "orthogonal": True,
    "orthonormal": False,
    "mass": "diagonal",
}


def count(degree, dim=3):
    """How many functions span the polynomials up to `degree` in `dim`
    variables."""
    result = 1
    for k in range(1, dim + 1):
        result = result * (degree + k) // k
    return result


def dense(tensor, transposed=False):
    """The values of a matrix as a dense array, in the orientation of the file it
    was read from where it was read `transposed`."""
    matrix = np.zeros(tensor.shape())
    for index, value in tensor.values().items():
        matrix[index] = float(value)
    return matrix.T if transposed else matrix


def directory(matricesDir, name):
    """The directory of the matrix files of a basis."""
    if name not in NAMES:
        raise ValueError(f"unknown basis {name}, expected one of {', '.join(NAMES)}")
    return matricesDir if name == NAMES[0] else os.path.join(matricesDir, name)


def shared(matricesDir, filename):
    """The path of a file that does not depend on the basis, given the directory
    of the basis: the directory of another basis lies in the one of the Dubiner
    basis, which holds those files for all of them."""
    if os.path.exists(os.path.join(matricesDir, "basis.json")):
        matricesDir = os.path.dirname(os.path.normpath(matricesDir))
    return os.path.join(matricesDir, filename)


class Basis:
    """What the generator may rely on of the basis whose matrices lie in a
    directory.

    hierarchical: the functions up to a degree come first and span the
    polynomials of that degree, so that a derivative, which lowers the degree,
    leaves the trailing ones zero. orthogonal: the mass matrices are diagonal.
    Where a basis is both, the projection onto a lower degree keeps the leading
    modes."""

    def __init__(self, matricesDir):
        path = os.path.join(matricesDir, "basis.json")
        if os.path.exists(path):
            with open(path) as file:
                description = json.load(file)
        else:
            description = _DUBINER
        self.name = description["family"]
        self.hierarchical = bool(description["hierarchical"])
        self.orthogonal = bool(description["orthogonal"])
        self.orthonormal = bool(description["orthonormal"])
        self.mass = description["mass"]
