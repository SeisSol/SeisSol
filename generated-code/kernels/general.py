# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

from yateto import Tensor
from yateto.input import parseJSONMatrixFile


def addStiffnessTensor(generator):
    stiffnessTensor = Tensor("stiffnessTensor", (3, 3, 3, 3))
    direction = Tensor("direction", (3,))
    christoffel = Tensor("christoffel", (3, 3))

    computeChristoffel = (
        christoffel["ik"] <= stiffnessTensor["ijkl"] * direction["j"] * direction["l"]
    )
    generator.add("computeChristoffel", computeChristoffel)


def includeMatrices(matricesDir):
    matrices = parseJSONMatrixFile("{}/sampling_directions.json".format(matricesDir))
    return set([matrices.samplingDirections])
